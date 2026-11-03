#!/usr/bin/env python3
"""
run_md.py - Automated protein-ligand molecular dynamics with OpenMM
            (Apple Silicon / OpenCL, single precision)

USAGE
-----
    1. Put `protein.pdb` and `ligand.sdf` in the folder you run this from.
    2. python run_md.py

INPUT REQUIREMENTS (the script checks what it can, but you must supply sane inputs)
    protein.pdb : receptor ONLY (no ligand, crystal waters or ions are stripped).
                  Should be a continuous chain; unmodelled internal gaps will make the
                  force field fail with a "No template found" error.
    ligand.sdf  : ONE molecule, 3D coordinates, EXPLICIT hydrogens, already in its
                  docked pose in the SAME coordinate frame as protein.pdb.
                  (An SDF straight out of a docking run that kept hydrogens is fine.)

ENVIRONMENT (conda-forge, osx-arm64)
    conda create -n openmm_md -c conda-forge python=3.11 \
        openmm openmmforcefields openff-toolkit ambertools pdbfixer rdkit numpy
    conda activate openmm_md
    python -m openmm.testInstallation      # should list 'OpenCL'

OUTPUTS
    complex_production.dcd   production trajectory (frame every DCD_INTERVAL_PS)
    system_solvated.pdb      topology/atom-order reference for the DCD (load both together)
    simulation.log           energy / temperature / volume / density / speed, all phases
    production.chk           checkpoint (platform-specific; restart on the same machine)
    final_state.xml          final positions + velocities
    final_frame.pdb          final frame (periodic box enforced)

PROTOCOL
    minimize -> NVT (heavy-atom restraints) -> NPT (heavy-atom restraints)
             -> unrestrained NPT production
"""

import os
import sys
import time

# =============================================================================
# USER SETTINGS - edit these
# =============================================================================
PROTEIN_FILE = "sEH_clean.pdb"
LIGAND_FILE = "F01_ligand.sdf"

PRODUCTION_NS = 1.0          # <-- production length in nanoseconds (edit me)
NVT_PS = 100.0               # NVT equilibration length (ps)
NPT_PS = 100.0               # NPT equilibration length (ps)

TEMPERATURE_K = 300.0
PRESSURE_ATM = 1.0
TIMESTEP_FS = 2.0            # 2 fs is safe with HBond constraints + rigid water
FRICTION_PER_PS = 1.0        # Langevin collision frequency

BOX_PADDING_NM = 1.0         # solute-to-box-edge padding
BOX_SHAPE = "cube"           # 'cube' | 'dodecahedron' | 'octahedron'
IONIC_STRENGTH_M = 0.15      # NaCl concentration (molar), system also neutralized
PH = 7.0                     # protonation state of the protein

# Positional restraints on protein + ligand heavy atoms during minimization/NVT/NPT.
# They are removed (k -> 0) for production. Prevents the pose from being
# distorted while the solvent relaxes. Set False for a fully unrestrained run.
USE_RESTRAINTS = True
RESTRAINT_K = 1000.0         # kJ/mol/nm^2

DCD_INTERVAL_PS = 10.0       # trajectory frame spacing
LOG_INTERVAL_PS = 5.0        # log / terminal print spacing
CHECKPOINT_INTERVAL_PS = 100.0

PROTEIN_FF = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
SMALL_MOLECULE_FF = "gaff-2.11"

OUT_TRAJ = "complex_production.dcd"
OUT_LOG = "simulation.log"
OUT_TOPOLOGY_PDB = "system_solvated.pdb"
OUT_CHECKPOINT = "production.chk"
OUT_STATE = "final_state.xml"
OUT_FINAL_PDB = "final_frame.pdb"


# =============================================================================
# STEP 0: input check - runs BEFORE the slow imports below
# =============================================================================
def fail(message):
    """Print a clear, human-readable error and exit."""
    bar = "=" * 72
    print(f"\n{bar}\nERROR\n{bar}\n{message}\n{bar}\n", file=sys.stderr)
    sys.exit(1)


def check_inputs():
    cwd = os.getcwd()
    missing = [f for f in (PROTEIN_FILE, LIGAND_FILE)
               if not os.path.isfile(os.path.join(cwd, f))]
    if missing:
        listing = ", ".join(sorted(os.listdir(cwd))[:25]) or "(empty folder)"
        fail(
            "Required input file(s) not found in the current working directory:\n"
            f"    folder : {cwd}\n"
            f"    missing: {', '.join(missing)}\n\n"
            "Place 'protein.pdb' and 'ligand.sdf' in that folder (or cd to the folder\n"
            "that contains them) and run the script again.\n\n"
            f"Files currently there: {listing}"
        )
    empty = [f for f in (PROTEIN_FILE, LIGAND_FILE)
             if os.path.getsize(os.path.join(cwd, f)) == 0]
    if empty:
        fail(f"These input file(s) exist but are empty: {', '.join(empty)}")


check_inputs()

# =============================================================================
# Heavy imports (only reached if the input files exist)
# =============================================================================
try:
    import numpy as np
    from openmm import (CustomExternalForce, LangevinMiddleIntegrator,
                        MonteCarloBarostat, Platform)
    from openmm.app import (CheckpointReporter, DCDReporter, HBonds, Modeller,
                            PDBFile, PME, Simulation, StateDataReporter)
    from openmm.unit import (atmosphere, femtosecond, kelvin, kilojoule_per_mole,
                             molar, nanometer, picosecond)
    from openmmforcefields.generators import SystemGenerator
    from pdbfixer import PDBFixer
    try:
        from openff.toolkit import Molecule
    except ImportError:  # older openff-toolkit
        from openff.toolkit.topology import Molecule
except ImportError as exc:
    fail(
        f"A required Python package is missing: {exc}\n\n"
        "Create the environment with:\n"
        "  conda create -n openmm_md -c conda-forge python=3.11 openmm \\\n"
        "      openmmforcefields openff-toolkit ambertools pdbfixer rdkit numpy\n"
        "  conda activate openmm_md"
    )

# =============================================================================
# Helpers
# =============================================================================
_LOG_FH = None
DT_PS = TIMESTEP_FS / 1000.0


def say(message=""):
    """Print to terminal and to simulation.log."""
    print(message, flush=True)
    if _LOG_FH is not None:
        _LOG_FH.write(message + "\n")
        _LOG_FH.flush()


def steps_for(ps):
    return max(1, int(round(ps / DT_PS)))


def prepare_protein():
    """Clean the receptor with PDBFixer: add missing atoms/hydrogens, drop heterogens."""
    fixer = PDBFixer(filename=PROTEIN_FILE)
    fixer.findMissingResidues()
    fixer.missingResidues = {}              # do NOT invent loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    # Removes crystal ligand, waters, ions and cofactors. If you need a structural
    # metal (e.g. Zn2+) it must be handled separately - it is dropped here.
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(PH)
    return fixer.topology, fixer.positions


def load_ligand():
    """Read ligand.sdf, validate it, compute AM1-BCC charges. Returns (molecule, positions)."""
    mols = Molecule.from_file(LIGAND_FILE, file_format="sdf",
                              allow_undefined_stereo=True)
    if isinstance(mols, (list, tuple)):
        if len(mols) == 0:
            fail(f"No molecules could be read from {LIGAND_FILE}.")
        if len(mols) > 1:
            say(f"  WARNING: {LIGAND_FILE} has {len(mols)} molecules; using only the first.")
        mol = mols[0]
    else:
        mol = mols

    if mol.n_conformers == 0:
        fail(f"{LIGAND_FILE} contains no 3D coordinates.")
    if not any(atom.atomic_number == 1 for atom in mol.atoms):
        fail(f"{LIGAND_FILE} has no explicit hydrogens. Add hydrogens (e.g. with Open Babel "
             "'obabel in.sdf -O ligand.sdf -h' or RDKit AddHs) at the docked pose first.")

    conf = mol.conformers[0]
    try:
        from openff.units.openmm import to_openmm
        positions = to_openmm(conf)
    except Exception:
        positions = conf                      # already an OpenMM quantity (old toolkits)

    xyz = np.asarray(positions.value_in_unit(nanometer), dtype=float)
    if np.ptp(xyz[:, 2]) < 1e-4:
        say("  WARNING: ligand z-coordinates are all identical - the SDF looks 2D, "
            "not a 3D docked pose.")

    say("  Computing AM1-BCC partial charges (AmberTools/sqm, can take 1-3 min)...")
    try:
        mol.assign_partial_charges(partial_charge_method="am1bcc")
    except Exception as exc:
        fail(f"AM1-BCC charge generation failed: {exc}\n\n"
             "Make sure AmberTools is installed in this environment:\n"
             "  conda install -c conda-forge ambertools")
    return mol, positions


def add_position_restraints(system, topology, positions_nm, n_solute):
    """Harmonic restraints on heavy atoms of the solute (protein + ligand) only."""
    force = CustomExternalForce("restraint_k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("restraint_k", RESTRAINT_K)   # kJ/mol/nm^2
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    count = 0
    for atom in topology.atoms():
        if atom.index >= n_solute:      # solvent/ions were appended after the solute
            break
        if atom.element is not None and atom.element.symbol != "H":
            p = positions_nm[atom.index]
            force.addParticle(atom.index, [p[0], p[1], p[2]])
            count += 1
    system.addForce(force)
    return count


def run_phase(simulation, name, n_steps, extra_reporters=()):
    """Run one MD phase with terminal + simulation.log reporting."""
    simulation.reporters.clear()
    simulation.currentStep = 0          # keeps progress % and report intervals per-phase
    say(f"\n# ===== {name}: {n_steps} steps = {n_steps * DT_PS:.1f} ps =====")
    interval = steps_for(LOG_INTERVAL_PS)
    for out in (sys.stdout, _LOG_FH):
        simulation.reporters.append(StateDataReporter(
            out, interval,
            step=True, time=True, potentialEnergy=True, temperature=True,
            volume=True, density=True, speed=True,
            progress=True, remainingTime=True, totalSteps=n_steps,
            separator="\t"))
    for reporter in extra_reporters:
        simulation.reporters.append(reporter)
    simulation.step(n_steps)


# =============================================================================
# Main pipeline
# =============================================================================
def main():
    global _LOG_FH
    t_start = time.time()
    _LOG_FH = open(OUT_LOG, "w")

    say("OpenMM protein-ligand MD pipeline")
    say(f"  protein={PROTEIN_FILE}  ligand={LIGAND_FILE}")
    say(f"  T={TEMPERATURE_K} K  P={PRESSURE_ATM} atm  dt={TIMESTEP_FS} fs  "
        f"NVT={NVT_PS} ps  NPT={NPT_PS} ps  production={PRODUCTION_NS} ns")

    # ---- 1. Protein ---------------------------------------------------------
    say("\n[1/7] Preparing protein (PDBFixer: missing atoms, hydrogens at pH "
        f"{PH})...")
    prot_top, prot_pos = prepare_protein()
    modeller = Modeller(prot_top, prot_pos)
    prot_xyz = np.asarray(modeller.positions.value_in_unit(nanometer), dtype=float)
    say(f"  Protein atoms (with H): {modeller.topology.getNumAtoms()}")

    # ---- 2. Ligand ----------------------------------------------------------
    say("\n[2/7] Loading and charging ligand...")
    ligand, lig_pos = load_ligand()
    lig_xyz = np.asarray(lig_pos.value_in_unit(nanometer), dtype=float)
    min_d = np.linalg.norm(prot_xyz[:, None, :] - lig_xyz[None, :, :], axis=2).min()
    say(f"  Ligand atoms: {ligand.n_atoms}   closest protein contact: {min_d * 10:.1f} A")
    if min_d > 0.5:
        say("  WARNING: ligand is >5 A from every protein atom. It is probably NOT in the "
            "protein's coordinate frame / docked pose. The run will continue but the "
            "result will be meaningless.")

    # ---- 3. Force field + combine ------------------------------------------
    say(f"\n[3/7] Building force field (ff14SB + TIP3P + {SMALL_MOLECULE_FF})...")
    system_generator = SystemGenerator(
        forcefields=PROTEIN_FF,
        small_molecule_forcefield=SMALL_MOLECULE_FF,
        molecules=[ligand],
        forcefield_kwargs={"constraints": HBonds, "rigidWater": True,
                           "removeCMMotion": True},
        periodic_forcefield_kwargs={"nonbondedMethod": PME,
                                    "nonbondedCutoff": 1.0 * nanometer},
    )
    modeller.add(ligand.to_topology().to_openmm(), lig_pos)
    n_solute = modeller.topology.getNumAtoms()
    say(f"  Solute (protein + ligand) atoms: {n_solute}")

    # ---- 4. Solvate + ions --------------------------------------------------
    say(f"\n[4/7] Solvating: TIP3P, {BOX_PADDING_NM} nm padding ({BOX_SHAPE}), "
        f"neutralize + {IONIC_STRENGTH_M} M NaCl...")
    modeller.addSolvent(
        system_generator.forcefield,
        model="tip3p",
        padding=BOX_PADDING_NM * nanometer,
        boxShape=BOX_SHAPE,
        ionicStrength=IONIC_STRENGTH_M * molar,
        positiveIon="Na+",
        negativeIon="Cl-",
        neutralize=True,
    )
    say(f"  Total atoms: {modeller.topology.getNumAtoms()}")
    with open(OUT_TOPOLOGY_PDB, "w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh)
    say(f"  Wrote {OUT_TOPOLOGY_PDB} (use with the DCD for analysis)")

    # ---- 5. System + platform ----------------------------------------------
    say("\n[5/7] Creating system and OpenCL simulation...")
    system = system_generator.create_system(modeller.topology)

    if USE_RESTRAINTS:
        n_res = add_position_restraints(
            system, modeller.topology,
            modeller.positions.value_in_unit(nanometer), n_solute)
        say(f"  Heavy-atom restraints on {n_res} atoms (k = {RESTRAINT_K} kJ/mol/nm^2)")

    integrator = LangevinMiddleIntegrator(
        TEMPERATURE_K * kelvin, FRICTION_PER_PS / picosecond, TIMESTEP_FS * femtosecond)

    try:
        platform = Platform.getPlatformByName("OpenCL")
    except Exception:
        available = [Platform.getPlatform(i).getName()
                     for i in range(Platform.getNumPlatforms())]
        fail("The OpenMM 'OpenCL' platform is not available (required for the M4 GPU).\n"
             f"Platforms found: {available}\n\n"
             "Reinstall OpenMM from conda-forge for osx-arm64 and check with:\n"
             "  python -m openmm.testInstallation")
    # Apple GPUs have no FP64, so 'mixed'/'double' precision is not possible.
    properties = {"Precision": "single"}
    try:
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    except Exception as exc:
        fail(f"Could not create an OpenCL context: {exc}")
    simulation.context.setPositions(modeller.positions)
    say(f"  Platform: {simulation.context.getPlatform().getName()}  "
        f"(Precision={properties['Precision']})")

    # ---- 6. Minimization ----------------------------------------------------
    say("\n[6/7] Energy minimization (until convergence)...")
    e0 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(tolerance=10 * kilojoule_per_mole / nanometer,
                              maxIterations=0)
    e1 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    say(f"  Potential energy: {e0.value_in_unit(kilojoule_per_mole):.1f} -> "
        f"{e1.value_in_unit(kilojoule_per_mole):.1f} kJ/mol")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * kelvin)

    # ---- 7. Equilibration + production -------------------------------------
    say("\n[7/7] MD protocol")
    t_md = time.time()

    # NVT
    run_phase(simulation, f"NVT equilibration @ {TEMPERATURE_K} K", steps_for(NVT_PS))

    # NPT: add the barostat now, then re-initialise the context keeping the state
    simulation.system.addForce(MonteCarloBarostat(
        PRESSURE_ATM * atmosphere, TEMPERATURE_K * kelvin, 25))
    simulation.context.reinitialize(preserveState=True)
    run_phase(simulation, f"NPT equilibration @ {PRESSURE_ATM} atm", steps_for(NPT_PS))

    # Production: release restraints
    if USE_RESTRAINTS:
        simulation.context.setParameter("restraint_k", 0.0)
    n_prod = steps_for(PRODUCTION_NS * 1000.0)
    extras = [
        DCDReporter(OUT_TRAJ, steps_for(DCD_INTERVAL_PS)),
        CheckpointReporter(OUT_CHECKPOINT, steps_for(CHECKPOINT_INTERVAL_PS)),
    ]
    t_prod = time.time()
    run_phase(simulation, f"PRODUCTION NPT ({PRODUCTION_NS} ns, unrestrained)",
              n_prod, extra_reporters=extras)
    prod_seconds = time.time() - t_prod

    # ---- Save final state ---------------------------------------------------
    simulation.reporters.clear()
    simulation.saveCheckpoint(OUT_CHECKPOINT)
    simulation.saveState(OUT_STATE)
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    with open(OUT_FINAL_PDB, "w") as fh:
        PDBFile.writeFile(simulation.topology, state.getPositions(), fh)

    total_min = (time.time() - t_start) / 60.0
    prod_speed = PRODUCTION_NS / (prod_seconds / 86400.0)
    say("\nDONE")
    say(f"  Production speed (average): {prod_speed:.1f} ns/day")
    say(f"  Total wall time: {total_min:.1f} min "
        f"(MD part {(time.time() - t_md) / 60.0:.1f} min)")
    say(f"  Trajectory: {OUT_TRAJ}   topology: {OUT_TOPOLOGY_PDB}   log: {OUT_LOG}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
        sys.exit(130)
    except Exception as exc:
        text = str(exc)
        hint = ""
        if "No template found" in text or "template" in text.lower():
            hint = ("\n\nHint: this usually means protein.pdb has chain gaps / missing "
                    "residues, a non-standard residue, or a leftover heterogen "
                    "(metal, cofactor). Inspect the residue named in the message above.")
        elif "NaN" in text:
            hint = ("\n\nHint: 'Particle coordinate is NaN' = the system blew up. Usually "
                    "a bad ligand pose (severe clashes) or a wrong ligand protonation/"
                    "geometry. Check the ligand in your docking viewer.")
        fail(f"{type(exc).__name__}: {text}{hint}")