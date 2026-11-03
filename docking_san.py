import os
import re
import csv
import shutil
import hashlib
import argparse
import subprocess
import urllib.request
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.markup import escape

console = Console()

# ==============================================================================
# 1. CANDIDATE MOLECULES
#    NOTE: screenshot me NXI-13 / NXI-14A / NXI-14B bhi dikh rahe hain.
#    Apni local file ka CANDIDATES dict yahan paste kar do (ye sirf original list hai).
# ==============================================================================
CANDIDATES = {
    "M9 (Stage-0 Parent)": "O=C(NC1=CC=C(C(NCCCCCCC(NO)=O)=O)C=C1)N[C@]2(C3)C[C@]4(C)C[C@](C[C@H]3C4)(C)C2",
    "NXI-01 (Stage-1 Orig)": "O=C(NC1=CC=C(OCC(O)C(=O)O)C=C1)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2",
    "NXI-05-Clean (Stage-3A)": "O=C(CCc1ccc(NC(=O)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2)cc1)Nc5ccccn5",
    "NXI-06-Clean (Stage-3B)": "CC1CCC(NC(=O)Nc2ccc(CCC(=O)Nc3ccccn3)cc2)CC1",
    "NXI-07-Clean (Stage-3C)": "CC1CCC(NC(=O)Nc2ccc(CCC(=O)CO)cc2)CC1",
    "NXI-08-Clean (Stage-3D)": "CC1CCC(NC(=O)Nc2ccc(CCC(=O)N(C)O)cc2)CC1",
}

# ==============================================================================
# 2. TARGETS
#    WARNING: 5EDU zebrafish (Danio rerio) HDAC6 CD2 hai, human numbering
#    (610/611/649/651/782) usse match nahi karti. Script ab found residues ka
#    naam print karti hai; mismatch ho to ligand (TSN) proximity fallback chalega.
#    Residue numbers RCSB page se verify kar lena.
# ==============================================================================
MAIN_TARGETS = {
    "sEH": {"pdb_id": "3I28", "chains": ["A"], "ligand_names": ["I28"],
            "pocket_residues": [335, 383, 466], "keep_zn": False, "box_size": [22, 22, 22]},
    "HDAC6": {"pdb_id": "5EDU", "chains": ["A"], "ligand_names": ["TSN", "TSA"],
              "pocket_residues": [610, 611, 649, 651, 782], "keep_zn": True, "box_size": [22, 22, 22]},
}

REVERSE_TARGETS = {
    "hERG_Heart": {"pdb_id": "7CN1", "chains": ["A", "B", "C", "D"], "ligand_names": [],
                   "pocket_residues": [624, 652, 656], "keep_zn": False, "box_size": [24, 24, 24], "type": "toxic"},
    "CYP3A4_Liver": {"pdb_id": "1TQN", "chains": ["A"], "ligand_names": ["HEM"],
                     "pocket_residues": [304, 309, 369], "keep_zn": False, "box_size": [24, 24, 24], "type": "toxic"},
    "CYP2D6_Liver": {"pdb_id": "2F9Q", "chains": ["A"], "ligand_names": ["HEM"],
                     "pocket_residues": [216, 301, 483], "keep_zn": False, "box_size": [24, 24, 24], "type": "toxic"},
    "COX2_Inflam": {"pdb_id": "5IKR", "chains": ["A"], "ligand_names": ["ID8"],
                    "pocket_residues": [120, 355, 523], "keep_zn": False, "box_size": [22, 22, 22], "type": "toxic"},
    "HSA_Albumin": {"pdb_id": "1AO6", "chains": ["A"], "ligand_names": [],
                    "pocket_residues": [199, 214, 218, 222, 242, 257], "keep_zn": False,
                    "box_size": [24, 24, 24], "type": "pk"},
    "AChE_Brain": {"pdb_id": "4EY7", "chains": ["A"], "ligand_names": ["E20"],
                   "pocket_residues": [86, 203, 286, 337], "keep_zn": False, "box_size": [22, 22, 22], "type": "cns_bonus"},
    "MAO_B_Brain": {"pdb_id": "2V5Z", "chains": ["A"], "ligand_names": ["SAF"],
                    "pocket_residues": [172, 198, 326, 435], "keep_zn": False, "box_size": [22, 22, 22], "type": "cns_bonus"},
    "BACE1_Brain": {"pdb_id": "4DJU", "chains": ["A"], "ligand_names": ["0KP"],
                    "pocket_residues": [32, 228], "keep_zn": False, "box_size": [22, 22, 22], "type": "cns_bonus"},
}

WORK_DIR = "clean_docking_data"
os.makedirs(WORK_DIR, exist_ok=True)

EXHAUSTIVENESS = 8
# Vina exhaustiveness se zyada threads use nahi kar sakta (warning aati hai), isliye cap
N_CPU = min(os.cpu_count() or 4, EXHAUSTIVENESS)


# ==============================================================================
# 3. HELPERS
# ==============================================================================
def check_tools():
    missing = [t for t in ("obabel", "vina") if shutil.which(t) is None]
    if missing:
        raise SystemExit(f"Missing executables in PATH: {', '.join(missing)}")


def purge_cache(patterns):
    """Ab sirf --rebuild flag dene par cache delete hoga (pehle har run me delete ho jata tha)."""
    import glob
    for pattern in patterns:
        for fpath in glob.glob(os.path.join(WORK_DIR, pattern)):
            try:
                os.remove(fpath)
            except OSError:
                pass


def zn_pdbqt_line(xyz, serial=99999):
    """Strict fixed-column PDBQT line (name 13-16, resName 18-20, chain 22, resSeq 23-26,
    xyz 31-54, charge 71-76, type 78-79). Purani hand-typed line ke columns khisak gaye the
    (atom name / residue name / chain / coordinates sahi columns par nahi the), isliye Vina
    ne parse error diya aur HDAC6 har ligand par 'Failed' aaya."""
    return (f"HETATM{serial:5d} {' ZN '}{' '}{'ZN ':<3s} A{999:4d}{' '}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.00:6.2f}{0.00:6.2f}    {2.000:+6.3f} Zn\n")


def valid_pdbqt(path, min_size=300):
    return os.path.exists(path) and os.path.getsize(path) > min_size


# ==============================================================================
# 4. RECEPTOR BUILDER
# ==============================================================================
def prepare_clean_receptor(target_name, cfg):
    pdb_id = cfg["pdb_id"]
    allowed_chains = set(cfg["chains"])
    lig_names = set(cfg.get("ligand_names", []))
    pocket_res_nums = set(cfg.get("pocket_residues", []))

    raw_pdb = os.path.join(WORK_DIR, f"{pdb_id}_raw.pdb")
    clean_pdb = os.path.join(WORK_DIR, f"{target_name}_clean.pdb")
    rec_pdbqt = os.path.join(WORK_DIR, f"{target_name}_receptor.pdbqt")           # with Zn (if any)
    rec_nozn = os.path.join(WORK_DIR, f"{target_name}_receptor_noZn.pdbqt")       # fallback
    center_file = os.path.join(WORK_DIR, f"{target_name}_center.txt")

    if not os.path.exists(raw_pdb) or os.path.getsize(raw_pdb) < 1000:
        console.print(f"[cyan]Downloading {target_name} ({pdb_id}) from RCSB PDB...[/]")
        req = urllib.request.Request(f"https://files.rcsb.org/download/{pdb_id}.pdb",
                                     headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as resp, open(raw_pdb, "wb") as out:
            out.write(resp.read())

    cache_ok = valid_pdbqt(rec_pdbqt, 2000) and os.path.exists(center_file)
    if cache_ok and cfg.get("keep_zn") and not valid_pdbqt(rec_nozn, 2000):
        # purane (buggy) run ka cache: noZn backup nahi hai -> receptor dobara banao
        console.print(f"[yellow]  {target_name}: stale cache mila (noZn backup missing), rebuild ho raha hai[/]")
        cache_ok = False
    if cache_ok:
        with open(center_file) as f:
            center = [float(x) for x in f.read().strip().split(",")]
        return rec_pdbqt, (rec_nozn if valid_pdbqt(rec_nozn, 2000) else None), center

    atom_lines, pocket_res_coords, found_res = [], [], {}
    lig_coords_by_instance, zn_list = {}, []

    with open(raw_pdb, errors="ignore") as f:
        for line in f:
            if line.startswith("ENDMDL"):      # NMR/cryo-EM multi-model: sirf model 1
                break
            if line.startswith("ATOM"):
                ch = line[21].strip() or "A"
                if ch not in allowed_chains:
                    continue
                if line[16] not in (" ", "A", "1"):
                    continue
                atom_lines.append(line)
                try:
                    res_num = int(line[22:26].strip())
                    if res_num in pocket_res_nums:
                        pocket_res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        found_res[(ch, res_num)] = line[17:20].strip()
                except ValueError:
                    pass
            elif line.startswith("HETATM"):
                rname = line[17:20].strip()
                ch = line[21].strip() or "A"
                res_seq = line[22:26].strip()
                try:
                    xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                except ValueError:
                    continue
                if rname in lig_names:
                    lig_coords_by_instance.setdefault(f"{rname}_{ch}_{res_seq}", []).append(xyz)
                elif rname == "ZN" and cfg.get("keep_zn"):
                    zn_list.append(xyz)

    if not atom_lines:
        raise RuntimeError(f"{target_name}: no ATOM lines for chains {allowed_chains}")

    # sanity: pocket residues mile ya nahi (5EDU numbering issue pakadne ke liye)
    if pocket_res_nums:
        console.print(f"[dim]  {target_name}: pocket residues found -> "
                      f"{sorted(set(f'{n}{r}' for (_, n), r in found_res.items()))}[/]")
        if not found_res:
            console.print(f"[bold red]  WARNING {target_name}: koi pocket residue match nahi hua, "
                          f"numbering check karo![/]")

    # center
    if pocket_res_coords:
        res_center = np.mean(pocket_res_coords, axis=0)
        center = res_center
        if lig_coords_by_instance:
            best_c, best_d = None, 1e9
            for coords in lig_coords_by_instance.values():
                c = np.mean(coords, axis=0)
                d = np.linalg.norm(c - res_center)
                if d < best_d:
                    best_c, best_d = c, d
            if best_d < 15.0:
                center = best_c
    elif lig_coords_by_instance:
        center = np.mean(next(iter(lig_coords_by_instance.values())), axis=0)
    else:
        center = np.mean([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in atom_lines], axis=0)
    center = np.array(center, dtype=float)

    if "HEM" in lig_names:
        center[2] += 3.0

    with open(clean_pdb, "w") as out:
        out.writelines(atom_lines)
        out.write("END\n")

    console.print(f"[yellow]Preparing {target_name:<13} center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})[/]")

    # -h: polar hydrogens add (pehle H add nahi ho rahe the)
    res = subprocess.run(["obabel", clean_pdb, "-xr", "-h", "-O", rec_pdbqt, "--partialcharge", "gasteiger"],
                         capture_output=True, text=True)
    if not valid_pdbqt(rec_pdbqt, 2000):
        raise RuntimeError(f"{target_name}: obabel receptor conversion failed: {res.stderr[-300:]}")

    clean_lines = []
    with open(rec_pdbqt, errors="ignore") as f:
        for l in f:
            if l.startswith(("ATOM", "HETATM")):
                l = re.sub(r"(?i)[-+]?\bnan\b", " 0.000", l)
                clean_lines.append(l)

    # fallback receptor (Zn ke bina) hamesha save
    with open(rec_nozn, "w") as f:
        f.writelines(clean_lines)

    if cfg.get("keep_zn") and zn_list:
        dists = [np.linalg.norm(np.array(z) - center) for z in zn_list]
        cat_zn = zn_list[int(np.argmin(dists))]
        clean_lines.append(zn_pdbqt_line(cat_zn))
        console.print(f"[dim]  catalytic Zn at ({cat_zn[0]:.1f}, {cat_zn[1]:.1f}, {cat_zn[2]:.1f}), "
                      f"{min(dists):.1f} A from center[/]")
    elif cfg.get("keep_zn"):
        console.print(f"[bold red]  WARNING {target_name}: Zn HETATM nahi mila[/]")

    with open(rec_pdbqt, "w") as f:
        f.writelines(clean_lines)
    with open(center_file, "w") as f:
        f.write(f"{center[0]:.3f},{center[1]:.3f},{center[2]:.3f}")

    return rec_pdbqt, (rec_nozn if cfg.get("keep_zn") else None), center.tolist()


# ==============================================================================
# 5. LIGAND BUILDER
# ==============================================================================
def lig_safe_name(name, smiles):
    # hash add kiya: SMILES badle to purana cached pdbqt reuse nahi hoga
    base = re.sub(r"[^a-zA-Z0-9_]", "_", name.split()[0])
    return f"{base}_{hashlib.md5(smiles.encode()).hexdigest()[:6]}"


def prepare_ligand_3d(name, smiles):
    safe = lig_safe_name(name, smiles)
    sdf = os.path.join(WORK_DIR, f"{safe}.sdf")
    pdbqt = os.path.join(WORK_DIR, f"{safe}.pdbqt")

    if valid_pdbqt(pdbqt):
        return pdbqt

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        console.print(f"[bold red]Invalid SMILES for {name}[/]")
        return None
    Chem.RemoveStereochemistry(mol)   # NOTE: stereo hata raha hai (original behaviour); chahiye to hata do
    mol_3d = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.maxIterations = 500
    params.randomSeed = 42
    cid = AllChem.EmbedMolecule(mol_3d, params)

    if cid != -1:
        try:
            AllChem.MMFFOptimizeMolecule(mol_3d, confId=cid, maxIters=500)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol_3d, confId=cid, maxIters=500)
            except Exception:
                pass
        w = Chem.SDWriter(sdf)
        w.write(mol_3d, confId=cid)
        w.close()
        subprocess.run(["obabel", sdf, "-O", pdbqt, "-p", "7.4", "--partialcharge", "gasteiger"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        subprocess.run(["obabel", f"-:{smiles}", "-O", pdbqt, "--gen3d", "-p", "7.4",
                        "--partialcharge", "gasteiger"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return pdbqt if valid_pdbqt(pdbqt) else None


# ==============================================================================
# 6. VINA ENGINE
# ==============================================================================
def parse_pose_score(pose_file):
    if not os.path.exists(pose_file) or os.path.getsize(pose_file) < 100:
        return None
    with open(pose_file, errors="ignore") as f:
        for line in f:
            if "REMARK VINA RESULT:" in line:
                m = re.search(r"REMARK VINA RESULT:\s+([-+]?\d+\.\d+)", line)
                if m:
                    return float(m.group(1))
    return None


def run_vina(lig_pdbqt, rec_pdbqt, center, box, out_pose):
    cmd = ["vina", "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{center[0]:.2f}", "--center_y", f"{center[1]:.2f}", "--center_z", f"{center[2]:.2f}",
           "--size_x", str(box[0]), "--size_y", str(box[1]), "--size_z", str(box[2]),
           "--exhaustiveness", str(EXHAUSTIVENESS), "--cpu", str(N_CPU), "--seed", "42",
           "--out", out_pose]
    res = subprocess.run(cmd, capture_output=True, text=True)
    sc = parse_pose_score(out_pose)
    if sc is None:
        for line in res.stdout.splitlines():
            m = re.match(r"^\s*1\s+([-+]?\d+\.\d+)", line)
            if m:
                sc = float(m.group(1))
                break
    # WARNING lines asli error nahi hoti, unhe skip karke last real error line dikhao
    lines = (res.stderr.strip() + "\n" + res.stdout.strip()).splitlines()
    real = [l for l in lines if l.strip() and "WARNING" not in l.upper()]
    return sc, (real[-1].strip() if real else "unknown error")


def dock_ligand(lig_pdbqt, rec, center, box, out_pose):
    cached = parse_pose_score(out_pose)
    if cached is not None:
        return cached, ""

    sc, err = run_vina(lig_pdbqt, rec["pdbqt"], center, box, out_pose)
    if sc is not None:
        return sc, ""

    if rec.get("fallback"):     # Zn ki wajah se fail -> bina Zn ke retry
        sc, err2 = run_vina(lig_pdbqt, rec["fallback"], center, box, out_pose)
        if sc is not None:
            return sc, "noZn"
        err = err2
    console.print(f"[dim red]({escape(err)})[/]", end=" ")   # escape: rich brackets ko markup na samjhe
    return None, ""


# ==============================================================================
# 7. MAIN
# ==============================================================================
def fmt(v):
    return "N/A" if v is None else f"{v:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", nargs="*", metavar="PATTERN",
                    help="cache delete karo, e.g. --rebuild HDAC6 hERG (bina arg = sab)")
    args = ap.parse_args()
    check_tools()
    if args.rebuild is not None:
        pats = [f"*{p}*" for p in args.rebuild] if args.rebuild else ["*"]
        purge_cache(pats)

    console.print(Panel.fit(
        "[bold cyan]STEP 2: Dual-Target (sEH & HDAC6) + 8-Receptor Off-Target Reverse Docking[/]",
        border_style="bright_blue"))

    all_targets = {**MAIN_TARGETS, **REVERSE_TARGETS}
    receptors = {}
    for t_name, cfg in all_targets.items():
        try:
            r_pdbqt, fb, center = prepare_clean_receptor(t_name, cfg)
            receptors[t_name] = {"pdbqt": r_pdbqt, "fallback": fb, "center": center, "box": cfg["box_size"]}
        except Exception as e:
            console.print(f"[bold red]Receptor {t_name} failed: {escape(str(e))}[/]")

    results = {}
    for name, smiles in CANDIDATES.items():
        console.print(f"\n[bold green]Docking Candidate:[/] {name}")
        lig = prepare_ligand_3d(name, smiles)
        scores = {t: None for t in all_targets}
        if lig is None:
            console.print("[bold red]   ligand prep failed[/]")
            results[name] = scores
            continue
        safe = lig_safe_name(name, smiles)
        for t_name, r in receptors.items():
            pose = os.path.join(WORK_DIR, f"{safe}_on_{t_name}.pdbqt")
            console.print(f"   - {t_name:<14} -> ", end="")
            sc, note = dock_ligand(lig, r, r["center"], r["box"], pose)
            scores[t_name] = sc
            if sc is None:
                console.print("[bold red]Failed[/]")
            else:
                extra = " [yellow](no Zn)[/]" if note == "noZn" else ""
                console.print(f"[bold cyan]{sc:.2f} kcal/mol[/]{extra}")
        results[name] = scores

    t_main = Table(title="\nDual-Target Potency, Balance & hERG Selectivity", expand=True, border_style="cyan")
    t_main.add_column("Candidate", style="bold white", width=23)
    for c in ("sEH (3I28)", "HDAC6 (5EDU)", "Dual Average", "Balance |sEH-HDAC6|", "hERG margin (hERG-sEH)"):
        t_main.add_column(c, justify="center")

    for name, sc in results.items():
        seh, hdac, herg = sc.get("sEH"), sc.get("HDAC6"), sc.get("hERG_Heart")
        s_str = f"[bold green]{seh:.2f}[/]" if seh is not None and seh <= -8.5 else fmt(seh)
        h_str = f"[bold green]{hdac:.2f}[/]" if hdac is not None and hdac <= -7.5 else fmt(hdac)
        if seh is not None and hdac is not None:
            diff = abs(seh - hdac)
            avg_s = f"[bold bright_green]{(seh + hdac) / 2:.2f}[/]"
            bal_s = f"[bold green]{diff:.2f} (Balanced)[/]" if diff <= 2.0 else f"[bold yellow]{diff:.2f}[/]"
        else:
            avg_s = bal_s = "N/A"
        if seh is not None and herg is not None:
            m = herg - seh
            mar_s = f"[bold green]+{m:.2f} (Safe)[/]" if m >= 2.0 else f"[bold red]{m:+.2f} (Narrow)[/]"
        else:
            mar_s = "N/A"
        t_main.add_row(name, s_str, h_str, avg_s, bal_s, mar_s)
    console.print(t_main)

    t_rev = Table(title="\n8-Receptor Reverse Docking Off-Target Panel (kcal/mol)", expand=True, border_style="magenta")
    t_rev.add_column("Candidate", style="bold white", width=22)
    for rt in REVERSE_TARGETS:
        t_rev.add_column(rt.replace("_", "\n"), justify="center")
    for name, sc in results.items():
        row = [name]
        for rt, cfg in REVERSE_TARGETS.items():
            val = sc.get(rt)
            if val is None:
                row.append("N/A")
            elif cfg["type"] == "toxic":
                row.append(f"[bold green]{val:.1f} (Safe)[/]" if val > -7.2 else f"[bold red]{val:.1f} (Risk)[/]")
            else:
                row.append(f"[bold cyan]{val:.1f}[/]")
        t_rev.add_row(*row)
    console.print(t_rev)

    with open("Stage3_Full_Docking_Report.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Candidate"] + list(all_targets.keys()))
        w.writeheader()
        for name, sc in results.items():
            w.writerow({"Candidate": name, **sc})
    console.print("[bold green]Saved scores to Stage3_Full_Docking_Report.csv[/]")


if __name__ == "__main__":
    main()