import csv
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

console = Console()

# ==============================================================================
# STAGE 0, STAGE 1, AND STAGE 3 CANDIDATES FOR COMPLETE COMPARISON
# ==============================================================================
CANDIDATES = {
    "NXI-07 (Baseline)": "CC1CCC(NC(=O)Nc2ccc(CCC(=O)CO)cc2)CC1",
    "NXI-12-Rigid": "O=C(/C=C/c1ccc(NC(=O)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2)cc1)CO",
    "NXI-13-Rigid": "O=C(Cc1ccc(NC(=O)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2)cc1)CO",
    "NXI-14A-Direct": "O=C(CO)c1ccc(NC(=O)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2)cc1",
    "NXI-14B-Indane": "O=C(CO)C1Cc2ccc(NC(=O)N[C@]3(C4)C[C@]5(C)C[C@H]4C5(C)C3)cc2C1"
}

# Initialize PAINS and Brenk Structural Toxicophore Filters
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
tox_catalog = FilterCatalog(params)


def calculate_cns_mpo_and_logbb(mw, clogp, tpsa, hbd):
    """
    Calculates:
    1. Clark's LogBB equation: LogBB = -0.0148 * TPSA + 0.152 * cLogP + 0.139
       (LogBB > 0.0: High Brain Penetration, LogBB > -0.3: Moderate/Good, < -1.0: Poor)
    2. Wager Pfizer CNS MPO Score (Scaled out of 6.0 based on MW, LogP, TPSA, HBD)
    """
    log_bb = (-0.0148 * tpsa) + (0.152 * clogp) + 0.139

    # Component scores (0.0 to 1.0 each)
    s_clogp = 1.0 if clogp <= 3.0 else max(0.0, 1.0 - (clogp - 3.0) / 2.0)
    s_mw = 1.0 if mw <= 360 else max(0.0, 1.0 - (mw - 360) / 140.0)
    s_tpsa = 1.0 if 40 <= tpsa <= 90 else (max(0.0, 1.0 - (tpsa - 90) / 30.0) if tpsa > 90 else tpsa / 40.0)
    s_hbd = 1.0 if hbd == 0 else max(0.0, 1.0 - (hbd / 3.5))
    
    # Estimate neutral/non-ionizable bonus (Carboxylic acids ionize at pH 7.4 -> heavy penalty)
    cns_mpo = (s_clogp + s_mw + s_tpsa + s_hbd) * 1.5  # Scaled to 6.0
    return log_bb, min(6.0, cns_mpo)


def main():
    console.print(Panel.fit(
        "[bold cyan]🔬 STEP 1: Deep ADMET, Advance BBB (LogBB & CNS MPO), and Toxicophore Audit[/]\n"
        "[dim]Validating Stage-3 Leads against M9 & NXI-01 using RDKit Medicinal Chemistry Filters[/]",
        border_style="bright_blue"
    ))

    table1 = Table(title="\n🧠 TABLE 1A: Basic & Advance Blood-Brain Barrier (BBB) Penetration Profile", expand=True, border_style="cyan")
    table1.add_column("Candidate", style="bold white", width=23)
    table1.add_column("MW\n(<500)", justify="center")
    table1.add_column("TPSA\n(<90 Å²)", justify="center")
    table1.add_column("cLogP\n(2.0-4.0)", justify="center")
    table1.add_column("HBD / HBA\n(≤3 / ≤7)", justify="center")
    table1.add_column("RotB\n(≤8)", justify="center")
    table1.add_column("Clark LogBB\n(>-0.3 Ideal)", justify="center")
    table1.add_column("Pfizer CNS MPO\n(≥4.0 / 6.0)", justify="center")
    table1.add_column("BOILED-Egg\nBBB Status", justify="center")

    table2 = Table(title="\n💊 TABLE 1B: ADMET, Drug-Likeness (QED), fsp3 & Structural Toxicophore Alerts", expand=True, border_style="magenta")
    table2.add_column("Candidate", style="bold white", width=23)
    table2.add_column("QED Score\n(>0.5 Ideal)", justify="center")
    table2.add_column("Fraction Csp3\n(>0.25 3D-ness)", justify="center")
    table2.add_column("Lipinski / Veber\nViolations", justify="center")
    table2.add_column("Ionization Trap\n(at pH 7.4)", justify="center")
    table2.add_column("PAINS / Brenk / NIH\nToxicophore Alerts", justify="center")
    table2.add_column("Final ADMET Verdict", justify="center")

    csv_rows = []

    for name, smiles in CANDIDATES.items():
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        tpsa = Descriptors.TPSA(mol)
        clogp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rotb = Lipinski.NumRotatableBonds(mol)
        qed = QED.qed(mol)
        fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)

        log_bb, cns_mpo = calculate_cns_mpo_and_logbb(mw, clogp, tpsa, hbd)

        # Check for Carboxylic Acid ionization trap at pH 7.4
        has_cooh = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))
        if has_cooh:
            cns_mpo -= 1.5  # Severe BBB penalty for anionic carboxylate at pH 7.4

        # BOILED-Egg BBB Yolk Criteria (TPSA < 79-90 and WLogP between 0.4 and 6.0)
        boiled_egg_bbb = (tpsa < 90.0) and (1.5 <= clogp <= 4.5) and (not has_cooh)

        # Check PAINS / Brenk / NIH Toxicophore Alerts
        matches = tox_catalog.GetMatches(mol)
        alert_names = [m.GetDescription() for m in matches]
        n_alerts = len(alert_names)

        # Lipinski & Veber Violations
        lip_v = sum([mw > 500, clogp > 5.0, hbd > 5, hba > 10])
        veb_v = sum([rotb > 10, tpsa > 140])

        # Format Table 1A Row
        tpsa_s = f"[bold green]{tpsa:.1f}[/]" if tpsa < 90 else f"[bold red]{tpsa:.1f}[/]"
        logp_s = f"[bold green]{clogp:.2f}[/]" if 2.0 <= clogp <= 4.2 else f"[bold yellow]{clogp:.2f}[/]"
        rot_s = f"[bold green]{rotb}[/]" if rotb <= 8 else f"[bold red]{rotb}[/]"
        logbb_s = f"[bold green]{log_bb:+.2f}[/]" if log_bb > -0.35 else f"[bold red]{log_bb:+.2f}[/]"
        mpo_s = f"[bold green]{cns_mpo:.2f} / 6.0[/]" if cns_mpo >= 4.0 else f"[bold red]{cns_mpo:.2f} / 6.0[/]"
        egg_s = "[bold bright_green]✔ Yolk (Permeable)[/]" if boiled_egg_bbb else "[bold red]✖ Impermeable[/]"

        table1.add_row(name, f"{mw:.1f}", tpsa_s, logp_s, f"{hbd} / {hba}", rot_s, logbb_s, mpo_s, egg_s)

        # Format Table 1B Row
        qed_s = f"[bold green]{qed:.2f}[/]" if qed >= 0.50 else f"[bold yellow]{qed:.2f}[/]"
        fsp3_s = f"[bold green]{fsp3:.2f}[/]" if fsp3 >= 0.30 else f"{fsp3:.2f}"
        ion_s = "[bold red]✖ Charged (-COO⁻)[/]" if has_cooh else "[bold green]✔ Neutral (Unionized)[/]"
        alert_s = "[bold green]0 Alerts (Clean)[/]" if n_alerts == 0 else f"[bold red]{n_alerts}: {', '.join(alert_names[:2])}[/]"
        verdict = "[bold bright_green]★ STAGE-3 GOLD[/]" if (boiled_egg_bbb and n_alerts == 0 and cns_mpo >= 4.0) else "[bold red]✖ REJECT[/]"

        table2.add_row(name, qed_s, fsp3_s, f"{lip_v} / {veb_v}", ion_s, alert_s, verdict)

        csv_rows.append({
            "Candidate": name, "MW": round(mw, 2), "TPSA": round(tpsa, 2),
            "cLogP": round(clogp, 2), "HBD": hbd, "HBA": hba, "RotB": rotb,
            "LogBB": round(log_bb, 3), "CNS_MPO": round(cns_mpo, 2),
            "QED": round(qed, 3), "fsp3": round(fsp3, 2),
            "Tox_Alerts": "; ".join(alert_names) if alert_names else "Clean"
        })

    console.print(table1)
    console.print(table2)

    with open("Stage3_ADMET_BBB_Audit.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_rows[0].keys())
        writer.writeheader()
        writer.writerows(csv_rows)
    console.print("[bold green]✔ Saved full ADMET & Advance BBB report to Stage3_ADMET_BBB_Audit.csv[/]")


if __name__ == "__main__":
    main()