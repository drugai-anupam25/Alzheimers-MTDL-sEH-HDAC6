import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.align import Align

def display_executive_summary():
    console = Console()
    
    # CSV file read karna (jo teri pichli script ne banayi thi)
    try:
        df = pd.read_csv("Off_Target_Profile.csv")
    except FileNotFoundError:
        console.print("[bold red]Error:[/bold red] 'Off_Target_Profile.csv' file nahi mili. Pehle inverse docking run karo.")
        return

    # Ek solid title aur table create karna
    table = Table(
        show_header=True, 
        header_style="bold cyan", 
        expand=True,
        border_style="bright_blue"
    )
    
    # Columns add karna
    table.add_column("Lead Molecule", style="bold magenta", justify="center", width=15)
    table.add_column("Clinical Off-Target Profile", style="bold green", width=30)
    table.add_column("Binding Affinity (kcal/mol)", justify="right", style="bold yellow", width=20)
    table.add_column("Status / Risk", justify="center", width=20)

    # Data ko table me format karna
    for index, row in df.iterrows():
        ligand = str(row['Ligand'])
        target = str(row['Off-Target'])
        affinity = float(row['Affinity (kcal/mol)'])
        
        # Risk analysis logic (Blind docking artifacts ko context dena)
        if "HSA_Plasma_Albumin" in target:
            status = "[bold blue]Safe (Transport)[/bold blue]"
        elif affinity < -9.5:
            status = "[bold red]High Interaction[/bold red]"
        else:
            status = "[bold white]Artifact / Low Risk[/bold white]"

        table.add_row(ligand, target, f"{affinity:.3f}", status)
        
        # Ek separator add karna har naye ligand ke baad
        if index < len(df) - 1 and df.iloc[index]['Ligand'] != df.iloc[index+1]['Ligand']:
            table.add_section()

    # Final output ko box me wrap karke print karna
    panel = Panel(
        Align.center(table), 
        title="[bold yellow]CSIR-CDRI Pipeline: Off-Target Blind Docking Validation[/bold yellow]", 
        subtitle="[dim]Apple M4 Local Simulation[/dim]",
        border_style="green"
    )
    
    console.print()
    console.print(panel)
    console.print()

if __name__ == "__main__":
    display_executive_summary()