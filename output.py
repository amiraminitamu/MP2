"""
Author: Amir Amini
Email: amiramini@tamu.edu
Affiliation: PhD Student, Texas A&M University
Date: May 2025

Description:
This is a small computational chemistry code written in Python to calculate 
MP2 (Møller–Plesset perturbation theory) correlation energy for a molecule 
and generate potential energy surface (PES) plots.

Note:
This code is strictly intended for educational and prototyping purposes. 
It is not optimized for production-scale calculations.

License: MIT
"""


import numpy as np
from datetime import datetime

def initialize_output(file="out.txt"):
    with open(file, "w") as f:
        f.write("********************************************\n")
        f.write("*************Author: Amir Amini*************\n")
        f.write("============================================\n")
        f.write("      H₂ Potential Energy Curve Report\n")
        f.write("      Method: RHF + MP2 (custom & PySCF)\n")
        f.write(f"      Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("============================================\n\n")

def write_step_output(file, distance, hf_energy, mp2_corr, total_energy, pyscf_mp2_total, mo_energies):
    with open(file, "a") as f:
        f.write(f"Bond Distance: {distance:.2f} Å\n")
        f.write("-" * 60 + "\n")
        f.write("Molecular Orbital Energies (ε):\n")
        for i, e in enumerate(mo_energies):
            f.write(f"  Orbital {i+1:2d}: {e: .10f} Hartree\n")
        f.write("\n")
        f.write("Energy Summary:\n")
        f.write(f"  HF Energy (custom):        {hf_energy: .10f} Hartree\n")
        f.write(f"  MP2 Correlation (custom):  {mp2_corr: .10f} Hartree\n")
        f.write(f"  MP2 Total Energy (custom): {total_energy: .10f} Hartree\n")
        f.write(f"  MP2 Total Energy (PySCF):  {pyscf_mp2_total: .10f} Hartree\n")
        f.write(f"  Δ(Custom − PySCF):         {total_energy - pyscf_mp2_total: .3e} Hartree\n")
        f.write("-" * 60 + "\n\n")

def write_table_header(file="out.txt"):
    with open(file, "a") as f:
        f.write("Summary Table:\n")
        f.write(f"{'Distance (Å)':>12} | {'HF':>14} | {'MP2 Corr':>14} | {'MP2 Total':>14} | {'PySCF MP2':>14} | {'Δ':>10}\n")
        f.write("-" * 80 + "\n")

def write_table_row(file, d, hf, corr, total, pyscf_total):
    delta = total - pyscf_total
    with open(file, "a") as f:
        f.write(f"{d:12.2f} | {hf:14.8f} | {corr:14.8f} | {total:14.8f} | {pyscf_total:14.8f} | {delta:10.2e}\n")

def write_summary(file, distances, total_custom, mp2_corrs, pyscf_totals):
    min_idx = int(np.argmin(total_custom))
    eq_distance = distances[min_idx]
    eq_energy = total_custom[min_idx]
    avg_diff = np.mean(np.abs(np.array(total_custom) - np.array(pyscf_totals)))

    with open(file, "a") as f:
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("Final Summary:\n")
        f.write(f"  Equilibrium bond length:       {eq_distance:.9f} Å\n")
        f.write(f"  Minimum MP2 total energy:      {eq_energy:.10f} Hartree\n")
        f.write(f"  MP2 energy gain at equilibrium:{mp2_corrs[min_idx]:.10f} Hartree\n")
        f.write(f"  Average |Δ(Custom − PySCF)|:   {avg_diff:.6e} Hartree\n")
        f.write("=" * 80 + "\n")

