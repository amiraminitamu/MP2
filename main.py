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
from pyscf import gto, scf, mp
from molecule import build_h2
from scf_solver import scf_cycle
from mp2_correction import compute_mp2_energy
from plot_utils import plot_energies
from output import initialize_output, write_step_output, write_table_header, write_table_row, write_summary

def main():
    distances = np.linspace(0.5, 2.0, 10)
    scf_steps = 100
    tolerance = 1e-10

    hf_energies = []
    mp2_corrections = []
    pyscf_mp2_totals = []

    output_file = "out.txt"
    initialize_output(output_file)

    for d in distances:
        mol = gto.Mole()
        mol.atom = build_h2(d)
        mol.spin = 0
        mol.basis = '6-31g'
        mol.build()

        mf = scf.RHF(mol)
        H = mf.get_hcore()
        S = mol.intor('int1e_ovlp')
        E_nuc = mol.energy_nuc()
        nbf = mol.nao_nr()

        hf_energy, mos, mo_energy = scf_cycle(H, S, E_nuc, mf, nbf, scf_steps, tolerance, mol)
        mp2_energy = compute_mp2_energy(mol, mos, mo_energy)
        total_energy = hf_energy + mp2_energy

        mf_pyscf = scf.RHF(mol).run()
        mp2_result = mp.MP2(mf_pyscf).run()
        pyscf_total = mp2_result.e_tot

        # Write detailed block and prepare for table
        write_step_output(output_file, d, hf_energy, mp2_energy, total_energy, pyscf_total, mo_energy)
        hf_energies.append(hf_energy)
        mp2_corrections.append(mp2_energy)
        pyscf_mp2_totals.append(pyscf_total)

    # Write summary table and final analysis
    write_table_header(output_file)
    for i in range(len(distances)):
        write_table_row(output_file, distances[i], hf_energies[i], mp2_corrections[i],
                        hf_energies[i] + mp2_corrections[i], pyscf_mp2_totals[i])
    write_summary(output_file, distances, [h + c for h, c in zip(hf_energies, mp2_corrections)],
                  mp2_corrections, pyscf_mp2_totals)

    # Plot the energies and save as plt.png
    plot_energies(distances, hf_energies, mp2_corrections, pyscf_mp2_totals)

if __name__ == "__main__":
    main()
