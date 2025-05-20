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


import matplotlib.pyplot as plt

def plot_energies(distances, hf_energies, mp2_corrections, pyscf_mp2_totals):
    total_energies = [hf + mp2 for hf, mp2 in zip(hf_energies, mp2_corrections)]

    plt.figure()
    plt.plot(distances, hf_energies, '-', label='Custom HF Energy')
    plt.plot(distances, total_energies, '-', label='Custom MP2 Total Energy')
    plt.plot(distances, pyscf_mp2_totals, '--', label='PySCF MP2 Total Energy')
    plt.xlabel('H–H Distance (Å)')
    plt.ylabel('Energy (Hartree)')
    plt.title('H₂ Potential Energy Curve (HF and MP2)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Save the figure
    plt.savefig("plt.png", dpi=300)
    print("Plot saved as plt.png")

    # Optional: Show the plot interactively
    plt.show()
