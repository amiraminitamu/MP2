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
from pyscf import ao2mo

def compute_mp2_energy(mol, mo_coeff, mo_energy):
    # Number of spatial molecular orbitals
    nmo = mo_coeff.shape[1]
    
    # Number of occupied spatial orbitals (assuming closed-shell RHF)
    nocc = mol.nelectron // 2
    
    # Total number of spin orbitals (each spatial orbital gives α and β)
    nso = 2 * nmo

    # Create spin-orbital energies: duplicate each spatial energy
    eps = np.zeros(nso)
    for p in range(nso):
        eps[p] = mo_energy[p // 2]

    # Compute two-electron integrals in AO basis, then transform to MO basis
    eri_ao = mol.intor('int2e')               # (μν|λσ) in AO basis
    eri_mo = ao2mo.full(mol, mo_coeff)        # (ij|kl) in MO basis, flattened
    eri_mo = ao2mo.restore(1, eri_mo, nmo)     # Restore full (ij|kl) tensor shape

    # Build antisymmetrized spin-orbital integrals:
    # <pq||rs> = <pr|qs> - <ps|qr>, with same-spin constraints
    eri_spin = np.zeros((nso, nso, nso, nso))
    for p in range(nso):
        for q in range(nso):
            for r in range(nso):
                for s in range(nso):
                    # Only same-spin contributions are non-zero
                    val1 = eri_mo[p//2, r//2, q//2, s//2] * (p%2 == r%2) * (q%2 == s%2)
                    val2 = eri_mo[p//2, s//2, q//2, r//2] * (p%2 == s%2) * (q%2 == r%2)
                    eri_spin[p, q, r, s] = val1 - val2

    # Compute MP2 correlation energy using spin-orbital formulation
    E_mp2 = 0.0
    for i in range(nocc * 2):                # Occupied spin orbitals
        for j in range(nocc * 2):
            for a in range(nocc * 2, nso):   # Virtual spin orbitals
                for b in range(nocc * 2, nso):
                    denom = eps[i] + eps[j] - eps[a] - eps[b]
                    if abs(denom) < 1e-10:   # Avoid division by near-zero
                        continue
                    E_mp2 += 0.25 * eri_spin[i, j, a, b]**2 / denom

    return E_mp2  # MP2 correlation energy (not total energy)
