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
import scipy as sp
from density import compute_density_matrix

def scf_cycle(H, S, E_nuc, mf, nbf, scf_steps, tolerance, mol):
    dm = mf.get_init_guess()
    E_last = 0.0

    for _ in range(scf_steps):
        n_occ = mol.nelectron // 2
        J = mf.get_j(dm=dm)
        K = mf.get_k(dm=dm)
        F = H + J - 0.5 * K
        S_inv_sqrt = sp.linalg.sqrtm(np.linalg.inv(S))
        F_ortho = S_inv_sqrt @ F @ S_inv_sqrt
        eigvals, eigvecs = np.linalg.eigh(F_ortho)
        mos = S_inv_sqrt @ eigvecs
        E_HF = 0.5 * np.einsum('ij,ij', dm, H + F) + E_nuc
        if abs(E_HF - E_last) < tolerance:
            break
        dm = compute_density_matrix(mos, n_occ)
        E_last = E_HF

    return E_HF, mos, eigvals

