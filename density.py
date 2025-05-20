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

def compute_density_matrix(mos, n_occ):
    nbasis = mos.shape[0]
    dm = np.zeros((nbasis, nbasis))
    for i in range(nbasis):
        for j in range(nbasis):
            for oo in range(n_occ):
                dm[i, j] += 2.0 * np.real(mos[i, oo] * mos[j, oo])
    return dm

