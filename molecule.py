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


def build_h2(distance):
    return f'''
    H 0.00000 0.00000 0.00000
    H {distance:.5f} 0.00000 0.00000
    '''

