#########################################################################################################
# Type of Program: Free Two-Body system under penetrable spatial confinement
# Author: Mr. Koustav Das Chakladar, Dr. Santanu Mondal, Dr. Jayanta Kr. Saha.
# Methodology: Lagrange-Laguerre-Mesh Method
# Reference: 1. Physical Review C 67, 044309 (2003)
#            2. Phys. Stat. Sol (b) 243, No. 5, 1095-1109 (2006)
#            3. Mach. Learn: Sci. Techno. 4 (2023) 015024
#########################################################################################################
import numpy as np
import sys
import json
from scipy.special import factorial
from scipy.special import roots_laguerre
from scipy.linalg  import eigh
from numba import njit, prange
import os
from modulecalculation import calculation
from datetime import datetime
Now=datetime.now()
######################################################################################################### 
# Controlling Variables
#########################################################################################################
N  = 100                             # Specify the degree of the Lagurre polynomials N
x, w , mu = roots_laguerre(N, True)  # Get the roots and weights
lamda =  w                           # lambda =  exp(xi)/ xi * L'n(xi)^2 
########################################################################################################
n0 = int(sys.argv[1])
Z = int(sys.argv[2])
l = int(sys.argv[3])
V0 = float(sys.argv[4])
D0 = float(sys.argv[5])
x0 = float(sys.argv[6])
V1 = float(sys.argv[7])
D1 = float(sys.argv[8])
x1 = float(sys.argv[9])
pk = int(sys.argv[10])


print("Running python")

# n0 = 1                               # Principle Quantum Number
# Z  = 1                               # Atomic Number
# l  = 0                               # Angular Momentum Quantum Number
# V0  = 4.0                            # Strength of  the First Barrier
# D0  = 2.0                            # Thickness of the First Barrier
# x0  = [3.0]                           # Position of  the First Barrier
# V1  = 4.0                            # Strength of  the Second Barrier
# x1  = [6.0,7.0]                          # Position of  the Second Barrier
# D1  = 2.0                         # Thickness of the Second Barrier
# pk   = 4                             # Sharpness parameter v(x) = V0 exp( - ((x-x0)/D0)^(2pk) )
first = True                         # Controll Parameter for First Barrier
second = True                       # Controll Parameter for Second Barrier
   
calculation(x, w, N, n0, Z,l,V0, x0,D0,V1, x1, D1, pk, first, second)

print("Execution Complete")
