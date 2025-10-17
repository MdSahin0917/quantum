import numpy as np
from scipy.special import factorial, laguerre
#########################################################################################################
# Kinetic Energy Matrix Element (Ref[1,2])
#########################################################################################################
def kinetic(i,j, x, n, N):
    if i != j:
        coff = (-1)**(i - j  + 1)
        coff = coff * (x[j]/x[i])**n * np.sqrt(x[i]/x[j]**3) 
        coff = coff * ((2*n - 3)*x[j] - (2*n -1)*x[i])/(x[i] - x[j])**2
    else:
        coff = (-12*n**2 + 24*n -8 + (4*N + 2)*x[i] - x[i]**2 )/(12*x[i]**2)
    return coff


#########################################################################################################
# Potential (Double Barrier)
######################################################################################################### 
def vfunc(x, Z, l, V0, x0, D0, V1, x1, D1, pk, first, second):
    if first == True and second == False:
        c = -Z/x + l*(l+1.0)/(2*x**2) + V0*np.exp(-(x-x0)**(2*pk)/D0**(2*pk))
    elif first == True and second == True:
        c = -Z/x + l*(l+1.0)/(2*x**2) + V0*np.exp(-(x-x0)**(2*pk)/D0**(2*pk)) + V1*np.exp(-(x-x1)**(2*pk)/D1**(2*pk))
    else:
        c = -Z/x + l*(l+1.0)/(2*x**2)
    return c

#########################################################################################################
# Potential Energy Matrix Element (Ref[1,2])
######################################################################################################### 
def potential(i,j, xarray, Z, l,  V0, x0, D0, V1, x1, D1, p, first, second):
    if i != j:
        coff = 0.0
    else:
        coff = vfunc(xarray[i], Z, l,  V0, x0, D0, V1, x1, D1, p, first, second)  
    return coff



