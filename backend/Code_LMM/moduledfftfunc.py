import numpy as np
from numba import njit, prange

@njit(parallel=True)
def dffunc(x,y,p_arr):
    result = np.zeros_like(p_arr, dtype=np.complex128)
    for ip in prange(len(p_arr)):
        pii = p_arr[ip]
        yp = 0.0
        for i in prange(len(x)):
            yp += y[i] *  np.exp(-1j*pii*x[i])
        result[ip] = yp
    return result / np.sqrt(2*np.pi)

@njit
def jl(r, p, l):
    if l == 0:
        return np.sin(p*r) / (p*r)
    elif l == 1:
        return (- np.cos(p*r)/(p*r)  + np.sin(p*r) / (p*r)**2 )
    elif l == 2:
        return (- np.sin(p*r)/(p*r)  - 3*np.cos(p*r) / (p*r)**2 + 3*np.sin(p*r) / (p*r)**3 )
    elif l == 3:
        return ( np.cos(p*r)/(p*r)  - 6*np.sin(p*r) / (p*r)**2 - 15*np.cos(p*r) / (p*r)**3 + 15*np.sin(p*r) / (p*r)**4)
    elif l == 4:
        return ( np.sin(p*r)/(p*r)  + 10*np.cos(p*r) / (p*r)**2 - 45*np.sin(p*r) / (p*r)**3 - 105*np.cos(p*r) / (p*r)**4 + 105*np.sin(p*r) / (p*r)**5)


        
        
@njit(parallel= True)
def phi(p, psiarray, xarray, l):
    results = np.zeros(len(xarray))
    for i in prange(0, len(xarray)): 
        xi = xarray[i]
        yi = psiarray[i]
        results[i] = np.sqrt(2/np.pi) * yi * xi * jl(xi, p, l)
 
    dx =  xarray[1] - xarray[0]
    value = dx/3 *( (results[0] + results[-1]) + (4 * np.sum(results[1:-1:2])) + (2 * np.sum(results[2:-2:2]) ) )
    return value


@njit(parallel= True)
def phiarray(parray, psiarray, xarray, l):
    results = np.zeros(len(parray))
    for i in prange(len(parray)):
        pii = parray[i]
        results[i] = phi(pii, psiarray, xarray, l)
    return results
    
if __name__=='__main__':
    x = np.linspace(0,100,10000)
    y = np.sin(x) + np.sin(3*x) + np.sin(2*x)
    p_arr = np.linspace(0,10,1000)
    yp = dffunc(x,y,p_arr)
    plt.plot(p_arr,np.absolute(yp))
    plt.xlim(0,2)
    plt.savefig('dftjpg.jpg')
    plt.show()


"""
@njit(parallel= True)
def phi(p, psiarray, xarray):
    body = 0.0
    for i in prange(0, len(psiarray)): 
        xi = xarray[i]
        yi = psiarray[i]
        body += yi * np.sin(p*xi)
    return np.sqrt(2/np.pi) * body / p


@njit(parallel= True)
def phiarray(parray, psiarray, xarray):
    results = np.zeros(len(parray))
    for i in prange(len(parray)):
        pii = parray[i]
        results[i] = phi(pii, psiarray, xarray)
    return results
"""
