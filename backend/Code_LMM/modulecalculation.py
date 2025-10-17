import numpy as np
from scipy.special import factorial
from scipy.special import roots_laguerre, laguerre
from scipy.linalg  import eigh
from scipy.integrate import simpson
from numba import njit, prange
from moduledfftfunc import phi, phiarray
from moduleLMM import kinetic, potential, vfunc
import os
from datetime import datetime
Now=datetime.now()

def calculation(x, w, N, n0,  Z, l, V0, x0,D0,V1, x1, D1, pk, first, second):
    n=1
    hamiltonian = np.zeros((N,N))


    ###################################################
    # Matrix Calculation
    ###################################################
    for i in range(0,N):
        for j in range(0,N):
            hamiltonian[i,j] = kinetic(i,j, x, n, N)/2.0 + potential(i,j, x, Z, l,  V0, x0, D0, V1, x1, D1, pk, first, second)

    ######################################################################################################
    # Property Calculation
    ######################################################################################################

    e, v = eigh(hamiltonian)           # Diagonalisation
    eg = e[n0-l-1] #* 27.21138             # Ground State Energy
    norm = np.sum(v.T[n0-l-1]**2)          # Normalisation
    xexp = np.sum(x * v.T[n0-l-1]**2)      # <x> expectation value of x
    x2exp = np.sum(x**2 * v.T[n0-l-1]**2)  # <x^2> expectation value of x^2
    x3exp = np.sum(x**3 * v.T[n0-l-1]**2)  # <x^3> expectation value of x^3
    xexp_1 = np.sum(v.T[n0-l-1]**2/x)      # <1/x> expectation value of 1/x
    xexp_2 = np.sum(v.T[n0-l-1]**2/x**2)   # <1/x^2> expectation value of 1/x^2
    xexp_3 = np.sum(v.T[n0-l-1]**2/x**2)   # <1/x^3> expectation value of 1/x^3


    CPE = -Z * xexp_1
    LPE =  l*(l+1.0) * xexp_2/2.0
    PE  = CPE + LPE
    KE =  eg - PE
    sd = np.sqrt(np.abs(x2exp - xexp**2))
    pf = sd/xexp
    alphaK = 4.0 * x2exp **2/9.0
    alphaB = 2.0 * (6.0 * x2exp **3 + 3.0 * x3exp **2 - 8.0 * xexp * x2exp * x3exp)
    alphaB = alphaB / (3.0 * (9.0 * x2exp - 8.0 * xexp**2)) 
    #####################################################################################################
    # Lagrange Function
    #####################################################################################################    
    def f(xj, i):
        if xj != x[i]:
            c = 0.0
        else:
            c = 1/np.sqrt(lamda[i])
        return c
    #####################################################################################################
    # Psi reconstruction with Lagrange function
    #####################################################################################################     
    # Psi reconstruction 
    def psi(xvar): 
        s = 0.0
        for i in range(0,N):
            s += v.T[n0-l-1][i] * f(xvar, i)
        return s
    slmdat = [2.53102424696929,
              2.09907862496785,
              2.04112500613393,
              2.02065962276832,
              2.01053680740952,
              2.00457769907142,
              2.00067684953873,
              1.99793460613024,
              1.99590577775833,
              1.99434607120380]

    # Shannon entropy
    lamda = [w[i]*np.exp(x[i]) for i in range(len(x)) if x[i] <600.0]
    Sr = -np.sum([lamda[i] * psi(x[i]) **2 * np.log(psi(x[i]) **2/x[i]**2) for i in range(len(lamda)) if psi(x[i]) **2/x[i]**2 > 1.0e-16]) + slmdat[l]

    # Disequilibrium
    Dr = np.sum([lamda[i] * psi(x[i]) **4/x[i]**2 for i in range(len(lamda))])
    
    # LMC complexity
    Cr = np.exp(Sr) * Dr

    #####################################################################################################
    #Continous Lagrange Function
    #####################################################################################################   
    def f_full(x, i , n, N, xarray):
        coeff = (-1)**i * np.sqrt(xarray[i])
        regul = (x/xarray[i])**n
        body  = laguerre(N)(x) / (x - xarray[i])
        weig  = np.exp(-x/2.0)
        return coeff * regul * body * weig
    #####################################################################################################
    #Psi with Continous Lagrange Function
    #####################################################################################################   
    def psi_full(xvar): 
        s = 0.0
        for i in range(0,N):
            s += v.T[n0-l-1][i] * f_full(xvar, i, n, N, x)
        return s

    xvar = np.linspace(1e-4, 100, 20000)
    dx = xvar[1] - xvar[0]
    psiarray = psi_full(xvar)
    FWHM0 = 2*D0*(np.log(2)**(1/(2*pk))) 
    FWHM1 = 2*D1*(np.log(2)**(1/(2*pk)))
    #####################################################################################################
    #Trapping, squeezing, swelling probability
    #####################################################################################################   
    if first == True and second == False:
        xi = x0 - FWHM0/2.0 
        xf = x0 + FWHM0/2.0

        # squeezing probability
        xR      = np.array([i for i in xvar if i <= xi]) # x < xi 
        psiR    = np.array([psiarray[i] for i in range(0, len(xvar)) if xvar[i] <= xi] ) 
        squeeze = simpson(psiR**2, x = xR)

        # Tunneling probability
        xR      = np.array([i for i in xvar if i > xi and i <= xf])  # xi < x < xf 
        psiR    = np.array([psiarray[i] for i in range(0, len(xvar)) if xvar[i] > xi and xvar[i] <= xf ] ) 
        tunnel = simpson(psiR**2, x = xR)


        # swelling probability
        xT      = np.array([i for i in xvar if i > xf])
        psiT    = np.array([psiarray[i] for i in range(0, len(xvar)) if xvar[i] > xf] ) 
        swelling = simpson(psiT**2, x = xT)
         
        
    elif first == True and second == True:
        x0i = x0 - FWHM0/2.0 
        x0f = x0 + FWHM0/2.0
        x1i = x1 - FWHM0/2.0 
        x1f = x1 + FWHM0/2.0

        # squeezing probability (x < x0i)
        xR      = np.array([i for i in xvar if i <= x0i])
        psiR    = np.array([psiarray[i] for i in range(0, len(xvar)) if xvar[i] <= x0i] ) 
        squeeze = simpson(psiR**2, x = xR)

        # swelling probability    (x > x1f)
        xT      = np.array([i for i in xvar if i > x1f])
        psiT    = np.array([psiarray[i] for i in range(0, len(xvar)) if xvar[i] > x1f] )
        swelling = simpson(psiT**2, x = xT)

        # tunnel probability - 1st (x0i < x < x0f)
        xTrap      = np.array([i for i in xvar if i > x0i and i <= x0f])
        psiTrap    = np.array([psiarray[i] for i in range(0, len(xvar)) if xvar[i] > x0i and xvar[i] <= x0f] ) 
        Tunnel1 = simpson(psiTrap**2, x = xTrap)

        # tunnel probability - 2nd  (x1i < x < x1f)
        xTrap      = np.array([i for i in xvar if i > x1i and i <= x1f])
        psiTrap    = np.array([psiarray[i] for i in range(0, len(xvar)) if xvar[i] > x1i and xvar[i] <= x1f] ) 
        Tunnel2 = simpson(psiTrap**2, x = xTrap)

        # trapping probability  (x0f < x < x1i)
#        xTrap      = np.array([i for i in xvar if i > x0f and i <= x1i])
#        psiTrap    = np.array([psiarray[i] for i in range(0, len(xvar)) if xvar[i] > x0f and xvar[i] <= x1i] ) 
#        Trap = simpson(psiTrap**2, x = xTrap)
        
        # Flatten arrays (ensures 1D)
        psiTrap = np.ravel(psiTrap)
        xTrap = np.ravel(xTrap)

	# Check lengths
        if len(psiTrap) != len(xTrap):
            print("Length mismatch:", len(psiTrap), len(xTrap))
        else:
            Trap = simpson(psiTrap**2, x=xTrap)

    else:
        trap = 0.0



    
    
 

    label={ 
        '0':'s',
        '1':'p',
        '2':'d',
        '3':'f',
        '4':'g',
    }
    element={ 
        '1':'Hydrogen',
        '2':'Helium',
        '3':'Lithium',
        '4':'Beryllium',
        '5':'Boron',
        '6':'Carbon',
        '7':'Nitrogen',
        '8':'Oxygen',
        '9':'Fluorine',
        '10':'Neon',
        '11':'Sodium',
        '12':'Potassium'


    }

    output = "output"

    if not os.path.exists(output):
        os.mkdir(output)
    
  
    cd_folder_name = output +'/'+f'{element[str(Z)]} in {n0}{label[str(l)]}'
    
    if not os.path.exists(f'{cd_folder_name}'):
        os.mkdir(cd_folder_name)

    if first == True and second == False:
        folder_name = f'({x0}_{V0}_{D0})'
    elif first == True and second == True:
        folder_name = f'({x0}_{V0}_{D0}), ({x1}_{V1}_{D1})'
    else:
        folder_name = 'free'
    
    full_name   = cd_folder_name +'/' +folder_name
    if not os.path.exists(full_name):
        os.mkdir(full_name)
     
    #
    f = open(f'{cd_folder_name}/{folder_name}/info.txt', 'w')
    f.write(f'''TITLE: Comprehensive study of the one-electron system ({element[str(Z)]}) in {n0}{label[str(l)]} state
        under penetrable barrier-type potential \n''')
    f.write('*************************************************************************************************** \n\n')
    f.write('AUTHOR                                       : Koustav Das Chakladar, Dr. Santanu Mondal, Dr. Jayanta Kumar Saha \n')
    f.write('DATE AND TIME                                : {n}\n'.format(n=Now))
    f.write('*************************************************************************************************** \n')
    f.write('SYSTEM SPECIFICATION (ATOMIC UNIT)           : \n\n')
    f.write('*************************************************************************************************** \n')
    f.write('ELEMENT                                      :{e}\n'.format(e=element[str(Z)]))
    f.write('ATOMIC NUMBER                                :{m} a.u.\n'.format(m=Z))
    f.write('PRINCIPLE QUANTUM NUMBER                     :{n}\n'.format(n=n0))
    f.write('ANGULAR MOMENTUM QUANTUM NUMBER (L)          :{L}\n'.format(L=l))
    f.write(f'STATE SPECIFICATION                         :{n0}{label[str(l)]}\n')
    if first == True and second == False:
        f.write(f'FIRST BARRIER HEIGHT                            :{V0} a.u.\n')
        f.write(f'POSITION OF THE FIRST BARRIER                   :{x0} a.u.\n')
        f.write(f'FIRST BARRIER FWHM                              :{FWHM0:.12e} a.u.\n')
    elif first == True and second == True:
        f.write(f'FIRST BARRIER HEIGHT                            :{V0} a.u.\n')
        f.write(f'POSITION OF THE FIRST BARRIER                   :{x0} a.u.\n')
        f.write(f'FIRST BARRIER FWHM                              :{FWHM0:.12e} a.u.\n')
        f.write(f'SECOND BARRIER HEIGHT                            :{V1} a.u.\n')
        f.write(f'POSITION OF THE SECOND BARRIER                   :{x1} a.u.\n')
        f.write(f'SECOND BARRIER FWHM                              :{FWHM1:.12e} a.u.\n')
    f.write('*************************************************************************************************** \n')
    f.write('STRUCTURAL PROPERTIES USING OPTIMIZED WAVE FUNCTION  (ATOMIC UNIT)   \n\n')
    f.write('*************************************************************************************************** \n')
    f.write('OVERLAP                                       :  {e:.8e} \n'.format(e=norm))
    f.write('ENERGY                                        :  {e:.8e} \n'.format(e=eg))
    f.write('KINETIC ENERGY                                :  {e:.8e}\n'.format(e=KE))
    f.write('POTENTIAL ENERGY                              : {e:.8e}\n'.format(e=PE))
    f.write('COULOMBIC POTENTIAL ENERGY                    : {e:.8e}\n'.format(e=CPE))
    f.write('CENTRIFUGAL POTENTIAL ENERGY                  : {e:.8e}\n'.format(e=-LPE))
    f.write('EXPECTATION VALUE OF r                        :  {e:.8e}\n'.format(e=xexp))
    f.write('EXPECTATION VALUE OF r^2                      :  {e:.8e}\n'.format(e=x2exp))
    f.write('EXPECTATION VALUE OF r^3                      :  {e:.8e}\n'.format(e=x3exp))
    f.write('EXPECTATION VALUE OF 1/r                      :  {e:.8e}\n'.format(e=xexp_1))
    f.write('EXPECTATION VALUE OF 1/r^2                    :  {e:.8e}\n'.format(e=xexp_2))
    f.write('EXPECTATION VALUE OF 1/r^3                    :  {e:.8e}\n'.format(e=xexp_3))
    f.write('SD(r) =sqrt(<r^2> - <r>^2)                    :  {e:.8e}\n'.format(e=sd))
    f.write('PEARSON CORRELATION OF r  = SD(r)/<r>         :  {e:.8e}\n'.format(e=pf))
    f.write('KIRKWOOD POLARISABILITY                       :  {e:.8e}\n'.format(e=alphaK))
    f.write('BUCKINGHUM POLARISABILITY                     :  {e:.8e}\n'.format(e=alphaB))

    if first == True and second == False:
        f.write('SQUEEZING PROBABILITY                         :  {e:.8e}\n'.format(e=squeeze))
        f.write('TUNNELING (1ST) PROBABILITY                   :  {e:.8e}\n'.format(e=tunnel))
        f.write('SWELLING  PROBABILITY                         :  {e:.8e}\n'.format(e=swelling))

         
        
    elif first == True and second == True:
        f.write('SQUEEZING PROBABILITY                         :  {e:.8e}\n'.format(e=squeeze))
        f.write('SWELLING  PROBABILITY                         :  {e:.8e}\n'.format(e=swelling))
        f.write('TRAPPING  PROBABILITY                         :  {e:.8e}\n'.format(e=Trap))
        f.write('TUNNELING (1ST) PROBABILITY                   :  {e:.8e}\n'.format(e=Tunnel1))
        f.write('TUNNELING (2ND) PROBABILITY                   :  {e:.8e}\n'.format(e=Tunnel2))

    else:
        trap = 0.0
    f.write('*************************************************************************************************** \n')
    f.write('INFORMATION THEORETIC MEASURES  (ATOMIC UNIT)   \n\n')
    f.write('*************************************************************************************************** \n')
    f.write('SHANON ENTROPY                                :  {e:.8e}\n'.format(e=Sr))
    f.write('DISEQUILIBRIUM                                :  {e:.8e}\n'.format(e=Dr))
    f.write('LMC COMPLEXITY                                :  {e:.8e}\n'.format(e=Cr))

    f.close()
    
    den = psiarray**2
    ybarr = np.array([vfunc(i, Z, l, V0, x0, D0, V1, x1, D1, pk, first, second) for i in xvar[1:]])


    f = open(f'{full_name}/densityPS.dat','w')
    for i in range(0,len(xvar)):
        f.write(f'{xvar[i]:.12e}     {den[i]:.12e}\n')
    f.close()
    
    f = open(f'{full_name}/barrier.dat','w')
    for i in range(0,len(xvar[1:])):
        f.write(f'{xvar[i]:.12e}     {ybarr[i]:.12e}\n')
    f.close()

    
    
    
