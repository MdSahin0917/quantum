
module Hlikefree
export A, Acage, Ashell
export Amom, overlap, kinetic, potential, radialmom
export energy, psi, phi
export Nquad, Radmom, Shannon, Disequilibrium
#########################################################################################################
# Basis Integral (0, infinity)
#########################################################################################################
function A(l::Int, a::Float64)::Float64
    coeff::Float64 = factorial(big(l))/a^(l+1)
    return coeff
end


#########################################################################################################
# Basis Integral (0, R)
#########################################################################################################
function Acage(l::Int, a::Float64, R::Float64)::Float64
    coeff::Float64 = factorial(big(l))/a^(l+1)
    sum::Float64 = 0
    for i::Int in 0:l
        sum += (a*R)^i/factorial(big(i))
    end
    return (1-sum*exp(-a*R))*coeff
end
##############################################################################################################
# Basis Integral (R1, R2)
##############################################################################################################
function Ashell(l::Int, a::Float64, R1::Float64, R2::Float64)::Float64
    A0::Float64 = (exp(-a*R1)- exp(-a*R2))/a
    A::Float64  = 0.0
    if l == 0
        return A0
    else
        for i::Int in range(1,l)
            A  = ( R1^i*exp(-a*R1)- R2^i*exp(-a*R2) )/a + i*A0/a
            A0 = A
        end
        return A
    end
end

##############################################################################################################
# Fourier Basis Integral (R1, R2)
##############################################################################################################
function Amom(p, l, a, q)
    P = sqrt(p^2 + a^2)
    xi = atan(p/a)
    term1  = 0.0
    for k in 0:q
        coeff  = 2 * factorial(q+k) * factorial(l - k -1) / ( (2*p)^(k+1) * factorial(k) * factorial(q-k) * P^(l-k) )
        body   = cos( pi*(k - q -1)/2.0  + (l - k)*xi  )
        term1 += coeff * body
    end
    return term1
end




#####################################################################################################
# Overlap
#####################################################################################################
function overlap(i, j, l, arr)
    ai, aj = arr[i], arr[j]
    s = A(2*l +2, ai + aj)
    return s
end

#####################################################################################################
# Kinetic energy
#####################################################################################################
function kinetic(i, j, l, arr)
    ai, aj = arr[i], arr[j]
    t1 = (aj^2 + ai^2)* A(2*l +2, ai + aj)
    t2 = -2 *(aj+ai)*(l+1) * A(2*l+1, ai + aj)
    t3 =  2* l*(l+1) * A(2*l, ai + aj)
    return -(t1+t2+t3)/4.0
end

#####################################################################################################
# potential energy
#####################################################################################################
function potential(i, j, l, z, v0, R, Delta, arr)
    ai, aj = arr[i], arr[j]
    col = -z * A(2*l +1, ai + aj)
    cen = (l*(l+1)/2.0) * A(2*l, ai + aj) 
    barr = v0 * Ashell(2*l+2, ai+aj, R, R+Delta)
    return col + cen + barr
end

#####################################################################################################
# RADIAL MOMENTS
#####################################################################################################
function radialmom(n, i, j, l, arr)
    ai, aj = arr[i], arr[j]
    s = A(2*l +2 +  n, ai + aj)
    return s
end

#####################################################################################################
# Energy
#####################################################################################################
using LinearAlgebra
function energy(N,l, z, v0, R, Delta, arr)
    S::Matrix{Float64} = zeros(N,N)
    T::Matrix{Float64} = zeros(N,N)
    V::Matrix{Float64} = zeros(N,N)
    
    for i::Int in 1:N
        for j::Int in 1:N
            S[i, j] = overlap(i, j, l, arr)
            T[i, j] = kinetic(i, j, l, arr)
            V[i, j] = potential(i, j, l, z, v0, R, Delta, arr)
        end    
    end
    result = eigen(T+V, S)
    ev = result.values
    C  = result.vectors
    return (real(ev),real(C))
end

#####################################################################################################
# POSITION SPACE WAVEFUNCTION
#####################################################################################################

function psi(r, l, arr, C, N)
    summ = 0.0
    for i in 1:N
        ai = arr[i]
        ci = C[i]
        summ += ci * exp(-ai * r)
    end
    return summ * r^(l+1)
end
#####################################################################################################
# MOMENTUM SPACE WAVEFUNCTION
#####################################################################################################
function phi(p, l, arr, C, N)
    summ = 0.0
    for i in 1:N
        ai = arr[i]
        ci = C[i]
        summ += ci * Amom(p, l+2, ai, l)
    end
    return summ * sqrt(2/pi)
end

#####################################################################################################
# NORMALIZATION
#####################################################################################################

function Nquad(denmom, roots, weights)
    summ = 0.0
    for i in range(1, length(roots))
        xi = roots[i]
        # pp = ((1 + xi)/(1 - xi))
        lambdai =  weights[i]
        summ += lambdai * (2/(1-xi)^2) * denmom[i]
    end
    return summ
end

#####################################################################################################
# RADIAL MOMENTS
#####################################################################################################


function Radmom(n, denmom, roots, weights)
    summ = 0.0
    for i in range(1, length(roots))
        xi = roots[i]
        pp = ((1 + xi)/(1 - xi))
        lambdai =  weights[i]
        summ += lambdai * (2/(1-xi) ^ 2) * denmom[i] * pp^n
    end
    return summ
end

#####################################################################################################
# INFORMATION ENTROPY
#####################################################################################################

function Shannon(denmom, roots, weights)
    summ = 0.0
    for i in range(1, length(roots))
        xi = roots[i]
        pp = ((1 + xi)/(1 - xi))
        lambdai =  weights[i]
        if denmom[i]/pp ^ 2 > 1e-12
            summ += lambdai * (2/(1-xi)^2) * denmom[i] * log(denmom[i]/pp^2)
        end
    end
    return -summ
end

#####################################################################################################
# DISEQUILIBRIUM
#####################################################################################################

function Disequilibrium(denmom, roots, weights)
    summ = 0.0
    for i in range(1, length(roots))
        xi = roots[i]
        pp = ((1 + xi)/(1 - xi))
        lambdai =  weights[i]
        if denmom[i]/pp ^ 2 > 1e-12
            summ += lambdai * (2/(1-xi)^2) * denmom[i] * (denmom[i]/pp^2)
        end
    end
    return summ 
end     
end 