

module calculation
export cal
include("Hlikefree.jl")
using .Hlikefree
using Dates
using NumericalIntegration
using LaTeXStrings
using Printf
function cal(R, Delta, v0, z, n, l, N, a0, aN, g, arr, slmdat)

    Dict1 = Dict( 0 => "s", 1 => "p", 2 => "d", 3 => "f", 4 => "g", 5 => "h")
    namedir = string(n) * Dict1[l]
    data_file_name = "$namedir/(V = $v0, R = $R , Delta = $Delta)/Correlated_data_(V = $v0, R = $R , Delta = $Delta).txt"
    ev, C = energy(N, l, z, v0, R, Delta, arr)

    function vpot(r, v0, R, Delta, l)
        if r >=R && r< R+Delta
            return v0 - z/r + l*(l+1)/ (2*r^2)
        else
            return - z/r + l*(l+1)/ (2*r^2)
        end
    end

    ###################################
    # READING ROOTS AND WEIGHTS
    ###################################
    f=open("backend/h-atom-in-att-cage/Gauss_Leg/mesh_200.txt","r")
    # Declaring empty array to store the data files from mominfo
    roots  =  []
    weights = []

    # reading each lines in input file
    for line in eachline(f)
        splt =  split(line, r"\s+")    
        
        # Evaluating the data by parsing as string
        d1 = eval(Meta.parse(splt[1])) 
        d2 = eval(Meta.parse(splt[2]))
        # putting the data from each lines in the declared empty array
        push!(roots,d1)
        push!(weights, d2)
    end
    #####################################################################
    # Calculation of Radial moments, kinetic energy, potential energy
    #####################################################################
    S = zeros(N,N)
    Ek = zeros(N,N)
    Ep = zeros(N,N)

    # matrix generation for kinetic energy, potential energy and radial moments
    for i::Int in 1:N
        for j::Int in 1:N
            S[i, j]      = overlap(i, j, l, arr)
            Ek[i, j]     = kinetic(i, j, l, arr)
            Ep[i, j]     = potential(i, j, l, z, v0, R, Delta, arr)
        end
    end

    # Linear Variational Parameters 
    v1 = real(C[1:N,n-l])
     


    # normalization constant
    Nc = v1' * S * v1

    # overlap check
    o = (v1'* S * v1)/ Nc


    # Kinetic energy
    kin = (v1'* Ek * v1)/ Nc


    # Potential energy
    pot = (v1'* Ep * v1)/ Nc

        

    # Continuous r mesh
    r = [0.000001 + 0.001*i for i in 0:300000]         

    limr = length(r)
    # Continuous Position space wavefunction (For Plottng)
    psiarrcont = [psi(i, l, arr, v1,N) for i in r]  

    # Gauss-Legendre Mesh in Position space
    rmesh = (ones(length(roots)) .+ roots) ./ (ones(length(roots)) .- roots)

    # Position space wavefunction on Gauss-Legendre Mesh
    psiarrmesh = [psi(i, l, arr, v1,N) for i in rmesh]


    dencont_r = psiarrcont .^2 
    denmesh_r = psiarrmesh .^2

    # Normalization 
    Norm = Nquad(denmesh_r,roots,weights)
    dencont_r = psiarrcont .^2 ./ Norm
    denmesh_r = psiarrmesh .^2 ./ Norm


    # Energy and Radial Moments
    en      = (ev[n-l])
    # rexpt   = (v1'* rmat * v1)/ Nc
    # r2expt  = (v1'* r2mat * v1)/ Nc
    # r3expt  = (v1'* r3mat * v1)/ Nc
    # r1nexpt = (v1'* r1nmat * v1)/ Nc
    # r2nexpt = (v1'* r2nmat * v1)/ Nc



    rexpt   =  sum(dencont_r .* r)       *0.001                  #Radmom(1, denmesh_r, roots, weights)
    r2expt  =  sum(dencont_r .* r .^2 )  *0.001                  #Radmom(2, denmesh_r, roots, weights)
    r3expt  =  sum(dencont_r .* r .^3 )  *0.001                  #Radmom(3, denmesh_r, roots, weights)
    r1nexpt =  sum(dencont_r ./ r     )  *0.001                  #Radmom(-1, denmesh_r, roots, weights)
    r2nexpt =  sum(dencont_r ./ r.^2  )  *0.001                  #Radmom(-2, denmesh_r, roots, weights)
    r3nexpt =  sum(dencont_r ./ r.^3  )  *0.001                  #Radmom(-3, denmesh_r, roots, weights)
    Sr      = Shannon(denmesh_r,roots,weights) + slmdat[l+1]
    Dr      = Disequilibrium(denmesh_r,roots,weights)

    swprob  = integrate([r[i] for i in 1:limr-1 if r[i] < R], [dencont_r[i] for i in 1:limr-1 if r[i] < R])
    tunprob = integrate([r[i] for i in 1:limr-1 if r[i] >= R && r[i] < R+Delta], [dencont_r[i] for i in 1:limr-1 if r[i] >= R && r[i] < R+Delta])
    swlprob = integrate([r[i] for i in 1:limr-1 if r[i] >= R + Delta],[dencont_r[i] for i in 1:limr-1 if r[i] >= R + Delta])

    totalprob = swprob + tunprob + swlprob

    varr = [vpot(i, v0, R, Delta, l) for i in r]
    vpotarr = varr .* dencont_r

    swprobvarr  = 0.001 * sum([vpotarr[i] for i in 1:limr-1 if r[i] < R])
    tunprobvarr = 0.001 * sum([vpotarr[i] for i in 1:limr-1 if r[i] >= R && r[i] < R+Delta])
    swlprobvarr = 0.001 * sum([vpotarr[i] for i in 1:limr-1 if r[i] >= R + Delta])

    totalprobvarr = swprobvarr + tunprobvarr + swlprobvarr

    cen = l*(l+1) * r2nexpt / 2.0
    std = sqrt(r2expt - rexpt^2)
    pr  = std / rexpt
    alphaK = 4.0 * r2expt ^2 /9.0
    alphaB = 2.0 * (6.0 * r2expt^3 + 3.0 * r3expt^2 - 8.0 * rexpt * r2expt * r3expt)
    alphaB = alphaB / (3.0 * (9.0 * r2expt - 8.0 * rexpt^2)) 



    #####################################################################
    # MOMENTUM SPACE PROPERTIES
    #####################################################################

    p = [0.0001 + 0.01*i for i in 0:20000]
    phiarrcont = [phi(i, l, arr, v1,N) for i in p]
    pmesh = (ones(length(roots)) .+ roots) ./ (ones(length(roots)) .- roots)
    phiarrmesh = [phi(i, l, arr, v1,N) for i in pmesh]
    dencont_p = phiarrcont .^2 .* p .^ 2
    denmesh_p = phiarrmesh .^2 .* pmesh .^ 2
    Norm = Nquad(denmesh_p,roots,weights)

    dencont_p = dencont_p ./ Norm
    denmesh_p = denmesh_p ./ Norm



    pexpt   =  Radmom(1, denmesh_p, roots, weights)
    p2expt  =  Radmom(2, denmesh_p, roots, weights)
    p3expt  =  Radmom(3, denmesh_p, roots, weights)
    p4expt  =  Radmom(4, denmesh_p, roots, weights)
    p1nexpt =  Radmom(-1, denmesh_p, roots, weights)
    p2nexpt =  Radmom(-2, denmesh_p, roots, weights)
    Sp      = Shannon(denmesh_p,roots,weights) + slmdat[l+1]
    Dp      = Disequilibrium(denmesh_p,roots,weights)

    Ip      =  4 * r2expt
    Ir      =  4 * p2expt

    St = Sp + Sr
    Iprod = Ip * Ir
    prodrp = r2expt * p2expt
    Ih = sqrt(prodrp)
    shape_r =  Dr * exp(Sr) 
    shape_p =  Dp * exp(Sp) 
    shape_t =  Dr * Dp * exp(St) 

    Jr = (1/(2 * pi * exp(1))) * exp(2.0 * Sr/3.0)
    Jp = (1/(2 * pi * exp(1))) * exp(2.0 * Sp/3.0)
    Jt = (1/(2 * pi * exp(1))) * exp(2.0 * St/3.0)

    Fr = Ir * Jr
    Fp = Ip * Jp
    Ft = Fr * Fp
    Dt = Dr * Dp
    
    Cr = Dr * Jr
    Cp = Dp * Jp
    Ct = Dr * Dp * Jt


    if !isdir(namedir)
        mkdir(namedir)
    end

    if !isdir("$namedir/(V = $v0, R = $R , Delta = $Delta)")
        mkdir("$namedir/(V = $v0, R = $R , Delta = $Delta)")
    end


    heading = "TITLE: Comprehensive Study of H-like Atom Under Penetrable Repulsive Barrier Potential  \n"
    Day = now()
    Date_and_time = Dates.format(Day , "yyyy-mm-dd HH:MM:SS")
    # File writing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    f = open(data_file_name, "w")
    write(f, "***************************************************************************************************\n")
    write(f, heading)
    write(f, "***************************************************************************************************\n")
    write(f, "AUTHOR                                   : Mr. Koustav D. Chakladar, Dr. Jayanta K. Saha \n")
    write(f, "DATE AND TIME                            : $Date_and_time \n\n")
    write(f, "***************************************************************************************************\n")
    write(f, "ALL QUANTITIES ARE IN ATOMIC UNIT IF NOT SPECIFIED OTHERWISE \n")
    write(f, "*************************************************************************************************** \n")
    write(f, "ATOMIC NUMBER               (Z)                 : $z \n")
    write(f, "ANGULAR MOMENTUM            (l)                 : $l \n")
    write(f, "NUMBER OF BASIS             (N)                 : $N \n")
    write(f, "INITIAL VALUE OF α          (α0)                : $a0 \n")
    write(f, "FINAL VALUE OF α            (αn)                : $aN \n")
    write(f, "POTENTIAL  OF THE BARRIER   (V)                 : $v0  \n")
    write(f, "POSITION   OF THE BARRIER   (R)                 : $R \n")
    write(f, "THICKNESS  OF THE BARRIER   (Δ)                 : $Delta \n")
    write(f, "*************************************************************************************************** \n")
    write(f , "OVERLAP                                :     $o \n")
    write(f , "ENERGY                                 :     $en \n")
    write(f , "KINETIC ENERGY                         :     $kin \n")
    write(f , "POTENTIAL ENERGY                       :     $pot \n")
    write(f , "CENTRIFUGAL KINETIC ENERGY             :     $cen \n")
    write(f, "EXPECTATION VALUE OF r                  :     $rexpt  \n")
    write(f, "EXPECTATION VALUE OF r^2                :     $r2expt \n")
    write(f, "EXPECTATION VALUE OF r^3                :     $r3expt \n")
    write(f, "EXPECTATION VALUE OF 1/r                :     $r1nexpt \n")
    write(f, "EXPECTATION VALUE OF 1/r^2              :     $r2nexpt \n")
    write(f, "EXPECTATION VALUE OF 1/r^3              :     $r3nexpt \n")
    write(f, "SD(r) =sqrt(<r^2> - <r>^2)              :     $std  \n")
    write(f,"PEARSON CORRELATION OF r =std(r)/<r>     :     $pr  \n")
    write(f, "EXPECTATION VALUE OF p                  :     $pexpt  \n")
    write(f, "EXPECTATION VALUE OF p^2                :     $p2expt \n")
    write(f, "EXPECTATION VALUE OF p^3                :     $p3expt \n")
    write(f, "EXPECTATION VALUE OF p^4                :     $p4expt \n")
    write(f, "EXPECTATION VALUE OF 1/p                :     $p1nexpt \n")
    write(f, "EXPECTATION VALUE OF 1/p^2              :     $p2nexpt \n")
    write(f, "KIRKWOOD POLARISABILITY                 :     $alphaK \n")
    write(f, "BUCKINGHUM POLARISABILITY               :     $alphaB \n")
    
    write(f,"SHANNON ENTROPY (POSITION SPACE)         :     $Sr \n")
    write(f,"SHANNON ENTROPY (MOMENTUM SPACE)         :     $Sp \n")
    write(f,"SHANNON ENTROPY (CONJUGATE SPACE)        :     $St \n")



    write(f,"FISHER ENTROPY  (POSITION SPACE)         :     $Ir \n")
    write(f,"FISHER ENTROPY  (MOMENTUM SPACE)         :     $Ip \n")
    write(f,"FISHER PRODUCT  (CONJUGATE SPACE)        :     $Iprod \n")



    write(f,"DISEQUILIBRIUM  (POSITION SPACE)         :     $Dr  \n")
    write(f,"DISEQUILIBRIUM  (MOMENTUM SPACE)         :     $Dp  \n")
    write(f,"DISEQUILIBRIUM PRODUCT (CONJUGATE SPACE) :     $Dt  \n")

    
    write(f,"SHAPE COMPLEXITY (POSITION SPACE)        :     $Cr \n")
    write(f,"SHAPE COMPLEXITY (MOMENTUM SPACE)        :     $Cp \n")
    write(f,"SHAPE COMPLEXITY (CONJUGATE SPACE)       :     $Ct \n")


    write(f,"LMC COMPLEXITY   (POSITION SPACE)        :     $shape_r \n")
    write(f,"LMC COMPLEXITY   (MOMENTUM SPACE)        :     $shape_p \n")
    write(f,"LMC COMPLEXITY   (CONJUGATE SPACE)       :     $shape_t \n")

    
    write(f,"EXPONENTIAL SHANNON (POSITION SPACE)     :     $Jr  \n")
    write(f,"EXPONENTIAL SHANNON (MOMENTUM SPACE)     :     $Jp  \n")
    write(f,"EXPONENTIAL SHANNON (CONJUGATE SPACE)    :     $Jt  \n")
    
   
    write(f,"FISHER-SHANNON (POSITION SPACE)          :     $Fr  \n")
    write(f,"FISHER-SHANNON (MOMENTUM SPACE)          :     $Fp  \n")
    write(f,"FISHER-SHANNON (CONJUGATE SPACE)         :     $Ft  \n")
    
    write(f,"HEISENBERG INFORMATION (CONJUGATE SPACE) :     $Ih  \n")

    write(f,"SQUEEZING  PROBABILITY                   :     $swprob  \n")
    write(f,"TUNNELLING PROBABILITY                   :     $tunprob  \n")
    write(f,"SWELLING PROBABILITY                     :     $swlprob  \n")
    write(f,"TOTAL  PROBABILITY                       :     $totalprob  \n")

    write(f,"REGION-I STRENGTH                        :     $swprobvarr  \n")
    write(f,"REGION-II STRENGTH                       :     $tunprobvarr  \n")
    write(f,"REGION-III STRENGTH                      :     $swlprobvarr  \n")
    write(f,"TOTAL POTENTIAL ENERGY                   :     $totalprobvarr  \n")
 
    close(f )                                                                                             



    #########################################################################
    # Plots of Density and Compton profile
    #########################################################################

      # Continuous r mesh
    power = [-4 + 0.001 * i for i in 0:10000 ]
    r = 10 .^ power



    # Continuous Position space wavefunction (For Plottng)
    psiarrcont = [psi(i, l, arr, v1,N) for i in r]  


    Norm = Nquad(denmesh_r,roots,weights)
    dencont_r = psiarrcont .^2 / Norm
    # position space density
    open("$namedir/(V = $v0, R = $R , Delta = $Delta)/one_particle_density_position_(V = $v0, R = $R , Delta = $Delta).dat", "w") do file
        lim = length(dencont_r)
        for i in 1:lim
            pii, dii = r[i], dencont_r[i]
            kprint = @sprintf("%0.10f                     %0.10E\n", pii, dii)
            write(file,kprint)
        end
    end

    # momentum space density

    open("$namedir/(V = $v0, R = $R , Delta = $Delta)/one_particle_density_momentum_(V = $v0, R = $R , Delta = $Delta).dat", "w") do file
        lim = length(dencont_p)
        for i in 1:lim
            pii, dii = p[i], dencont_p[i]
            kprint = @sprintf("%0.10f                     %0.10E\n", pii, dii)
            write(file,kprint)
        end
    end



    f = open("$namedir/(V = $v0, R = $R , Delta = $Delta)/compton_(V = $v0, R = $R , Delta = $Delta).dat", "w")
    # compton profile
    for i in 1:3000
        q = p[i]
        lim = length(p)
        pcomparr = p[i:lim]
        dencomparr  = dencont_p[i:lim] ./ pcomparr
        comp  = integrate(pcomparr,dencomparr, SimpsonEven())
        kprint = @sprintf("%0.10f                     %0.10E\n", q, comp)
        write(f,kprint)
    end
    close(f)

    open("$namedir/(V = $v0, R = $R , Delta = $Delta)/Barrier_(V = $v0, R = $R , Delta = $Delta).dat", "w") do file
        lim = length(r)
        for i in 1:lim
            pii, dii = r[i], vpot(r[i], v0, R, Delta, l)
            kprint = @sprintf("%0.10f                     %0.10E\n", pii, dii)
            write(file,kprint)
        end
    end
end
end
