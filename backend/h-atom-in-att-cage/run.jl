##############################
# Run
if length(ARGS) < 5
    error("Usage: julia run.jl N l n v0")
end
# get arguments from Python
v0 = parse(Float64, ARGS[1])
n = parse(Int, ARGS[2])
l = parse(Int, ARGS[3])
N = parse(Int, ARGS[4])
maxR = parse(Int, ARGS[5])
println("Running julia with v0 = ", v0)


###################################
include("Hlikefree.jl")
include("calculation.jl")
using .Hlikefree
using .calculation
using Dates
using NumericalIntegration
using LaTeXStrings
using Printf

##############################
#ALL PARAMETERS
###################################
R= [0.1+0.1*i for i in 0:maxR]
Delta = 5.0

z = 1

a0 = 0.01 
aN = 500.0
g = (aN/a0)^(1/N)
arr = [a0 * g^i for i in 0:N-1]


# angular shannon entropy
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

#println(R)
for i in R
    i = round(i, digits=5)
    cal(i, Delta, v0, z, n, l, N, a0, aN, g, arr, slmdat)
end





