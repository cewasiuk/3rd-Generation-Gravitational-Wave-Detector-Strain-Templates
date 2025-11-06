using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("Annihilation.jl")

z = 1 + 2im;

@test isapprox(expinti(z), 1.0421677081649356 + im * 3.7015014259378742)