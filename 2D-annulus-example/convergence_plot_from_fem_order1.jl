# Copyright (c) 2024 Ashwin. S. Nayak, Andrés Prieto, Daniel Fernández Comesaña
#
# This file is part of pressure-projections
#
# SPDX-License-Identifier:  MIT

"""
Computation and plot of the error curves obtained from the post-processing projection techniques
when the displacement field is computed from the first-order Raviart-Thomas FE discretizations
"""
# Load the package to use the generate_mesh function
include("meshing.jl")

# Load the computing.jl file to compute the relative errors
include("../computing.jl")

# Define the range of N values: logspace between 1e1 and 1e4
N_values = Integer.(round.(10.0.^range(log10(25), stop=log10(600), length=15)))
k=50.0 # wave number (2, 20, 100)
orderFE = 1 # order of the finite element
orderFE_postprocessing = 1 # order of the finite element for postprocessing

# Define the exact solution# Compute the analytical solution and the displacement error
p0 = 1.0 # Pressure amplitude
Zs = -1im * hankelh1(0, 2.0 * k) / hankelh1(1, 2.0 * k) # Impedance of the absorbing boundary
pex(x) = p0 * hankelh1(0, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k) # =-divergence(uex(x))
uex(x) = VectorValue(-x[1]/(sqrt(x[1]^2+x[2]^2)) / k * p0 * hankelh1(1, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k),
                        -x[2]/(sqrt(x[1]^2+x[2]^2)) / k * p0 * hankelh1(1, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k))
grad_pex(x) = VectorValue(-x[1]/(sqrt(x[1]^2+x[2]^2)) * k * p0 * hankelh1(1, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k),
                         -x[2]/(sqrt(x[1]^2+x[2]^2)) * k * p0 * hankelh1(1, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k)) # =k^2*uex(x)
exact_solution = Dict("ExactPressure"=>pex, "ExactDisplacement"=>uex, "ExactGradientPressure"=>grad_pex)

# Dictionaonary  with the boundary data associated to each boundary tag
uD(x) = VectorValue(0.0, 0.0)
pD(x) = p0
boundary_data = Dict("RigidBoundary"=>uD, "WallPressure"=>pD, "AbsorbingBoundary"=>Zs)

# Define the arrays to store the complex-valued relative errors
error_L2_u = zeros(size(N_values, 1))
error_H1_u = zeros(size(N_values, 1))
error_L2_p = zeros(size(N_values, 1))
error_H1_p = zeros(size(N_values, 1))
error_L2_proj_DG = zeros(size(N_values, 1))
error_H1_proj_DG = zeros(size(N_values, 1))
error_L2_proj_CG = zeros(size(N_values, 1))
error_H1_proj_CG = zeros(size(N_values, 1))
error_L2_proj_H1 = zeros(size(N_values, 1))
error_H1_proj_H1 = zeros(size(N_values, 1))
error_L2_proj_H2_DG = zeros(size(N_values, 1))
error_H1_proj_H2_DG = zeros(size(N_values, 1))

# Perform sweeping in h_values with a fixed k value
for j in 1:size(N_values, 1)
    # Compute the relative errors
    error_L2_u[j], error_H1_u[j], error_L2_p[j], error_H1_p[j], error_L2_proj_DG[j], error_H1_proj_DG[j], error_L2_proj_CG[j], error_H1_proj_CG[j], error_L2_proj_H1[j], error_H1_proj_H1[j], error_L2_proj_H2_DG[j], error_H1_proj_H2_DG[j] = compute_from_fem(k, N_values[j], orderFE, orderFE_postprocessing, exact_solution, boundary_data)
end
h_values = 1.0./N_values

# Save errors to a JDL2 file
using JLD2
save("convergence_data_from_fem_k=$(Integer(k))_orderFE=$(orderFE)_orderFEpostprocessing=$(orderFE_postprocessing).jld2",
    "h_values", h_values, "k", k, "orderFE", orderFE, "orderFE_postprocessing", orderFE_postprocessing,
    "error_L2_u", error_L2_u, "error_H1_u", error_H1_u, "error_L2_p", error_L2_p, "error_H1_p", error_H1_p,
    "error_L2_proj_DG", error_L2_proj_DG, "error_H1_proj_DG", error_H1_proj_DG, "error_L2_proj_CG", error_L2_proj_CG, "error_H1_proj_CG", error_H1_proj_CG,
    "error_L2_proj_H1", error_L2_proj_H1, "error_H1_proj_H1", error_H1_proj_H1, "error_L2_proj_H2_DG", error_L2_proj_H2_DG, "error_H1_proj_H2_DG", error_H1_proj_H2_DG)

# Plot the relative errors L2 and H1 with latex labels and log scale in Makie
# include("plot_fem_order1.jl")
