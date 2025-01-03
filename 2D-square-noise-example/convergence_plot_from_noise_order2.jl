# Copyright (c) 2024 Ashwin. S. Nayak, Andrés Prieto, Daniel Fernández Comesaña
#
# This file is part of pressure-projections
#
# SPDX-License-Identifier:  MIT

"""
Computation and plot of the error curves obtained from the post-processing projection techniques (second-order Raviart-Thomas FE discretizations)
when the displacement field is computed from a displacemtn field in the presence of noise
"""
# Load the computing.jl file to use the generate_mesh function
include("computing_with_noise.jl")

# Define the range of N values: logspace between 1e1 and 1e4
SNR_values = Integer.(round.(10.0.^range(log10(1), stop=log10(1e4), length=15)))
k=25.0 # wave number (2, 20, 100)
orderFE = 2 # order of the finite element used to interpolate the exact solution
orderFE_postprocessing = 2 # order of the finite element for postprocessing
PPW = 10 # Points per Wavelength in th mesh

# Define the exact solution# Compute the analytical solution and the displacement error
p0 = 1.0 # Pressure amplitude
pex(x) = p0 * exp(1im*k*x[2]) # =-divergence(uex(x))
uex(x) = VectorValue(0.0, 1im * p0 / k * exp(1im*k*x[2]))
grad_pex(x) = VectorValue(0.0, 1im * p0 * k * exp(1im*k*x[2])) # =k^2*uex(x)
exact_solution = Dict("ExactPressure"=>pex, "ExactDisplacement"=>uex, "ExactGradientPressure"=>grad_pex)

# Define the arrays to store the complex-valued relative errors
error_L2_u = zeros(size(SNR_values, 1))
error_H1_u = zeros(size(SNR_values, 1))
error_L2_proj_DG = zeros(size(SNR_values, 1))
error_H1_proj_DG = zeros(size(SNR_values, 1))
error_L2_proj_CG = zeros(size(SNR_values, 1))
error_H1_proj_CG = zeros(size(SNR_values, 1))
error_L2_proj_H1 = zeros(size(SNR_values, 1))
error_H1_proj_H1 = zeros(size(SNR_values, 1))
error_L2_proj_H2_DG = zeros(size(SNR_values, 1))
error_H1_proj_H2_DG = zeros(size(SNR_values, 1))

# Compute noise field (normalized at [0,1])
N, model, V_RT, noise_field, uh_interp = compute_noise(k, orderFE, PPW, exact_solution)

# Perform sweeping in SNR_values with a fixed k value
for j in 1:size(SNR_values, 1)
    # Compute the relative errors
    error_L2_u[j], error_H1_u[j], error_L2_proj_DG[j], error_H1_proj_DG[j], error_L2_proj_CG[j], error_H1_proj_CG[j], error_L2_proj_H1[j], error_H1_proj_H1[j], error_L2_proj_H2_DG[j], error_H1_proj_H2_DG[j] = compute_pressure_from_noise(k, N, SNR_values[j], orderFE_postprocessing, model, V_RT, noise_field, uh_interp, false)
end

# Save errors to a JDL2 file
using JLD2
save("results/convergence_data_from_noise_k=$(Integer(k))_orderFEpostprocessing=$(orderFE_postprocessing).jld2",
    "SNR_values", SNR_values, "k", k, "orderFE_postprocessing", orderFE_postprocessing,
    "error_L2_u", error_L2_u, "error_H1_u", error_H1_u,
    "error_L2_proj_DG", error_L2_proj_DG, "error_H1_proj_DG", error_H1_proj_DG, "error_L2_proj_CG", error_L2_proj_CG, "error_H1_proj_CG", error_H1_proj_CG,
    "error_L2_proj_H1", error_L2_proj_H1, "error_H1_proj_H1", error_H1_proj_H1, "error_L2_proj_H2_DG", error_L2_proj_H2_DG, "error_H1_proj_H2_DG", error_H1_proj_H2_DG)

# Plot the relative errors L2 and H1 with latex labels and log scale in Makie
# include("plot_noise_order2.jl")
