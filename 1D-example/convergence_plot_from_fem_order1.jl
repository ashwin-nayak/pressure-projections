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

# Define the range of N values: logspace between 1e1 and 500
N_values = Integer.(round.(10.0.^range(1, stop=log10(2000), length=15)))
k = 50.0 # wavenumber
orderFE = 1 # order of the finite element
orderFE_postprocessing = 1 # order of the finite element for postprocessing

# Define the exact solution# Compute the analytical solution and the displacement error
p0 = 1.0; # Pressure amplitude
Zs = 1.0; # Impedance of the absorbing boundary
pex(x) = p0*exp(1im*k*x[1]) # =-divergence(uex(x))
uex(x) = VectorValue(1im/k*p0*exp(1im*k*x[1]), 0.0)
grad_pex(x) = VectorValue(1im*k*p0*exp(1im*k*x[1]), 0.0) # =k^2*uex(x)
exact_solution = Dict("ExactPressure"=>pex, "ExactDisplacement"=>uex, "ExactGradientPressure"=>grad_pex)

# Dictionanary  with the boundary data associated to each boundary tag
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

# Perform sweeping in 2*pi/k.*N_values with a fixed k value
for j in 1:size(N_values, 1)
    # Compute the relative errors
    error_L2_u[j], error_H1_u[j], error_L2_p[j], error_H1_p[j], error_L2_proj_DG[j], error_H1_proj_DG[j], error_L2_proj_CG[j], error_H1_proj_CG[j], error_L2_proj_H1[j], error_H1_proj_H1[j], error_L2_proj_H2_DG[j], error_H1_proj_H2_DG[j] = compute_from_fem(k, N_values[j], orderFE, orderFE_postprocessing, exact_solution, boundary_data)
end

# Plot the relative errors L2 and H1 with latex labels and log scale in Makie
using CairoMakie, LaTeXStrings
fig = Figure(size = (600, 400), fonts = (; regular="CMU Serif")) ## probably you need to install this font in your system
ax = Axis(fig[1, 1], xlabel = L"$\lambda/h$", ylabel = L"Relative error (%)$$", xtickalign = 1,
xticksize = 10, ytickalign = 1, yticksize = 10, xscale=log10, yscale=log10, # xticks = vcat([2*pi/k.*N_values[end]], 10 .^range(-3, stop=-1, length=3)), yticks = 10 .^range(-5, stop=2, length=8)
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(9),
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9))
# Error u
scatterlines!(2*pi/k.*N_values, error_L2_u, color = :black, linestyle = :dash ,  label = L"$\mathrm{L}^2$-error $\mathbf{u}_{h,\mathrm{RT}_{0}}$", marker = '△', markersize = 10)
scatterlines!(2*pi/k.*N_values, error_H1_u, color = :black, label = L"$\mathrm{H(div)}$-error $\mathbf{u}_{h,\mathrm{RT}_{0}}$", marker = '△', markersize = 10)
# Error p
scatterlines!(2*pi/k.*N_values, error_L2_p, color = :green, linestyle = :dash, label = L"$\mathrm{L}^2$-error $p_{h,\mathrm{CG}_{1}}$", marker = '◇', markersize = 10)
scatterlines!(2*pi/k.*N_values, error_H1_p, color = :green, label = L"$\mathrm{H}^1$-error $p_{h,\mathrm{CG}_{1}}$", marker = '◇', markersize = 10)
# Error L2 DG projection
scatterlines!(2*pi/k.*N_values, error_L2_proj_DG, color = :orange, linestyle = :dash, label = L"$\mathrm{L}^2$-error $\Pi_{h,\mathrm{DG}_{0}}^{\mathrm{L}^2}$", marker = '▷', markersize = 10)
scatterlines!(2*pi/k.*N_values, error_H1_proj_DG, color = :orange, label = L"$\mathrm{H}^1$-error $\Pi_{h,\mathrm{DG}_{0}}^{\mathrm{L}^2}$", marker = '▷', markersize = 10)
# Error L2 CG projection
scatterlines!(2*pi/k.*N_values, error_L2_proj_CG, color = :blue, linestyle = :dash, label = L"$\mathrm{L}^2$-error $\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{L}^2}$", marker = '□', markersize = 10)
scatterlines!(2*pi/k.*N_values, error_H1_proj_CG, color = :blue, label = L"$\mathrm{H}^1$-error $\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{L}^2}$", marker = '□', markersize = 10)
# Error H1 CG projection
scatterlines!(2*pi/k.*N_values, error_L2_proj_H1, color = :red, linestyle = :dash, label = L"$\mathrm{L}^2$-error $\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{H}^1}$", marker = '○', markersize = 10)
scatterlines!(2*pi/k.*N_values, error_H1_proj_H1, color = :red, label = L"$\mathrm{H}^1$-error $\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{H}^1}$", marker = '○', markersize = 10)
# Error H2 DG projection
scatterlines!(2*pi/k.*N_values, error_L2_proj_H2_DG, color = :purple, linestyle = :dash, label = L"$\mathrm{L}^2$-error $\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{H}^1}$", marker = '☆', markersize = 10)
scatterlines!(2*pi/k.*N_values, error_H1_proj_H2_DG, color = :purple, label = L"$\mathrm{H}^1$-error $\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{H}^1}$", marker = '☆', markersize = 10)
# axislegend(; position = :rb, nbanks = 1, framecolor = (:grey, 0.5))
Legend(fig[1, 2], ax, halign = :right, valign = :top, labelsize = 17, markersize = 10, markerstrokewidth = 0.5)
xlims!(ax, 2*pi/k.*N_values[end], 2*pi/k.*N_values[1])
ylims!(ax, 5e-3, 2e2)

# Plots triangles for order of convergence
order=2
x0=10.0^1.5; y0=1e-1
x1=10.0; y1=y0*(x0/x1)^order
lines!([(x0, y0), (x1, y0), (x1, y1), (x0, y0)], color=:black, linewidth=0.5)
text!([(10.0^((log10(x0)+log10(x1))/2), 10.0^(0.8*log10(y0)+0.2*log10(y1)))], text=L"O(h^2)", color=:black, fontsize=10)
# Triangle inverted
order=1
x0=1e2; y0=1e1
x1=10.0^1.5; y1=y0*(x1/x0)^order
lines!([(x0, y0), (x1, y0), (x0, y1), (x0, y0)], color=:black, linewidth=0.5)
text!([(10.0^(0.95*log10(x0)+0.05*log10(x1)), 10.0^(0.5*log10(y1)+0.5*log10(y0)))], text=L"O(h)", color=:black, fontsize=10)

fig
save("convergence_plot_from_fem_k=$(Integer(k))_orderFE=$(orderFE)_orderFEpostprocessing=$(orderFE_postprocessing).pdf", fig)
