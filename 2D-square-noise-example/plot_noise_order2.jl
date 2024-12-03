# Copyright (c) 2024 Ashwin. S. Nayak, Andrés Prieto, Daniel Fernández Comesaña
#
# This file is part of pressure-projections
#
# SPDX-License-Identifier:  MIT

"""
Plot the error curves obtained from the post-processing projection techniques (second-order Raviart-Thomas FE discretizations)
when the displacement field is computed from a displacemtn field in the presence of noise
"""
# Plot from fem order 1
k=25.0
orderFE_postprocessing = 2

# Load variables from JDL2 file
using JLD2
file = "results/convergence_data_from_noise_k=$(Integer(k))_orderFEpostprocessing=$(orderFE_postprocessing).jld2"
vars = load(file)
SNR_values = vars["SNR_values"]
error_L2_u = vars["error_L2_u"]
error_H1_u = vars["error_H1_u"]
error_L2_proj_DG = vars["error_L2_proj_DG"]
error_H1_proj_DG = vars["error_H1_proj_DG"]
error_L2_proj_CG = vars["error_L2_proj_CG"]
error_H1_proj_CG = vars["error_H1_proj_CG"]
error_L2_proj_H1 = vars["error_L2_proj_H1"]
error_H1_proj_H1 = vars["error_H1_proj_H1"]
error_L2_proj_H2_DG= vars["error_L2_proj_H2_DG"]
error_H1_proj_H2_DG = vars["error_H1_proj_H2_DG"]

# Plot the relative errors L2 and H1 with latex labels and log scale in Makie
using CairoMakie, LaTeXStrings
fig = Figure(size = (600, 400), fonts = (; regular="CMU Serif")) ## probably you need to install this font in your system
ax = Axis(fig[1, 1], xlabel = L"SNR$$", ylabel = L"Relative error (%)$$", xtickalign = 1,
xticksize = 10, ytickalign = 1, yticksize = 10, xscale=log10, yscale=log10, xticks = 10 .^range(-1, stop=4, length=6), yticks = 10 .^range(-3, stop=3, length=7),
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(9),
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9))
# Error u
scatterlines!(SNR_values, error_L2_u, color = :black, linestyle = :dash ,  label = L"$\mathrm{L}^2$-error $\mathbf{u}_{\mathrm{SNR}}$", marker = '△', markersize = 10)
scatterlines!(SNR_values, error_H1_u, color = :black, label = L"$\mathrm{H(div)}$-error $\mathbf{u}_{\mathrm{SNR}}$", marker = '△', markersize = 10)
# Error L2 DG projection
scatterlines!(SNR_values, error_L2_proj_DG, color = :orange, linestyle = :dash, label = L"$\mathrm{L}^2$-error $\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{L}^2}$", marker = '▷', markersize = 10)
scatterlines!(SNR_values, error_H1_proj_DG, color = :orange, label = L"$\mathrm{H}^1$-error $\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{L}^2}$", marker = '▷', markersize = 10)
# Error L2 CG projection
scatterlines!(SNR_values, error_L2_proj_CG, color = :blue, linestyle = :dash, label = L"$\mathrm{L}^2$-error $\Pi_{h,\mathrm{CG}_{2}}^{\mathrm{L}^2}$", marker = '□', markersize = 10)
scatterlines!(SNR_values, error_H1_proj_CG, color = :blue, label = L"$\mathrm{H}^1$-error $\Pi_{h,\mathrm{CG}_{2}}^{\mathrm{L}^2}$", marker = '□', markersize = 10)
# Error H1 CG projection
scatterlines!(SNR_values, error_L2_proj_H1, color = :red, linestyle = :dash, label = L"$\mathrm{L}^2$-error $\Pi_{h,\mathrm{CG}_{2}}^{\mathrm{H}^1}$", marker = '○', markersize = 10)
scatterlines!(SNR_values, error_H1_proj_H1, color = :red, label = L"$\mathrm{H}^1$-error $\Pi_{h,\mathrm{CG}_{2}}^{\mathrm{H}^1}$", marker = '○', markersize = 10)
# Error H2 DG projection
scatterlines!(SNR_values, error_L2_proj_H2_DG, color = :purple, linestyle = :dash, label = L"$\mathrm{L}^2$-error $\Pi_{h,\mathrm{DG}_{2}}^{\mathrm{H}^2}$", marker = '☆', markersize = 10)
scatterlines!(SNR_values, error_H1_proj_H2_DG, color = :purple, label = L"$\mathrm{H}^1$-error $\Pi_{h,\mathrm{DG}_{2}}^{\mathrm{H}^2}$", marker = '☆', markersize = 10)
# axislegend(; position = :rb, nbanks = 1, framecolor = (:grey, 0.5))
Legend(fig[1, 2], ax, halign = :right, valign = :top, labelsize = 17, markersize = 10, markerstrokewidth = 0.5)
xlims!(ax, SNR_values[1], SNR_values[end])
ylims!(ax, 0.05, 5e2)

# Plots triangles for order of convergence, oriented to the right
order=1
x0=5.0; y0=0.3
x1=40.0; y1=y0*(x1/x0)^order
lines!([(x0, y0), (x1, y0), (x0, y1), (x0, y0)], color=:black, linewidth=0.5)
text!([(10.0^((0.95*log10(x0)+0.05*log10(x1))), 10.0^(0.95*log10(y0)+0.05*log10(y1)))], text=L"O(\mathrm{SNR})", color=:black, fontsize=10)

fig
save("results/convergence_plot_from_noise_k=$(Integer(k))_orderFEpostprocessing=$(orderFE_postprocessing).pdf", fig)
