"""
Plot the error curves obtained from the post-processing projection techniques
when the displacement field is computed from the interpolation of the exact solution in first-order FE discretizations
"""

# Load variables from JDL2 file
using JLD2
file = "convergence_data_from_exact_k=50_orderFEpostprocessing=1.jld2"
vars = load(file)
k = vars["k"]
orderFE_postprocessing = vars["orderFE_postprocessing"] 
N_values = 1.0./vars["h_values"]
error_L2_u = vars["error_L2_u"]
error_H1_u = vars["error_H1_u"]
error_L2_p = vars["error_L2_p"]
error_H1_p = vars["error_H1_p"]
error_L2_proj_DG = vars["error_L2_proj_DG"]
error_H1_proj_DG = vars["error_H1_proj_DG"]
error_L2_proj_CG = vars["error_L2_proj_CG"]
error_H1_proj_CG = vars["error_H1_proj_CG"]
error_L2_proj_H1 = vars["error_L2_proj_H1"]
error_H1_proj_H1 = vars["error_H1_proj_H1"]
error_L2_proj_H2_DG = vars["error_L2_proj_H2_DG"]
error_H1_proj_H2_DG = vars["error_H1_proj_H2_DG"]

# Plot the relative errors L2 and H1 with latex labels and log scale in Makie
using CairoMakie, LaTeXStrings
fig = Figure(size = (600, 400), fonts = (; regular="CMU Serif")) ## probably you need to install this font in your system
ax = Axis(fig[1, 1], xlabel = L"$\lambda/h$", ylabel = L"Relative error (%)$$", xtickalign = 1,
xticksize = 10, ytickalign = 1, yticksize = 10, xscale=log10, yscale=log10, #xticks = vcat([2*pi/k.*N_values[end]], 10 .^range(-3, stop=-1, length=3)), yticks = 10 .^range(-5, stop=2, length=8),
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(9), 
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9))
# Error u
scatterlines!(2*pi/k.*N_values, error_L2_u, color = :black, linestyle = :dash ,  label = L"$\mathrm{L}^2$-error $I_{h,\mathrm{RT}_{0}}$", marker = '△', markersize = 10)
scatterlines!(2*pi/k.*N_values, error_H1_u, color = :black, label = L"$\mathrm{H(div)}$-error $I_{h,\mathrm{RT}_{0}}$", marker = '△', markersize = 10)
# Error p
scatterlines!(2*pi/k.*N_values, error_L2_p, color = :green, linestyle = :dash, label = L"$\mathrm{L}^2$-error $I_{h,\mathrm{CG}_{1}}$", marker = '◇', markersize = 10)
scatterlines!(2*pi/k.*N_values, error_H1_p, color = :green, label = L"$\mathrm{H}^1$-error $I_{h,\mathrm{CG}_{1}}$", marker = '◇', markersize = 10)
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
ylims!(ax, 3e-2, 2e2) 

# Plots triangles for order of convergence
order=2
x0=10.0^1.4; y0=2e-1
x1=10.0; y1=y0*(x0/x1)^order
lines!([(x0, y0), (x1, y0), (x1, y1), (x0, y0)], color=:black, linewidth=0.5)
text!([(10.0^((log10(x0)+log10(x1))/2), 10.0^(0.8*log10(y0)+0.2*log10(y1)))], text=L"O(h^2)", color=:black, fontsize=10)
# Triangle inverted
order=1
x0=10.0^1.4; y0=4e1
x1=10.0; y1=y0*(x1/x0)^order
lines!([(x0, y0), (x1, y0), (x0, y1), (x0, y0)], color=:black, linewidth=0.5)
text!([(10.0^(0.95*log10(x0)+0.05*log10(x1)), 10.0^(0.5*log10(y1)+0.5*log10(y0)))], text=L"O(h)", color=:black, fontsize=10)

# Save the figure to a PDF file
fig
save("convergence_plot_from_exact_k=$(Integer(k))_orderFEpostprocessing=$(orderFE_postprocessing).pdf", fig)