"""
Plotting the performance curves (CPU time and memory) obtained from the post-processing projection techniques
when the displacement field is computed from the first-order Raviart-Thomas FE discretizations
"""
# Plot from fem order 1
orderFE = 1

# Load variables from JDL2 file
using JLD2
file = "CPUtime_data_from_fem_k=50_orderFE=1_orderFEpostprocessing=1.jld2"
vars = load(file)
k = vars["k"]
orderFE_postprocessing = vars["orderFE_postprocessing"]
N_values = 1.0./vars["h_values"]
# CPUtimes from ns to s
CPUtime_u = vars["CPUtime_u"]/1e9
CPUtime_p = vars["CPUtime_p"]/1e9
CPUtime_proj_DG = vars["CPUtime_proj_DG"]/1e9
CPUtime_proj_CG = vars["CPUtime_proj_CG"]/1e9
CPUtime_proj_H1 = vars["CPUtime_proj_H1"]/1e9
CPUtime_proj_H2_DG = vars["CPUtime_proj_H2_DG"]/1e9

# Check if std is larger than mean values and truncate it to mean values
CPUtime_u[findall(CPUtime_u[:,2].>CPUtime_u[:,1]),2] .= CPUtime_u[findall(CPUtime_u[:,2].>CPUtime_u[:,1]),1]*0.99
CPUtime_p[findall(CPUtime_p[:,2].>CPUtime_p[:,1]),2] .= CPUtime_p[findall(CPUtime_p[:,2].>CPUtime_p[:,1]),1]*0.99
CPUtime_proj_DG[findall(CPUtime_proj_DG[:,2].>CPUtime_proj_DG[:,1]),2] .= CPUtime_proj_DG[findall(CPUtime_proj_DG[:,2].>CPUtime_proj_DG[:,1]),1]*0.99
CPUtime_proj_CG[findall(CPUtime_proj_CG[:,2].>CPUtime_proj_CG[:,1]),2] .= CPUtime_proj_CG[findall(CPUtime_proj_CG[:,2].>CPUtime_proj_CG[:,1]),1]*0.99
CPUtime_proj_H1[findall(CPUtime_proj_H1[:,2].>CPUtime_proj_H1[:,1]),2] .= CPUtime_proj_H1[findall(CPUtime_proj_H1[:,2].>CPUtime_proj_H1[:,1]),1]*0.99
CPUtime_proj_H2_DG[findall(CPUtime_proj_H2_DG[:,2].>CPUtime_proj_H2_DG[:,1]),2] .= CPUtime_proj_H2_DG[findall(CPUtime_proj_H2_DG[:,2].>CPUtime_proj_H2_DG[:,1]),1]*0.99

# Check if sts is NaN and then set it to 0.01*mean values
CPUtime_u[findall(isnan.(CPUtime_u[:,2])),2] .= CPUtime_u[findall(isnan.(CPUtime_u[:,2])),1]*0.01
CPUtime_p[findall(isnan.(CPUtime_p[:,2])),2] .= CPUtime_p[findall(isnan.(CPUtime_p[:,2])),1]*0.01
CPUtime_proj_DG[findall(isnan.(CPUtime_proj_DG[:,2])),2] .= CPUtime_proj_DG[findall(isnan.(CPUtime_proj_DG[:,2])),1]*0.01
CPUtime_proj_CG[findall(isnan.(CPUtime_proj_CG[:,2])),2] .= CPUtime_proj_CG[findall(isnan.(CPUtime_proj_CG[:,2])),1]*0.01
CPUtime_proj_H1[findall(isnan.(CPUtime_proj_H1[:,2])),2] .= CPUtime_proj_H1[findall(isnan.(CPUtime_proj_H1[:,2])),1]*0.01
CPUtime_proj_H2_DG[findall(isnan.(CPUtime_proj_H2_DG[:,2])),2] .= CPUtime_proj_H2_DG[findall(isnan.(CPUtime_proj_H2_DG[:,2])),1]*0.01

# Plot CPU times with latex labels and log scale in Makie
using CairoMakie, LaTeXStrings
fig = Figure(size = (600, 400), fonts = (; regular="CMU Serif")) ## probably you need to install this font in your system
ax = Axis(fig[1, 1], xlabel = L"$\lambda/h$", ylabel = L"CPU time (s)$$", xtickalign = 1,
xticksize = 10, ytickalign = 1, yticksize = 10, xscale=log10, yscale=log10, # xticks = vcat([2*pi/k.*N_values[end]], 10 .^range(-3, stop=-1, length=3)), yticks = 10 .^range(-5, stop=2, length=8)
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(9),
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9))
# Error u
#scatterlines!(2*pi/k.*N_values, CPUtime_u[:,1], color = :black, label = L"$\mathbf{u}_{h,\mathrm{RT}_{0}}$", marker = '△', markersize = 10)
#errorbars!(2*pi/k.*N_values, CPUtime_u[:,1], CPUtime_u[:,2], color = :black, whiskerwidth = 10)
#band!(2*pi/k.*N_values, CPUtime_u[:,1]+CPUtime_u[:,2], CPUtime_u[:,1]-CPUtime_u[:,2], color = (:black, 0.5))
# Error p
#scatterlines!(2*pi/k.*N_values, CPUtime_p[:,1], color = :green, label = L"$p_{h,\mathrm{CG}_{1}}$", marker = '◇', markersize = 10)
#errorbars!(2*pi/k.*N_values, CPUtime_p[:,1], CPUtime_p[:,2], color = :green, whiskerwidth = 10)
#band!(2*pi/k.*N_values, CPUtime_p[:,1]+CPUtime_p[:,2], CPUtime_p[:,1]-CPUtime_p[:,2], color = (:green, 0.5))
# Error L2 DG projection
scatterlines!(2*pi/k.*N_values, CPUtime_proj_DG[:,1], color = :orange, label = L"$\Pi_{h,\mathrm{DG}_{0}}^{\mathrm{L}^2}$", marker = '▷', markersize = 10)
#errorbars!(2*pi/k.*N_values, CPUtime_proj_DG[:,1], CPUtime_proj_DG[:,2], color = :orange, whiskerwidth = 10)
band!(2*pi/k.*N_values, CPUtime_proj_DG[:,1]+CPUtime_proj_DG[:,2], CPUtime_proj_DG[:,1]-CPUtime_proj_DG[:,2], color = (:orange, 0.5))
# Error L2 CG projection
scatterlines!(2*pi/k.*N_values, CPUtime_proj_CG[:,1], color = :blue, label = L"$\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{L}^2}$", marker = '□', markersize = 10)
#errorbars!(2*pi/k.*N_values, CPUtime_proj_CG[:,1], CPUtime_proj_CG[:,2], color = :blue, whiskerwidth = 10)
band!(2*pi/k.*N_values, CPUtime_proj_CG[:,1]+CPUtime_proj_CG[:,2], CPUtime_proj_CG[:,1]-CPUtime_proj_CG[:,2], color = (:blue, 0.5))
# Error H1 CG projection
scatterlines!(2*pi/k.*N_values, CPUtime_proj_H1[:,1], color = :red, label = L"$\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{H}^1}$", marker = '○', markersize = 10)
#errorbars!(2*pi/k.*N_values, CPUtime_proj_H1[:,1], CPUtime_proj_H1[:,2], color = :red, whiskerwidth = 10)
band!(2*pi/k.*N_values, CPUtime_proj_H1[:,1]+CPUtime_proj_H1[:,2], CPUtime_proj_H1[:,1]-CPUtime_proj_H1[:,2], color = (:red, 0.5))
# Error H2 DG projection
scatterlines!(2*pi/k.*N_values, CPUtime_proj_H2_DG[:,1], color = :purple, label = L"$\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{H}^1}$", marker = '☆', markersize = 10)
#errorbars!(2*pi/k.*N_values, CPUtime_proj_H2_DG[:,1], CPUtime_proj_H2_DG[:,2], color = :purple, whiskerwidth = 10)
band!(2*pi/k.*N_values, CPUtime_proj_H2_DG[:,1]+CPUtime_proj_H2_DG[:,2], CPUtime_proj_H2_DG[:,1]-CPUtime_proj_H2_DG[:,2], color = (:purple, 0.5))
# axislegend(; position = :rb, nbanks = 1, framecolor = (:grey, 0.5))
Legend(fig[1, 2], ax, halign = :right, valign = :top, labelsize = 17, markersize = 10, markerstrokewidth = 0.5)
xlims!(ax, 2*pi/k.*N_values[end], 2*pi/k.*N_values[1])
ylims!(ax, 4e-4, 2e1)

# Plots triangles for order of convergence
order=2
x0=10.0^1.3; y0=1e-3
x1=10.0^0.8; y1=y0*(x0/x1)^order
lines!([(x0, y0), (x1, y0), (x0, y1), (x0, y0)], color=:black, linewidth=0.5)
text!([(10.0^(0.8*log10(x0)+0.2*log10(x1)), 10.0^(0.8*log10(y0)+0.2*log10(y1)))], text=L"O(h^2)", color=:black, fontsize=10)

fig
save("CPUtime_plot_from_fem_k=$(Integer(k))_orderFE=$(orderFE)_orderFEpostprocessing=$(orderFE_postprocessing).pdf", fig)

# Plot memory and allocations with latex labels and log scale in Makie
Memory_u = vars["Memory_u"]
Memory_p = vars["Memory_p"]
Memory_proj_DG = vars["Memory_proj_DG"]
Memory_proj_CG = vars["Memory_proj_CG"]
Memory_proj_H1 = vars["Memory_proj_H1"]
Memory_proj_H2_DG = vars["Memory_proj_H2_DG"]

# Plot Memory with latex labels and log scale in Makie with values in the first column of the data
fig = Figure(size = (600, 400), fonts = (; regular="CMU Serif")) ## probably you need to install this font in your system
ax = Axis(fig[1, 1], xlabel = L"$\lambda/h$", ylabel = L"Memory (Mb)$$", xtickalign = 1,
xticksize = 10, ytickalign = 1, yticksize = 10, xscale=log10, yscale=log10, # xticks = vcat([2*pi/k.*N_values[end]], 10 .^range(-3, stop=-1, length=3)), yticks = 10 .^range(-5, stop=2, length=8)
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(9),
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9))
# Memmory for u in MegaBytes
#scatterlines!(2*pi/k.*N_values, Memory_u[:,1]/1e6, color = :black, label = L"$\mathbf{u}_{h,\mathrm{RT}_{0}}$", marker = '△', markersize = 10)
#scatterlines!(2*pi/k.*N_values, Memory_p[:,1]/1e6, color = :green, label = L"$p_{h,\mathrm{CG}_{1}}$", marker = '◇', markersize = 10)
scatterlines!(2*pi/k.*N_values, Memory_proj_DG[:,1]/1e6, color = :orange, label = L"$\Pi_{h,\mathrm{DG}_{0}}^{\mathrm{L}^2}$", marker = '▷', markersize = 10)
scatterlines!(2*pi/k.*N_values, Memory_proj_CG[:,1]/1e6, color = :blue, label = L"$\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{L}^2}$", marker = '□', markersize = 10)
scatterlines!(2*pi/k.*N_values, Memory_proj_H1[:,1]/1e6, color = :red, label = L"$\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{H}^1}$", marker = '○', markersize = 10)
scatterlines!(2*pi/k.*N_values, Memory_proj_H2_DG[:,1]/1e6, color = :purple, label = L"$\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{H}^1}$", marker = '☆', markersize = 10)
# axislegend(; position = :rb, nbanks = 1, framecolor = (:grey, 0.5))
Legend(fig[1, 2], ax, halign = :right, valign = :top, labelsize = 17, markersize = 10, markerstrokewidth = 0.5)
xlims!(ax, 2*pi/k.*N_values[end], 2*pi/k.*N_values[1])
ylims!(ax, 2e0, 7e3)

# Plots triangles for order of convergence
order=2
x0=10.0^1.5; y0=4e0
x1=10.0^1.0; y1=y0*(x0/x1)^order
lines!([(x0, y0), (x1, y0), (x0, y1), (x0, y0)], color=:black, linewidth=0.5)
text!([(10.0^(0.8*log10(x0)+0.2*log10(x1)), 10.0^(0.8*log10(y0)+0.2*log10(y1)))], text=L"O(h^2)", color=:black, fontsize=10)

fig
save("Memory_plot_from_fem_k=$(Integer(k))_orderFE=$(orderFE)_orderFEpostprocessing=$(orderFE_postprocessing).pdf", fig)

# Plot memory and allocations with latex labels and log scale in Makie
Memory_u = vars["Memory_u"]
Memory_p = vars["Memory_p"]
Memory_proj_DG = vars["Memory_proj_DG"]
Memory_proj_CG = vars["Memory_proj_CG"]
Memory_proj_H1 = vars["Memory_proj_H1"]
Memory_proj_H2_DG = vars["Memory_proj_H2_DG"]

# Plot Allocations with latex labels and log scale in Makie with values in the first column of the data
fig = Figure(size = (600, 400), fonts = (; regular="CMU Serif")) ## probably you need to install this font in your system
ax = Axis(fig[1, 1], xlabel = L"$\lambda/h$", ylabel = L"Allocations$$", xtickalign = 1,
xticksize = 10, ytickalign = 1, yticksize = 10, xscale=log10, yscale=log10, # xticks = vcat([2*pi/k.*N_values[end]], 10 .^range(-3, stop=-1, length=3)), yticks = 10 .^range(-5, stop=2, length=8)
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(9),
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9))
# Memmory for u in MegaBytes
#scatterlines!(2*pi/k.*N_values, Memory_u[:,2], color = :black, label = L"$\mathbf{u}_{h,\mathrm{RT}_{0}}$", marker = '△', markersize = 10)
#scatterlines!(2*pi/k.*N_values, Memory_p[:,2], color = :green, label = L"$p_{h,\mathrm{CG}_{1}}$", marker = '◇', markersize = 10)
scatterlines!(2*pi/k.*N_values, Memory_proj_DG[:,2], color = :orange, label = L"$\Pi_{h,\mathrm{DG}_{0}}^{\mathrm{L}^2}$", marker = '▷', markersize = 10)
scatterlines!(2*pi/k.*N_values, Memory_proj_CG[:,2], color = :blue, label = L"$\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{L}^2}$", marker = '□', markersize = 10)
scatterlines!(2*pi/k.*N_values, Memory_proj_H1[:,2], color = :red, label = L"$\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{H}^1}$", marker = '○', markersize = 10)
scatterlines!(2*pi/k.*N_values, Memory_proj_H2_DG[:,2], color = :purple, label = L"$\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{H}^1}$", marker = '☆', markersize = 10)
# axislegend(; position = :rb, nbanks = 1, framecolor = (:grey, 0.5))
Legend(fig[1, 2], ax, halign = :right, valign = :top, labelsize = 17, markersize = 10, markerstrokewidth = 0.5)
xlims!(ax, 2*pi/k.*N_values[end], 2*pi/k.*N_values[1])

# # Plots triangles for order of convergence
# order=2
# x0=10.0^1.5; y0=3e0
# x1=10.0^1.0; y1=y0*(x0/x1)^order
# lines!([(x0, y0), (x1, y0), (x0, y1), (x0, y0)], color=:black, linewidth=0.5)
# text!([(10.0^(0.8*log10(x0)+0.2*log10(x1)), 10.0^(0.8*log10(y0)+0.2*log10(y1)))], text=L"O(h^2)", color=:black, fontsize=10)

fig
save("Allocations_plot_from_fem_k=$(Integer(k))_orderFE=$(orderFE)_orderFEpostprocessing=$(orderFE_postprocessing).pdf", fig)
