"""
Compute the solution of displacement-based using first-order Raviart-Thomas discretization
with the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry

# Load special function package (only required for the 2D simulations)
using SpecialFunctions

# Load the plotting packages
using GridapMakie, CairoMakie, FileIO, LaTeXStrings

# Load benchmarking package
using BenchmarkTools

# Load the meshing.jl file to use the generate_mesh function
include("projecting.jl")

function compute_from_fem_with_CPUtime(k, N, orderFE, orderFE_postprocessing,  boundary_data, plot_flag=false)

    # Generate the mesh for a given mesh size
    model = generate_mesh(N)

    # Compute the solution of displacement-based formulation using first-order Raviart-Thomas discretization
    uh = RT_solve(k, orderFE, model, boundary_data)
    benchmark_u = @benchmark uh = RT_solve($k, $orderFE, $model, $boundary_data)

    # Compute the solution of pressure-based formulation using Lagrange discretization
    benchmark_p = @benchmark ph = CG_solve($k, $orderFE_postprocessing, $model, $boundary_data) 

    # Compute the pressure projections
    benchmark_proj_DG = @benchmark ph_L2_DG = L2_DG_projection($uh, $orderFE_postprocessing, $model)
    benchmark_proj_CG = @benchmark ph_L2_CG = L2_CG_projection($uh, $orderFE_postprocessing, $model)
    benchmark_proj_H1 = @benchmark ph_H1 = H1_projection($k, $uh, $orderFE_postprocessing, $model)
    benchmark_proj_H2_DG = @benchmark ph_H2 = H2_projection_DG($k, $uh, $orderFE_postprocessing, $model, $N)


    # Plot using GridapMakie some fields only if plot_flag is true
    if plot_flag==true
        Ω = Triangulation(model) # Computational domain
        # Plot mesh
        fig, ax, plt = wireframe(Ω, color=:black, linewidth=1.0)
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        # Axis equal
        ax.aspect = DataAspect()
        # Save the figure to pdf format with 300 dpi
        save("./results/mesh_N=$(N).pdf", fig)

        # Plot pressure field
        fig , ax , plt = plot(Ω, real(ph_H2))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        # cbar = Colorbar(fig[2,1], plt, label=L"\mathrm{Re}(p_{\mathrm{ex}})", vertical=false)
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{H}^1})")
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_p_H2_N=$(N).pdf", fig)
    end
    
    # # Write the results to a VTK file with fileneame "results_N.vtu" with N the mesh size
    # writevtk(Ω,"./results/results_N=$(N)_order=$(orderFE).vtu", cellfields=["Re(uh)"=>real(uh), "Im(uh)"=>imag(uh),
    #                                         "Re(uex)"=>real(uex ∘ xp), "Im(uex)"=>imag(uex ∘ xp),
    #                                         "Re(ph)"=>real(ph), "Im(ph)"=>imag(ph),
    #                                         "Re(pex)"=>real(pex ∘ xp), "Im(pex)"=>imag(pex ∘ xp)])
    return benchmark_u, benchmark_p, benchmark_proj_DG, benchmark_proj_CG, benchmark_proj_H1, benchmark_proj_H2_DG
end
