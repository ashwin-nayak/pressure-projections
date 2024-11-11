"""
Computation of the different post-processing pressure projections solution and measurement of their errors
from a displacement field, which is either obtained from a Raviart-Thomas discretization (implemented in Gridap)
or from piecewise polynomial interpolations on a triangular mesh
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry

# Load special function package (only required for the 2D simulations)
using SpecialFunctions

# Load the plotting packages
using GridapMakie, CairoMakie, FileIO, LaTeXStrings

# Load the meshing.jl file to use the generate_mesh function
include("projecting.jl")

# Compute the L^2-relative error in the entire computational domain
function L2_error(fh, fex, dΩ, xp)
    error = fh - fex
    return 100 * sqrt(abs(sum(∫(error ⋅ error)*dΩ)/sum(∫((fex ∘ xp) ⋅ (fex ∘ xp))*dΩ)))
end

# Compute the L^2-relative error in the entire computational domain
function L2_error_vector(fh, fex, dΩ, xp)
    error = fh - fex
    return 100 * sqrt(abs(sum(∫(error' ⋅ error)*dΩ)/sum(∫((fex ∘ xp)' ⋅ (fex ∘ xp))*dΩ)))
end

# Compute the energy-norm-relative error in the entire computational domain
function energy_error(fh, fex, dfh, dfex, dΩ, xp, k)
    error = fh - fex
    error_grad = dfh - dfex
    return 100 * sqrt(abs(sum(∫(k^2*error ⋅ error + error_grad' ⋅ error_grad)*dΩ)/sum(∫(k^2*(fex ∘ xp) ⋅ (fex ∘ xp) + (dfex ∘ xp)' ⋅ (dfex ∘ xp))*dΩ)))
end

# Compute the energy-norm-relative error in the entire computational domain
function energy_error_vector(fh, fex, dfh, dfex, dΩ, xp, k)
    error = fh - fex
    error_grad = dfh - dfex
    return 100 * sqrt(abs(sum(∫(k^2*error' ⋅ error + error_grad ⋅ error_grad)*dΩ)/sum(∫(k^2*(fex ∘ xp)' ⋅ (fex ∘ xp) + (dfex ∘ xp) ⋅ (dfex ∘ xp))*dΩ)))
end

# Compute all the projections of the pressure field using the displacement field from the Raviart-Thomas discretization
function compute_from_fem(k, N, orderFE, orderFE_postprocessing, exact_solution, boundary_data)

    # Get the exact solution from the dictinoray
    pex = exact_solution["ExactPressure"]
    grad_pex = exact_solution["ExactGradientPressure"]
    uex = exact_solution["ExactDisplacement"]

    # Generate the mesh for a given mesh size
    model = generate_mesh(N)

    # Compute the solution of displacement-based formulation using first-order Raviart-Thomas discretization
    uh = RT_solve(k, orderFE, model, boundary_data)

    # Compute the solution of pressure-based formulation using Lagrange discretization
    ph = CG_solve(k, orderFE_postprocessing, model, boundary_data)

    # Compute the pressure projections
    ph_L2_DG = L2_DG_projection(uh, orderFE_postprocessing, model)
    ph_L2_CG = L2_CG_projection(uh, orderFE_postprocessing, model)
    ph_H1 = H1_projection(k, uh, orderFE_postprocessing, model)
    ph_H2 = H2_projection_DG(k, uh, orderFE_postprocessing, model, N)

    # Define the measure for computing the error with high-oscillatory solutions
    degree_ex = 16*orderFE_postprocessing # Quadrature degree

    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain
    dΩex = Measure(Ω, degree_ex)
    xp = get_physical_coordinate(Ω)

    # Compute the L^2-relative error in the entire computational domain
    error_L2_u = L2_error_vector(uh, uex, dΩex, xp)

    # Compute the energy-norm-relative error in the entire computational domain
    error_energy_u = energy_error_vector(uh, uex, -divergence(uh), pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain for ph_L2_DG
    error_L2_proj_DG = L2_error(ph_L2_DG, pex, dΩex, xp)

    # Compute the H^1-relative error in the entire computational domain for ph_L2_DG
    error_energy_proj_DG = energy_error(ph_L2_DG, pex, ∇(ph_L2_DG), grad_pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain for ph_L2_CG
    error_L2_proj_CG = L2_error(ph_L2_CG, pex, dΩex, xp)

    # Compute the H^1-relative error in the entire computational domain for ph_L2_CG
    error_energy_proj_CG = energy_error(ph_L2_CG, pex, ∇(ph_L2_CG), grad_pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain for ph_H1
    error_L2_proj_H1 = L2_error(ph_H1, pex, dΩex, xp)

    # Compute the H1-relative error in the entire computational domain for ph_H1
    error_energy_proj_H1 = energy_error(ph_H1, pex, ∇(ph_H1), grad_pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain for ph_H1
    error_L2_proj_H2_DG = L2_error(ph_H2, pex, dΩex, xp)

    # Compute the H1-relative error in the entire computational domain for ph_H1
    error_energy_proj_H2_DG = energy_error(ph_H2, pex, ∇(ph_H2), grad_pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain
    error_L2_p = L2_error(ph, pex, dΩex, xp)

    # Compute the energy-norm-relative error in the entire computational domain
    error_energy_p = energy_error(ph, pex, ∇(ph), grad_pex, dΩex, xp, k)

    # # Write the results to a VTK file with fileneame "results_N.vtu" with N the mesh size
    # writevtk(Ω,"./results/results_N=$(N)_order=$(orderFE).vtu", cellfields=["Re(uh)"=>real(uh), "Im(uh)"=>imag(uh),
    #                                         "Re(uex)"=>real(uex ∘ xp), "Im(uex)"=>imag(uex ∘ xp),
    #                                         "Re(ph)"=>real(ph), "Im(ph)"=>imag(ph),
    #                                         "Re(pex)"=>real(pex ∘ xp), "Im(pex)"=>imag(pex ∘ xp)])
    return error_L2_u, error_energy_u, error_L2_p, error_energy_p, error_L2_proj_DG, error_energy_proj_DG, error_L2_proj_CG, error_energy_proj_CG, error_L2_proj_H1, error_energy_proj_H1, error_L2_proj_H2_DG, error_energy_proj_H2_DG
end

# Compute all the projections of the pressure field using the displacement field from the interpolation of the exact solution
function compute_from_exact(k, N, orderFE_postprocessing, exact_solution, plot_flag=false)

    # Get the exact solution from the dictinoray
    pex = exact_solution["ExactPressure"]
    grad_pex = exact_solution["ExactGradientPressure"]
    uex = exact_solution["ExactDisplacement"]

    # Generate the mesh for a given mesh size
    model = generate_mesh(N)

    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain
    xp = get_physical_coordinate(Ω)

    # Compute the solution of displacement-based formulation using first-order Raviart-Thomas discretization
    reffe = ReferenceFE(raviart_thomas, Float64, orderFE_postprocessing-1)
    V_RT = TestFESpace(model, reffe, conformity=:Hdiv, vector_type=Vector{ComplexF64})
    uh = interpolate_everywhere(uex, V_RT)

    # Compute the solution of pressure-based formulation using Lagrange discretization
    reffe = ReferenceFE(lagrangian, Float64, orderFE_postprocessing)
    V_CG = TestFESpace(model, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
    ph = interpolate_everywhere(pex, V_CG)

    # Compute the pressure projections from exact solution
    ph_L2_DG = L2_DG_projection(uh, orderFE_postprocessing, model)
    ph_L2_CG = L2_CG_projection(uh, orderFE_postprocessing, model)
    ph_H1 = H1_projection(k, uh, orderFE_postprocessing, model)
    ph_H2 = H2_projection_DG(k, uh, orderFE_postprocessing, model, N)

    # Define the measure for computing the error with high-oscillatory solutions
    degree_ex = 16*orderFE_postprocessing # Quadrature degree

    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain
    dΩex = Measure(Ω, degree_ex)
    xp = get_physical_coordinate(Ω)

    # Compute the L^2-relative error in the entire computational domain
    error_L2_u = L2_error_vector(uh, uex, dΩex, xp)

    # Compute the energy-norm-relative error in the entire computational domain
    error_energy_u = energy_error_vector(uh, uex, -divergence(uh), pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain for ph_L2_DG
    error_L2_proj_DG = L2_error(ph_L2_DG, pex, dΩex, xp)

    # Compute the H^1-relative error in the entire computational domain for ph_L2_CG
    error_energy_proj_DG = energy_error(ph_L2_DG, pex, ∇(ph_L2_DG), grad_pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain for ph_L2_CG
    error_L2_proj_CG = L2_error(ph_L2_CG, pex, dΩex, xp)

    # Compute the H^1-relative error in the entire computational domain for ph_L2_CG
    error_energy_proj_CG = energy_error(ph_L2_CG, pex, ∇(ph_L2_CG), grad_pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain for ph_H1
    error_L2_proj_H1 = L2_error(ph_H1, pex, dΩex, xp)

    # Compute the H1-relative error in the entire computational domain for ph_H1
    error_energy_proj_H1 = energy_error(ph_H1, pex, ∇(ph_H1), grad_pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain for ph_H1
    error_L2_proj_H2_DG = L2_error(ph_H2, pex, dΩex, xp)

    # Compute the H1-relative error in the entire computational domain for ph_H1
    error_energy_proj_H2_DG = energy_error(ph_H2, pex, ∇(ph_H2), grad_pex, dΩex, xp, k)

    # Compute the L^2-relative error in the entire computational domain
    error_L2_p = L2_error(ph, pex, dΩex, xp)

    # Compute the energy-norm-relative error in the entire computational domain
    error_energy_p = energy_error(ph, pex, ∇(ph), grad_pex, dΩex, xp, k)

    # Plot using GridapMakie some fields only if plot_flag is true
    if plot_flag==true
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
        fig , ax , plt = plot(Ω, real(pex ∘ xp))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        # cbar = Colorbar(fig[2,1], plt, label=L"\mathrm{Re}(p_{\mathrm{ex}})", vertical=false)
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(p_{\mathrm{ex}})")
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_pex_N=$(N).pdf", fig)
    end

    # # Write the results to a VTK file with fileneame "results_N.vtu" with N the mesh size
    # writevtk(Ω,"./results/results_N=$(N)_order=$(orderFE).vtu", cellfields=["Re(uh)"=>real(uh), "Im(uh)"=>imag(uh),
    #                                         "Re(uex)"=>real(uex ∘ xp), "Im(uex)"=>imag(uex ∘ xp),
    #                                         "Re(ph)"=>real(ph), "Im(ph)"=>imag(ph),
    #                                         "Re(pex)"=>real(pex ∘ xp), "Im(pex)"=>imag(pex ∘ xp)])
    return error_L2_u, error_energy_u, error_L2_p, error_energy_p, error_L2_proj_DG, error_energy_proj_DG, error_L2_proj_CG, error_energy_proj_CG, error_L2_proj_H1, error_energy_proj_H1, error_L2_proj_H2_DG, error_energy_proj_H2_DG
end
