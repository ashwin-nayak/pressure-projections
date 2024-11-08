"""
Computation of the different post-processing pressure projections solution and measurement of their errors
from a displacement field, which is either obtained from a Raviart-Thomas discretization (implemented in Gridap) but in the presence of a Gaussian noise field on a triangular mesh
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry

# Load special function package
using SpecialFunctions

# Load Random package for generating random numbers and modify immutable fields
using Random

# Load the plotting packages
using GridapMakie, CairoMakie, FileIO, LaTeXStrings

# Load the projecting and computing functions
include("meshing.jl")
include("../projecting.jl")
include("../computing.jl")

function compute_noise(k, orderFE, PPW, exact_solution)

    # Generate meshes for a given wavenumber convergence_data_from_noise_k
    λ = 2.0 * π / k # Wavelength
    L = 1.0 # Domain length

    # Compute the solution of displacement-based formulation using first-order Raviart-Thomas discretization
    N_noise = Int(floor(6 * L / λ)) # Mesh size for the noise field
    model_noise = generate_mesh(N_noise)
    reffe = ReferenceFE(raviart_thomas, Float64, orderFE-1)
    V_RT_noise = TestFESpace(model_noise, reffe, conformity=:Hdiv, vector_type=Vector{ComplexF64})

    # Generate the noise field
    rng = MersenneTwister(1234) # seed for the random number generator
    values_noise = randn(rng, ComplexF64, V_RT_noise.nfree)
    noise = FEFunction(V_RT_noise, AbstractVector{ComplexF64}(values_noise))
    Ω = Triangulation(model_noise) # Computational domain
    xp = get_physical_coordinate(Ω)
    #noise_interp = Interpolable(noise)
    noise_interp = noise ∘ xp

    # Compute the solution of displacement-based formulation using first-order Raviart-Thomas discretization
    N = Int(floor(PPW * L / λ)) # Mesh size for the fluid domain
    model = generate_mesh(N)
    V_RT = TestFESpace(model, reffe, conformity=:Hdiv, vector_type=Vector{ComplexF64})
    noise_field = interpolate_everywhere(noise_interp, V_RT)

    # Get the exact solution as an interpolable function on the Raviart-Thomas discrete space
    uex = exact_solution["ExactDisplacement"]
    uh_interp = interpolate_everywhere(uex, V_RT)

    # Return the noise field, the interpolated displacement field, and the discrete space
    return N, model, V_RT, noise_field, uh_interp
end

function compute_pressure_from_noise(k, N, SNR, orderFE_postprocessing, model, V_RT, noise_field, uh_interp, plot_flag=false)

    # Get the exact solution from the dictinoray
    pex = exact_solution["ExactPressure"]
    grad_pex = exact_solution["ExactGradientPressure"]
    uex = exact_solution["ExactDisplacement"]

    # Compute displacement field with noise
    uh = FEFunction(V_RT, AbstractVector{ComplexF64}(uh_interp.free_values .* (1.0 .+ noise_field.free_values ./ SNR)))

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
    
    # # Write the results to a VTK file with fileneame "results_N.vtu" with N the mesh size
    # writevtk(Ω,"./results/results_N=$(N)_SNR=$(SNR).vtu", cellfields=["Re(uh)"=>real(uh), "Im(uh)"=>imag(uh),
    #         "Re(uex)"=>real(uex ∘ xp), "Im(uex)"=>imag(uex ∘ xp),
    #         "Re(u_interp)"=>real(uh_interp), "Im(u_interp)"=>imag(uh_interp),
    #         "Re(noise)"=>real(noise_field), "Im(noise)"=>imag(noise_field),
    #         "Re(ph_L2_DG)"=>real(ph_L2_DG), "Im(ph_L2_DG)"=>imag(ph_L2_DG),
    #         "Re(ph_L2_CG)"=>real(ph_L2_CG), "Im(ph_L2_CG)"=>imag(ph_L2_CG),
    #         "Re(ph_H1)"=>real(ph_H1), "Im(ph_H1)"=>imag(ph_H1),
    #         "Re(ph_H2)"=>real(ph_H2), "Im(ph_H2)"=>imag(ph_H2),
    #         "Re(pex)"=>real(pex ∘ xp), "Im(pex)"=>imag(pex ∘ xp)])
    
    # Plot using GridapMakie some fields only if plot_flag is true
    if plot_flag==true
        # Plot displacement field with noise
        fig , ax , plt = plot(Ω, real(uh_interp ⋅ VectorValue(0.0, 1.0)))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        xlims!(ax, 0.0, 1.0)
        ylims!(ax, 0.0, 1.0) 
        plt.colorrange=(-5e-2,5e-2)
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(\mathbf{u}_{\mathrm{ex}})\cdot\mathbf{e}_{y}", ticks=[-5e-2, -0.025, 0.0, 0.025, 5e-2])
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_uh_interp_N=$(N)_SNR=$(SNR).pdf", fig)

        # Plot displacement field with noise
        fig , ax , plt = plot(Ω, real(uh ⋅ VectorValue(0.0, 1.0)))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        xlims!(ax, 0.0, 1.0)
        ylims!(ax, 0.0, 1.0) 
        plt.colorrange=(-5e-2,5e-2)
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(\mathbf{u}_{\mathrm{SNR}})\cdot\mathbf{e}_{y}", ticks=[-5e-2, -0.025, 0.0, 0.025, 5e-2])
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_uh_N=$(N)_SNR=$(SNR).pdf", fig)

        # Plot pressure field (exact): pex
        fig , ax , plt = plot(Ω, real(pex ∘ xp))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        xlims!(ax, 0.0, 1.0)
        ylims!(ax, 0.0, 1.0)
        plt.colorrange=(-1.0,1.0)
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(p_{\mathrm{ex}})", ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
        cbar.limits = (-1.0, 1.0)
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_pex_N=$(N)_SNR=$(SNR).pdf", fig)
        
        # Plot pressure field with noise: ph_L2_DG
        fig , ax , plt = plot(Ω, real(ph_L2_DG))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        xlims!(ax, 0.0, 1.0)
        ylims!(ax, 0.0, 1.0)
        plt.colorrange=(-1.0,1.0)
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(\Pi_{h,\mathrm{DG}_{0}}^{\mathrm{L}^2})", ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
        cbar.limits = (-1.0, 1.0)
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_ph_L2_DG_N=$(N)_SNR=$(SNR).pdf", fig)

        # Plot pressure field with noise: ph_L2_CG
        fig , ax , plt = plot(Ω, real(ph_L2_CG))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        xlims!(ax, 0.0, 1.0)
        ylims!(ax, 0.0, 1.0)
        plt.colorrange=(-1.0,1.0)
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{L}^2})", ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
        cbar.limits = (-1.0, 1.0)
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_ph_L2_CG_N=$(N)_SNR=$(SNR).pdf", fig)

        # Plot pressure field with noise: ph_H1
        fig , ax , plt = plot(Ω, real(ph_H1))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        xlims!(ax, 0.0, 1.0)
        ylims!(ax, 0.0, 1.0)
        plt.colorrange=(-1.0,1.0)
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{H}^1})", ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
        cbar.limits = (-1.0, 1.0)
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_ph_H1_N=$(N)_SNR=$(SNR).pdf", fig)

        # Plot pressure field with noise: ph_H2
        fig , ax , plt = plot(Ω, real(ph_H2))
        # Set x- and y-axis labels
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        xlims!(ax, 0.0, 1.0)
        ylims!(ax, 0.0, 1.0)
        plt.colorrange=(-1.0,1.0)
        # Axis equal
        ax.aspect = DataAspect()
        # Colorbar with limits and fixed ticks
        cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(\Pi_{h,\mathrm{DG}_{1}}^{\mathrm{H}^1})", ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
        cbar.limits = (-1.0, 1.0)
        # Save the figure to pdf format with 300 dpi
        save("./results/plot_ph_H2_N=$(N)_SNR=$(SNR).pdf", fig)

    end

    return error_L2_u, error_energy_u, error_L2_proj_DG, error_energy_proj_DG, error_L2_proj_CG, error_energy_proj_CG, error_L2_proj_H1, error_energy_proj_H1, error_L2_proj_H2_DG, error_energy_proj_H2_DG
end
