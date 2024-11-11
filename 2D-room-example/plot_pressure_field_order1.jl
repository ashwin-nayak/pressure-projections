"""
Plot the displacement field and the post-processed (H^1-projection) pressure field for an specific triangular mesh (N), a fixed FE order, and a given wavenumber k
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry

# Load special function package (only required for the 2D simulations)
using SpecialFunctions

# Load the plotting packages
using GridapMakie, CairoMakie, FileIO, LaTeXStrings

# Load the package to use the generate_mesh function
include("meshing.jl")

# Load the meshing.jl file to use the generate_mesh function
include("../projecting.jl")

# Define the range of N values: logspace between 1e1 and 1e4
N = 400 # mesh value of refinement
k=50.0 # wave number (2, 20, 100)
orderFE = 1 # order of the finite element
orderFE_postprocessing = 1 # order of the finite element for postprocessing

# Dictionaonary  with the boundary data associated to each boundary tag
p0 = 1.0 # Pressure amplitude
Zs = 1.0 # Specific surface impedance
uD(x) = VectorValue(0.0, 0.0)
pD(x) =p0 * (x[2]>1.5) * (x[2]<2.0)
boundary_data = Dict("RigidBoundary"=>uD, "WallPressure"=>pD, "AbsorbingBoundary"=>Zs)

# Generate the mesh for a given mesh size
model = generate_mesh(N)

# Compute the solution of displacement-based formulation using first-order Raviart-Thomas discretization
uh = RT_solve(k, orderFE, model, boundary_data)

# Compute the pressure projections
ph_H1 = H1_projection(k, uh, orderFE_postprocessing, model)

# Plot mesh
Ω = Triangulation(model) # Computational domain
fig, ax, plt = wireframe(Ω, color=:black, linewidth=1.0)
# Set x- and y-axis labels
ax.xlabel = L"x"
ax.ylabel = L"y"
# Axis equal
ax.aspect = DataAspect()
# Save the figure to pdf format with 300 dpi
save("./results/mesh_N=$(N).pdf", fig)
# Plot pressure field

fig , ax , plt = plot(Ω, real(ph_H1))
# Set x- and y-axis labels
ax.xlabel = L"x"
ax.ylabel = L"y"
# Axis equal
ax.aspect = DataAspect()
# Colorbar with limits and fixed ticks
cbar = Colorbar(fig[1,2], plt, label=L"\mathrm{Re}(\Pi_{h,\mathrm{CG}_{1}}^{\mathrm{H}^1})", tellheight=true)
rowsize!(fig.layout, 1, ax.scene.px_area[].widths[2])
# Save the figure to pdf format with 300 dpi
save("./results/plot_p_H1_N=$(N).pdf", fig)
