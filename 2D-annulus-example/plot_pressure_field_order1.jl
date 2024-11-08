"""
Compute the solution of displacement-based using first-order Raviart-Thomas discretization
with the Finite Element Method in Gridap.jl
"""

# Load the package to use the generate_mesh function
include("meshing.jl")

# Load the computing.jl file to compute the relative errors
include("../computing.jl")

# Define the range of N values: logspace between 1e1 and 500
N_value = 25
k = 50.0 # wavenumber
orderFE_postprocessing = 1 # order of the finite element for postprocessing

# Define the exact solution# Compute the analytical solution and the displacement error
p0 = 1.0 # Pressure amplitude
pex(x) = p0 * hankelh1(0, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k) # =-divergence(uex(x))
uex(x) = VectorValue(-x[1]/(sqrt(x[1]^2+x[2]^2)) / k * p0 * hankelh1(1, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k), 
                        -x[2]/(sqrt(x[1]^2+x[2]^2)) / k * p0 * hankelh1(1, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k))
grad_pex(x) = VectorValue(-x[1]/(sqrt(x[1]^2+x[2]^2)) * k * p0 * hankelh1(1, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k), 
                         -x[2]/(sqrt(x[1]^2+x[2]^2)) * k * p0 * hankelh1(1, k*sqrt(x[1]^2+x[2]^2)) / hankelh1(0, k)) # =k^2*uex(x)
exact_solution = Dict("ExactPressure"=>pex, "ExactDisplacement"=>uex, "ExactGradientPressure"=>grad_pex)

# Compute the relative errors and plot
error_L2_u, error_H1_u, error_L2_p, error_H1_p, error_L2_proj_DG, error_H1_proj_DG, error_L2_proj_CG, error_H1_proj_CG, error_L2_proj_H1, error_H1_proj_H1, error_L2_proj_H2_DG, error_H1_proj_H2_DG = compute_from_exact(k, N_value, orderFE_postprocessing, exact_solution, true)