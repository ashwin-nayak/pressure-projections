"""
Plot the interpolated values of the noisy pressure field for an specific triangular mesh (N)
a fixed FE order, and a given wavenumber k
"""
# Load the computing.jl file to use the generate_mesh function
include("computing_with_noise.jl")

# Define the range of N values: logspace between 1e1 and 1e4
SNR_value = 2.0
k=25.0 # wave number (2, 20, 100)
orderFE = 1 # order of the finite element used to interpolate the exact solution
orderFE_postprocessing = 1 # order of the finite element for postprocessing
PPW = 20 # Points per Wavelength in th mesh

# Define the exact solution# Compute the analytical solution and the displacement error
p0 = 1.0 # Pressure amplitude
pex(x) = p0 * exp(1im*k*x[2]) # =-divergence(uex(x))
uex(x) = VectorValue(0.0, 1im * p0 / k * exp(1im*k*x[2]))
grad_pex(x) = VectorValue(0.0, 1im * p0 * k * exp(1im*k*x[2])) # =k^2*uex(x)
exact_solution = Dict("ExactPressure"=>pex, "ExactDisplacement"=>uex, "ExactGradientPressure"=>grad_pex)

# Compute noise field (normalized at [0,1])
N, model, V_RT, noise_field, uh_interp = compute_noise(k, orderFE, PPW, exact_solution)

# Perform sweeping in SNR_values with a fixed k value
error_L2_u, error_H1_u, error_L2_proj_DG, error_H1_proj_DG, error_L2_proj_CG, error_H1_proj_CG, error_L2_proj_H1, error_H1_proj_H1, error_L2_proj_H2_DG, error_H1_proj_H2_DG = compute_pressure_from_noise(k, N, SNR_value, orderFE_postprocessing, model, V_RT, noise_field, uh_interp, true)