"""
Benchmarking and plot of the performance curves (CPU time and memory) obtained from the post-processing projection techniques
when the displacement field is computed from the second-order Raviart-Thomas FE discretizations
"""
# Load the packages
using BenchmarkTools

# Load the package to use the generate_mesh function
include("meshing.jl")

# Load the computing.jl file to compute the CPU times
include("computing_with_CPUtime.jl")

# Define the range of N values: logspace between 1e1 and 1e4
N_values = Integer.(round.(10.0.^range(log10(10), stop=log10(300), length=15)))
k=50.0 # wave number (2, 20, 100)
orderFE = 2 # order of the finite element
orderFE_postprocessing = 2 # order of the finite element for postprocessing

# Dictionaonary  with the boundary data associated to each boundary tag
p0 = 1.0 # Pressure amplitude
Zs = 1.0 # Specific surface impedance
uD(x) = VectorValue(0.0, 0.0)
pD(x) = p0
boundary_data = Dict("RigidBoundary"=>uD, "WallPressure"=>pD, "AbsorbingBoundary"=>Zs)

# Define the arrays to store the CPU times
CPUtime_u = zeros(length(N_values), 2)
Memory_u = zeros(length(N_values), 2)
CPUtime_p = zeros(length(N_values), 2)
Memory_p = zeros(length(N_values), 2)
CPUtime_proj_DG = zeros(length(N_values), 2)
Memory_proj_DG = zeros(length(N_values), 2)
CPUtime_proj_CG = zeros(length(N_values), 2)
Memory_proj_CG = zeros(length(N_values), 2)
CPUtime_proj_H1 = zeros(length(N_values), 2)
Memory_proj_H1 = zeros(length(N_values), 2)
CPUtime_proj_H2_DG = zeros(length(N_values), 2)
Memory_proj_H2_DG = zeros(length(N_values), 2)

# Perform sweeping in h_values with a fixed k value
for j in 1:size(N_values, 1)
    # Compute the CPU times
    benchmark_u, benchmark_p, benchmark_proj_DG, benchmark_proj_CG, benchmark_proj_H1, benchmark_proj_H2_DG  = compute_from_fem_with_CPUtime(k, N_values[j], orderFE, orderFE_postprocessing, boundary_data, false)
    # Compute the average and BenchmarkTools.std CPU times
    CPUtime_u[j, 1] = BenchmarkTools.mean(benchmark_u).time
    CPUtime_u[j, 2] = BenchmarkTools.std(benchmark_u).time
    CPUtime_p[j, 1] = BenchmarkTools.mean(benchmark_p).time
    CPUtime_p[j, 2] = BenchmarkTools.std(benchmark_p).time
    CPUtime_proj_DG[j, 1] = BenchmarkTools.mean(benchmark_proj_DG).time
    CPUtime_proj_DG[j, 2] = BenchmarkTools.std(benchmark_proj_DG).time
    CPUtime_proj_CG[j, 1] = BenchmarkTools.mean(benchmark_proj_CG).time
    CPUtime_proj_CG[j, 2] = BenchmarkTools.std(benchmark_proj_CG).time
    CPUtime_proj_H1[j, 1] = BenchmarkTools.mean(benchmark_proj_H1).time
    CPUtime_proj_H1[j, 2] = BenchmarkTools.std(benchmark_proj_H1).time
    CPUtime_proj_H2_DG[j, 1] = BenchmarkTools.mean(benchmark_proj_H2_DG).time
    CPUtime_proj_H2_DG[j, 2] = BenchmarkTools.std(benchmark_proj_H2_DG).time
    # Compute the average and BenchmarkTools.std of allocated Memory
    Memory_u[j, 1] = benchmark_u.memory
    Memory_u[j, 2] = benchmark_u.allocs
    Memory_p[j, 1] = benchmark_p.memory
    Memory_p[j, 2] = benchmark_p.allocs
    Memory_proj_DG[j, 1] = benchmark_proj_DG.memory
    Memory_proj_DG[j, 2] = benchmark_proj_DG.allocs
    Memory_proj_CG[j, 1] = benchmark_proj_CG.memory
    Memory_proj_CG[j, 2] = benchmark_proj_CG.allocs
    Memory_proj_H1[j, 1] = benchmark_proj_H1.memory
    Memory_proj_H1[j, 2] = benchmark_proj_H1.allocs
    Memory_proj_H2_DG[j, 1] = benchmark_proj_H2_DG.memory
    Memory_proj_H2_DG[j, 2] = benchmark_proj_H2_DG.allocs
end
h_values = 1.0./N_values

# Save errors to a JDL2 file
using JLD2
save("CPUtime_data_from_fem_k=$(Integer(k))_orderFE=$(orderFE)_orderFEpostprocessing=$(orderFE_postprocessing).jld2", 
    "h_values", h_values, "k", k, "orderFE", orderFE, "orderFE_postprocessing", orderFE_postprocessing, 
    "CPUtime_u", CPUtime_u, "Memory_u", Memory_u, "CPUtime_p", CPUtime_p, "Memory_p", Memory_p,
    "CPUtime_proj_DG", CPUtime_proj_DG, "Memory_proj_DG", Memory_proj_DG, "CPUtime_proj_CG", CPUtime_proj_CG, "Memory_proj_CG", Memory_proj_CG,
    "CPUtime_proj_H1", CPUtime_proj_H1, "Memory_proj_H1", Memory_proj_H1, "CPUtime_proj_H2_DG", CPUtime_proj_H2_DG, "Memory_proj_H2_DG", Memory_proj_H2_DG)

# Plot the CPU times with latex labels and log scale in Makie
# include("plot_fem_order1.jl")