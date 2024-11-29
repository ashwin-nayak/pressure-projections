# Read JLD2 data
using JLD2
# data = load("results/convergence_data_from_exact_k=50_orderFEpostprocessing=1.jld2")
# data = load("results/convergence_data_from_exact_k=50_orderFEpostprocessing=2.jld2")
data = load("results/convergence_data_from_fem_k=50_orderFE=1_orderFEpostprocessing=1.jld2")
# data = load("results/convergence_data_from_fem_k=50_orderFE=2_orderFEpostprocessing=2.jld2")

# Save errors to a CSV file
using DelimitedFiles
csv_header = [
    "2pibyk_N_values", "h_values",
    "error_L2_u", "error_H1_u",
    "error_L2_p", "error_H1_p",
    "error_L2_proj_DG", "error_H1_proj_DG",
    "error_L2_proj_CG", "error_H1_proj_CG",
    "error_L2_proj_H1", "error_H1_proj_H1",
    "error_L2_proj_H2_DG", "error_H1_proj_H2_DG"
]
N_values = 1.0./data["h_values"]

csv_data = hcat(
    2*pi/data["k"].*N_values, data["h_values"],
    data["error_L2_u"], data["error_H1_u"],
    data["error_L2_p"], data["error_H1_p"],
    data["error_L2_proj_DG"], data["error_H1_proj_DG"],
    data["error_L2_proj_CG"], data["error_H1_proj_CG"],
    data["error_L2_proj_H1"], data["error_H1_proj_H1"],
    data["error_L2_proj_H2_DG"], data["error_H1_proj_H2_DG"]
)

# open("results/convergence_data_from_exact_k=$(Integer(data["k"]))_orderFEpostprocessing=$(Integer(data["orderFE_postprocessing"])).csv", "w") do io
#     writedlm(io, [csv_header])
#     writedlm(io, csv_data)
# end

open("results/convergence_data_from_fem_k=$(Integer(data["k"]))_orderFE=$(data["orderFE"])_orderFEpostprocessing=$(data["orderFE_postprocessing"]).csv", "w") do io
    writedlm(io, [csv_header])
    writedlm(io, csv_data)
end
