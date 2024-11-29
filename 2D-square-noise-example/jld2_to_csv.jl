# Read JLD2 data
using JLD2
# data = load("results/convergence_data_from_noise_k=25_orderFEpostprocessing=1.jld2")
data = load("results/convergence_data_from_noise_k=25_orderFEpostprocessing=2.jld2")


# Save errors to a CSV file
using DelimitedFiles
csv_header = [
    "SNR_values",
    "error_L2_u", "error_H1_u",
    "error_L2_proj_DG", "error_H1_proj_DG",
    "error_L2_proj_CG", "error_H1_proj_CG",
    "error_L2_proj_H1", "error_H1_proj_H1",
    "error_L2_proj_H2_DG", "error_H1_proj_H2_DG"
]

csv_data = hcat(
    data["SNR_values"],
    data["error_L2_u"], data["error_H1_u"],
    data["error_L2_proj_DG"], data["error_H1_proj_DG"],
    data["error_L2_proj_CG"], data["error_H1_proj_CG"],
    data["error_L2_proj_H1"], data["error_H1_proj_H1"],
    data["error_L2_proj_H2_DG"], data["error_H1_proj_H2_DG"]
)

open("results/convergence_data_from_noise_k=$(Integer(data["k"]))_orderFEpostprocessing=$(data["orderFE_postprocessing"]).csv", "w") do io
    writedlm(io, [csv_header])
    writedlm(io, csv_data)
end
