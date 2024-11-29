# Read JLD2 data
using JLD2
# data = load("results/CPUtime_data_from_fem_k=50_orderFE=1_orderFEpostprocessing=1.jld2")
data = load("results/CPUtime_data_from_fem_k=50_orderFE=2_orderFEpostprocessing=2.jld2")
N_values = 1.0./data["h_values"]

# Save errors to a CSV file
using DelimitedFiles

# CPUTime Plot

csv_header = [
    "2pibyk_N_values",
    "CPUtime_proj_DG_Mean",
    "CPUtime_proj_DG_StdMax",
    "CPUtime_proj_DG_StdMin",
    "CPUtime_proj_CG_Mean",
    "CPUtime_proj_CG_StdMax",
    "CPUtime_proj_CG_StdMin",
    "CPUtime_proj_H1_Mean",
    "CPUtime_proj_H1_StdMax",
    "CPUtime_proj_H1_StdMin",
    "CPUtime_proj_H2_DG_Mean",
    "CPUtime_proj_H2_DG_StdMax",
    "CPUtime_proj_H2_DG_StdMin"
]

# CPUtimes from ns to s
CPUtime_u = data["CPUtime_u"]/1e9
CPUtime_p = data["CPUtime_p"]/1e9
CPUtime_proj_DG = data["CPUtime_proj_DG"]/1e9
CPUtime_proj_CG = data["CPUtime_proj_CG"]/1e9
CPUtime_proj_H1 = data["CPUtime_proj_H1"]/1e9
CPUtime_proj_H2_DG = data["CPUtime_proj_H2_DG"]/1e9

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

csv_data = hcat(
    2*pi/data["k"].*N_values,
    CPUtime_proj_DG[:,1],
    CPUtime_proj_DG[:,1]+CPUtime_proj_DG[:,2],
    CPUtime_proj_DG[:,1]-CPUtime_proj_DG[:,2],
    CPUtime_proj_CG[:,1],
    CPUtime_proj_CG[:,1]+CPUtime_proj_CG[:,2],
    CPUtime_proj_CG[:,1]-CPUtime_proj_CG[:,2],
    CPUtime_proj_H1[:,1],
    CPUtime_proj_H1[:,1]+CPUtime_proj_H1[:,2],
    CPUtime_proj_H1[:,1]-CPUtime_proj_H1[:,2],
    CPUtime_proj_H2_DG[:,1],
    CPUtime_proj_H2_DG[:,1]+CPUtime_proj_H2_DG[:,2],
    CPUtime_proj_H2_DG[:,1]-CPUtime_proj_H2_DG[:,2]
)

open("results/cputime_k=$(Integer(data["k"]))_orderFEpostprocessing=$(data["orderFE_postprocessing"]).csv", "w") do io
    writedlm(io, [csv_header])
    writedlm(io, csv_data)
end


# Memory Plot
csv_header = [
    "2pibyk_N_values",
    "Memory_proj_DG", "Memory_proj_CG",
    "Memory_proj_H1", "Memory_proj_H2_DG"
]

csv_data = hcat(
    2*pi/data["k"].*N_values,
    data["Memory_proj_DG"][:,1]/1e6, data["Memory_proj_CG"][:,1]/1e6,
    data["Memory_proj_H1"][:,1]/1e6, data["Memory_proj_H2_DG"][:,1]/1e6
)

open("results/memory_k=$(Integer(data["k"]))_orderFEpostprocessing=$(data["orderFE_postprocessing"]).csv", "w") do io
    writedlm(io, [csv_header])
    writedlm(io, csv_data)
end
