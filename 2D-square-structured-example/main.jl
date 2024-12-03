# Copyright (c) 2024 Ashwin. S. Nayak, Andrés Prieto, Daniel Fernández Comesaña
#
# This file is part of pressure-projections
#
# SPDX-License-Identifier:  MIT

mkpath("meshes")
mkpath("results")
include("convergence_plot_from_fem_order1.jl")
include("plot_fem_order1.jl")
include("convergence_plot_from_fem_order2.jl")
include("plot_fem_order2.jl")
include("convergence_plot_from_exact_order1.jl")
include("plot_exact_order1.jl")
include("convergence_plot_from_exact_order2.jl")
include("plot_exact_order2.jl")
include("plot_pressure_field_order1.jl")
