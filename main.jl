# Copyright (c) 2024 Ashwin. S. Nayak, Andrés Prieto, Daniel Fernández Comesaña
#
# This file is part of pressure-projections
#
# SPDX-License-Identifier:  MIT

# List the subdirectories with main.jl
subdirs = [
    "1D-example",
    "2D-annulus-example",
    "2D-square-structured-example",
    "2D-square-unstructured-example",
    "2D-square-noise-example",
    "2D-room-example"
]

# Loop through each subdirectory and run main.jl
for subdir in subdirs
    sub_main_path = joinpath(subdir, "main.jl")
    # Check if main.jl exists in the specified subdirectory
    if isfile(sub_main_path)
        println("Running main.jl in subdirectory: $subdir")
        # Change to the subdirectory and execute main.jl in the subdirectory context
        cd(subdir) do
            Base.include(Main, sub_main_path)
        end
    else
        println("No main.jl found in subdirectory: $subdir, skipping.")
    end
end
