PressureProjections
=============

Code supplement for "Accuracy of post-processing projections for displacement based finite element simulations in room acoustics" 

A Julia implementation (based on Gridap package) of different $L^2$ and $H^1$-projections for computing the acoustic pressure from the particle displacement field. This repository contains the prototype codes to reproduce all the numerical results included as benchmark cases in the manuscript "Accuracy of post-processing projections for displacement based finite element simulations in room acoustics". 

Details:
  - Weak formulation for pressure $L^2$ and $H^1$-projections.
  - Continuous and discontinuous Galerkin discretizations for computing the pressure field
  - First and second-order finite element discretizations are analyzed
  - Error convergence and computational performance is illustrated with different benchmarks


Authors:
  - [Ashwin S. Nayak](https://orcid.org/0000-0002-9855-2377)
  - [Andrés Prieto](https://orcid.org/0000-0002-4399-6878)
  - [Daniel Fernández-Comesaña](https://orcid.org/0000-0003-3286-6637)

Affiliation:
  - [Centro of Mathematical Research and Technology of Galicia](https://citmaga.gal/en/home), A Coruña, Spain.

Repository:
  - https://github.com/ashwin-nayak/PressureProjections

License:
  - MIT, see [`LICENSE.txt`](LICENSE.txt).

<!-- Table of Contents -->
<details open>
<summary>Table of Contents</summary>
<ul>
  <li><a href="#getting-started">Getting Started</a></li>
  <li><a href="#run-using-docker-container">Run using docker container</a></li>
  <li><a href="#run-in-julia">Run benchmarks in Julia</a>
    <ul>
      <li><a href="#2D-convergence-examples">2D convergence examples</a></li>
      <li><a href="#2D-square-noise-example">2D square noise example</a></li>
      <li><a href="#2D-room-example">2D room example</a></li>
    </ul>
  </li>
</ul>
</details>

## Getting Started

The library is structured into the following folders, which correspond to each benchmark scenarios:

| Folder              | Information                                                                     |
|---------------------|---------------------------------------------------------------------------------|
| `1D-example`                      | One-dimensional planewave propagation in an elongated rectangle |
| `2D-annulus-example`              | Two-dimensional radial wave propagation in a quarter of an annular domain |
| `2D-square-structured-example`    | Two-dimensional planewave propagation in the unit square using a structured mesh |
| `2D-square-unstructured-example`  | Two-dimensional planewave propagation in the unit square using an unstructured mesh |
| `2D-square-noise-example`         | Two-dimensional planewave propagation in the unit square in the presence of noise |
| `2D-room-example`                 | Two-dimensional transmitted field in the interior of the cross-section of an auditorium |

To use the code, you must first set the environment and dependent Julia packages.
Instructions provided here are provided for standard UNIX distributions, but maybe easily adopted (but not tested) in other operating systems.

You may either choose to use,
  - [Build and run in containers](#run-using-docker-container): Recommended to primarily reproduce results
  - [Run in Julia](#run-in-julia): Recommended for a quick start browsing on the numerical results.

## Run using docker container

The project offers a [`Dockerfile`](Dockerfile) to automate the configuration and execution by,
setting relevant dependencies, hosting the source code, building executables and running the demos.

Necessary tools:
  - [Docker Desktop](https://docs.docker.com/desktop/) specifically Docker Engine (tested with `26.1.3`).

Build the `PressureProjections` image which contains the environment, copies the source code and builds executables

```bash
docker build -f Dockerfile --tag PressureProjections .
```

Run the container by sharing the main folder within and user id information to manage the generated file permissions,
```bash
docker run --rm -u $(id -u):$(id -g) -v ${PWD}/main:julia
```
This should run an instance of the Julia REPL where each benchmark can be executed. Each benchmark folder has an intentionally empty meshes and results folder where the meshes and the plots are storage.

> Ensure that the directory to store the `main` folder for this repository exists before running, since this needs to be mounted as shared volume within the host and container. This also preserves user permissions on the files that are created within the container.

## Run in Julia

If a previous installation of Julia is available, and the requirements of the [`Manifest.toml`](Manifest.toml) and [`Project.toml`](Project.toml) are satisfied, then each benchmark can be run independently in its respective folder (where they are available) as follows.

#### 2D convergence examples
On the benchmark cases without noise, the computation of the error curves and the projected pressure fields are computed running:

* Computation of the error curves from the first-order finite element discretizations using the finite element approximation of the displacement field
  ```bash
  julia convergence_plot_from_fem_order1.jl
  ```
* Computation of the error curves from the second-order finite element discretizations using the finite element approximation of the displacement field
  ```bash
  julia convergence_plot_from_fem_order2.jl
  ```
* Computation of the error curves from the first-order finite element discretizations using the interpolated values from the exact solution
  ```bash
  julia convergence_plot_from_exact_order1.jl
  ```
* Computation of the error curves from the second-order finite element discretizations using the interpolated values from the exact solution
  ```bash
  julia convergence_plot_from_exact_order2.jl
  ```
* Plot the interpolated values (first-order approximation) from the exact solution in a given simplicial mesh
  ```bash
  julia plot_pressure_field_order1.jl
  ```

Additionally, the two-dimensional benchmarks have specific purpose Julia scripts to plot the error curves, where the data is read from pre-computed values (already available in JLD2 files included in this repository):

* Plot the error curves from the first-order finite element discretizations using the finite element approximation of the displacement field
  ```bash
  julia plot_from_fem_order1.jl
  ```
* Plot the error curves from the second-order finite element discretizations using the finite element approximation of the displacement field
  ```bash
  julia plot_from_fem_order2.jl
  ```
* Plot the error curves from the first-order finite element discretizations using the interpolated values from the exact solution
  ```bash
  julia plot_from_exact_order1.jl
  ```
* Plot the error curves from the second-order finite element discretizations using the interpolated values from the exact solution
  ```bash
  julia plot_from_exact_order2.jl
  ```

#### 2D square noise example
The benchmark case with noisy values have some specific scripts to compute the performance of the pressure projections approaches:

* Computation of the error curves from the first-order finite element discretizations using the finite element approximation of the displacement field
  ```bash
  julia convergence_plot_from_noise_order1.jl
  ```
* Computation of the error curves from the second-order finite element discretizations using the finite element approximation of the displacement field
  ```bash
  julia convergence_plot_from_noise_order2.jl
  ```
* Plot the (first-order approximation) from pressure projections in the presence of noisy values and using a given simplicial mesh
  ```bash
  julia plot_solutions_with_noise_order1.jl
  ```
* Plot the error curves from the first-order finite element discretizations using noisy values
  ```bash
  julia plot_noise_order1.jl
  ```
* Plot the error curves from the second-order finite element discretizations using noisy values
  ```bash
  julia plot_noise_order2.jl
  ```

#### 2D room example
Finally, the auditorium benchmark case have similar scripts as those described above, but instead of plotting error curves, it measures and plots the CPU time and the allocation memory required for each projection.

