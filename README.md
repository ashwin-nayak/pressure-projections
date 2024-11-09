Pressure Projections
====================

Code supplement for article "**Accuracy of post-processing projections for displacement based finite element simulations in room acoustics**".

This repository contains scripts to reproduce all the numerical results included as benchmark cases in the article.
The implementation is in the **Julia** programming language, primarily using [_Gridap_](https://github.com/gridap/Gridap.jl) package.

Details:
  - Weak formulation for pressure $L^2$ and $H^1$-projections.
  - Continuous and discontinuous Galerkin discretizations for computing the pressure field
  - First and second-order finite element discretizations are analyzed
  - Error convergence and computational performance is illustrated with different benchmarks


Authors:
  - **Ashwin S. Nayak**
  [![orcid](https://img.shields.io/badge/%20-orcid-black?logo=orcid&style=plastic)](https://orcid.org/0000-0002-9855-2377)
  [![github](https://img.shields.io/badge/%20-github-black?logo=github&style=plastic)](https://github.com/ashwin-nayak)
  - **Andrés Prieto**
  [![orcid](https://img.shields.io/badge/%20-orcid-black?logo=orcid&style=plastic)](https://orcid.org/0000-0002-4399-6878)
  [![github](https://img.shields.io/badge/%20-github-black?logo=github&style=plastic)](https://github.com/maprieto)
  - **Daniel Fernández-Comesaña**
  [![orcid](https://img.shields.io/badge/%20-orcid-black?logo=orcid&style=plastic)](https://orcid.org/0000-0003-3286-6637)
  [![github](https://img.shields.io/badge/%20-github-black?logo=github&style=plastic)](https://github.com/fernandez-microflown)

Affiliation:
  - [Centre for Mathematical Research and Technology of Galicia (CITMaGa)](https://citmaga.gal/en/home), A Coruña, Spain.

Repository:
  - https://github.com/ashwin-nayak/pressure-projections

License:
  - MIT, see [`LICENSE.txt`](LICENSE.txt).

<!-- Table of Contents -->
<details open>
<summary>Table of Contents</summary>
<ul>
  <li><a href="#getting-started">Getting Started</a></li>
  <li><a href="#run-using-docker-container">Run using docker container</a></li>
  <li><a href="#compile-and-run">Run benchmarks in Julia</a>
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

To run the scripts, you must first set the environment and dependent Julia packages.
Instructions provided here are provided for standard UNIX distributions, but maybe easily adopted (but not tested) in other operating systems.

You may either choose to use,
  - [Build and run in containers](#run-using-docker-container): Recommended to primarily reproduce results
  - [Compile and Run](#compile-and-run): Recommended for customizations and further development.

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

## Compile and Run

### Configure Dependencies

Setup the right versions and proper environment to run the scripts by following along the links,

1. [Install Julia (`v1.11`)](https://docs.julialang.org/en/v1.11/)
2. [Install Dependencies](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project)

### Run Scripts

If a previous installation of Julia is available, and the dependencies are installed, then each benchmark can be run independently in its respective folder.
