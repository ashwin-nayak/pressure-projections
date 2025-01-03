Pressure Projections
====================

Code supplement for article "**Accuracy of post-processing projections for displacement based finite element simulations in room acoustics**".

This repository contains scripts to reproduce all the numerical results included as benchmark cases in the article.
The implementation is in the **Julia** programming language, primarily using [_Gridap_](https://github.com/gridap/Gridap.jl) package.

Details:
  - Weak formulation for pressure $L^2$ and $H^1$-projections
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
  <li><a href="#run-in-julia">Run benchmarks in Julia</a></li>
  <li><a href="#development-container-for-developers">Development container (for developers)</a></li>
</ul>
</details>

## Getting Started

First, clone the repository or download the contents locally on your computer, and ensure you are working within project directory,
```bash
git clone git@github.com:ashwin-nayak/pressure-projections.git
cd pressure-projections
```

The library is structured into the following folders, which correspond to each benchmark scenarios:

| Folder              | Information                                                                     |
|---------------------|---------------------------------------------------------------------------------|
| [`1D-example`](1D-example)                      | One-dimensional planewave propagation in an elongated rectangle |
| [`2D-annulus-example`](2D-annulus-example)              | Two-dimensional radial wave propagation in a quarter of an annular domain |
| [`2D-square-structured-example`](2D-square-structured-example)    | Two-dimensional planewave propagation in the unit square using a structured mesh |
| [`2D-square-unstructured-example`](2D-square-unstructured-example)  | Two-dimensional planewave propagation in the unit square using an unstructured mesh |
| [`2D-square-noise-example`](2D-square-noise-example)         | Two-dimensional planewave propagation in the unit square in the presence of noise |
| [`2D-room-example`](2D-room-example)                 | Two-dimensional transmitted field in the interior of the cross-section of an auditorium |

To run the scripts, you must first set the environment and dependent Julia packages.
Instructions provided here are provided for standard UNIX distributions, but maybe easily adopted (but not tested) in other operating systems.

You may either choose to use,
  - [Build and run in containers](#run-using-docker-container): Recommended to primarily reproduce results
  - [Compile and Run](#run-in-julia): Recommended for customizations and further development.

> [!WARNING]
> The code execution requires a large amount of free RAM memory (>40 GB) and a runtime of about ~20 h (due to large problem sizes for highly resolved meshes).

## Run using docker container

The project offers a [`Dockerfile`](Dockerfile) to automate the configuration and execution by,
setting relevant dependencies, hosting the source code, building executables and running the demos.

Necessary tools:
  - [Docker Desktop](https://docs.docker.com/desktop/) specifically Docker Engine (tested with `27.3.1`).

Build the `pressure-projections` image which contains the environment, copies the source code and builds executables

```bash
docker build -f Dockerfile --tag pressure-projections .
```

Run a disposable container by sharing the project folder within,
```bash
docker run --rm \
  -u $(id -u):$(id -g) \
  -v $PWD:/src \
  -w /src \
  pressure-projections \
    julia \
      --project=/pressure-projections \
      --startup-file=no \
      --history-file=no \
      main.jl
```

This should run an instance of the Julia REPL where each benchmark can be executed.

> [!TIP]
> The flags ensure correct user permissions are set on the container-generated output files. Additionally, the container information is cleared at the end of the run.

## Run in Julia

Setup the right versions and proper environment to run the scripts by following along the links,

1. [Install Julia (`v1.11.0`)](https://docs.julialang.org/en/v1.11.0/) (not described here)
2. [Install Dependencies](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project)


> [!TIP]
> If [`juliaup`](https://github.com/JuliaLang/juliaup) is available (recommended), the right version can be installed simply by,
>```bash
>juliaup add 1.11.0
>```

> [!TIP]
> If you already have installed Julia before, verify the version you have by,
> ```shell
>julia -e 'using InteractiveUtils; versioninfo()'
>```

Install relevant dependencies to ensure proper environment,
```bash
julia +1.11.0 --project=. -e "using Pkg; Pkg.instantiate();"
```

Each benchmark can be run independently in its respective folder.
Alternatively, an automated script [`main.jl`](main.jl) is provided to run **all the benchmarks**,
```bash
julia +1.11.0 --project=. main.jl
```

## Development Container (for developers)

If you intend to develop the source code without modifying/installing any dependencies on your host computer, you can make use of [development containers](https://containers.dev/) for setting up the requisite environment.

A specification file [`.devcontainer.json`](.devcontainer.json) is provided for building the devcontainer and can be utilized by supporting editors.

Necessary tools :
- [Docker Engine](https://docs.docker.com/engine/install/) (tested with `27.3.1`),
- [Visual Studio Code](https://code.visualstudio.com/) with [Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension (support since 2019. Tested with version `1.94.2`).

Start VS Code, run the "Dev Containers: Open Folder in Container..." command from the Command Palette (F1), and select the project folder.
While done for the first time, this should start building the container and can take a while (~5 min) and subsequently open the project within the container.
The build is cached for subsequent runs and should be fairly quick thereon.

Proceed to [running the benchmarks](#run-in-julia).
