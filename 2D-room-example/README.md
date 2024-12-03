# Two-dimensional Auditorium Case : Computational Performance evaluation

Acoustic propagation within a 2D room scenario.

- Custom room domain, $\Omega$ (see `meshing.jl` or the article for details)
- Boundary Conditions,
  - Rigid Boundaries on top and left,
  - Bottom Staircase boundary with impedance $Z=1.0$,
  - Pressure $p=1.0$ on a segment of right boundary $y=(1.5,2)$.

> [!IMPORTANT]
> No exact solution is known apriori for this scenario and hence, no analytical errors can be computed.

The computational time required to compute the pressure field using the first and second-order approximations are computed.

## Run

The following instructions run the example **assuming you are in the current directory**.

An automated script [`main.jl`](main.jl) is provided to run scripts in a particular order,

```bash
julia +1.11.0 --project=.. main.jl
```
This should generate the meshes of various resolutions in the `meshes` folder. All output fields are written within `.jld2` files and visualizations in `.pdf` files within the `results` folder.
