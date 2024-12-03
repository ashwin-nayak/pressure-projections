1D Example
==========

One-dimensional planewave propagation in an elongated rectangle.


- Domain, $\Omega = {(0,L)\times(0,h)}$
  - Length, $L=1$
  - Width, $h$ is the mesh-width.
- Boundary conditions
  - Rigid boundary, $\boldsymbol{u} \cdot \boldsymbol{n} = 0$ at $y=0$ and $y=h$,
  - Pressure $p = 1.0$ at $x=0$,
  - Impedance $Z = 1.0$ at $x=L$.

The boundary conditions are expressly chosen such that the exact solution can be computed analytically and the errors can be computed.

The exact solution here is given as,

$$
\begin{aligned}
    \boldsymbol{u}\_{\mathrm{ex}} (x,y) &= {\frac{{p}\_{0}}{k^2}} \exp (\mathrm{i} k x) \boldsymbol{e}\_{x}, \\
    {p}\_{\mathrm{ex}} (x,y)  &= {p}\_{0} \exp(\mathrm{i} k x).
\end{aligned}
$$

## Run

The following instructions run the example **assuming you are in the current directory**.

An automated script [`main.jl`](main.jl) is provided to run scripts in a particular order,

```bash
julia +1.11.0 --project=.. main.jl
```
This should generate the meshes of various resolutions in the `meshes` folder. All output fields are written within `.jld2` files and visualizations in `.pdf` files within the `results` folder.
