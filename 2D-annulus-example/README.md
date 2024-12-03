# Two-dimensional radial problem

Radial wave propagation in an annular domain.

- Domain, $\Omega= \\\{(x,y)\in\mathbb{R}^2: x>0,\ y>0,\ 1<x^2+y^2<4\\\}$
  - i.e., between radius, $1< R < 2$.
- Boundary Conditions
  - Rigid boundary, $\boldsymbol{u} \cdot \boldsymbol{n} = 0$ at $x=0$ and $y=0$,
  - Pressure $p=1.0$ at $R=1$,
  - Impedance $Z=-\mathrm{i} \mathrm{H}^{(1)}\_{0}(2k)/\mathrm{H}^{(1)}\_{1}(2k)$ at $R=2$.

The boundary conditions are expressly chosen such that the exact solution can be computed analytically and the errors can be computed.

The exact solution here is given as,

$$
\begin{aligned}
\boldsymbol{u}\_{\mathrm{ex}} (x,y) &= -\frac{p_0}{k} \frac{\mathrm{H}^{(1)}\_{1}(k\sqrt{x^2+y^2})}{\mathrm{H}^{(1)}\_{0}(k)}\boldsymbol{e}\_{r}, \\
{p}\_{\mathrm{ex}} (x,y)  &= p\_0 \frac{\mathrm{H}^{(1)}\_{0}(k\sqrt{x^2+y^2})}{\mathrm{H}^{(1)}\_{0}(k)}.
\end{aligned}
$$

Errors are computed for a range of mesh resolutions for interpolated exact displacement field and FEM-solver computed displacement field.
Errors are also computed for different projection techniques proposed for one and two orders of approximation.

## Run

The following instructions run the example **assuming you are in the current directory**.

An automated script [`main.jl`](main.jl) is provided to run scripts in a particular order,

```bash
julia +1.11.0 --project=.. main.jl
```
This should generate the meshes of various resolutions in the `meshes` folder. All output fields are written within `.jld2` files and visualizations in `.pdf` files within the `results` folder.
