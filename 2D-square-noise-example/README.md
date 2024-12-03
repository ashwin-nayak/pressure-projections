# Unit Square Domain with Noisy Data

2D Plane wave propagation with noise in unit square domain.

- Domain, $\Omega=(0,1) \times (0,1)$
- Boundary Conditions
  - Rigid boundary, $\boldsymbol{u} \cdot \boldsymbol{n} = 0$ at $x=0$ and $x=1$,
  - Pressure $p=1.0$ at $y=0$,
  - Impedance $Z=1.0$ at $y=1$.

The boundary conditions are expressly chosen such that the exact solution can be computed analytically and the errors can be computed.

The exact solution here is given as,

$$
\begin{aligned}
    \boldsymbol{u}\_{\mathrm{ex}} (x,y) &= {\frac{\mathrm{i} {p}\_{0}}{k}} \exp (\mathrm{i} k x) \boldsymbol{e}\_{y}, \\
    {p}\_{\mathrm{ex}} (x,y)  &= {p}\_{0} \exp(\mathrm{i} k y).
\end{aligned}
$$

However, a noise field $\chi_{h}(x,y)$ is added into this displacement field when included, i.e.
$$
\boldsymbol{u}_\mathrm{SNR}(x,y) = \boldsymbol{u}_\mathrm{ex}(x,y)\left(1+\frac{\chi_{h}(x,y)}{\mathrm{SNR}}\right).
$$

Here, $\chi_{h}(x,y)$ is the interpolated Gaussian noise $\mathcal{N}(0,1)$, while the parameter $\mathrm{SNR}$ is the Signal-to-Noise ratio and varied between different runs.

Errors are computed for a range of mesh resolutions for interpolated exact displacement field.
Errors are also computed for different projection techniques proposed for one and two orders of approximation.

## Run

The following instructions run the example **assuming you are in the current directory**.

An automated script [`main.jl`](main.jl) is provided to run scripts in a particular order,

```bash
julia +1.11.0 --project=.. main.jl
```
This should generate the meshes of various resolutions in the `meshes` folder. All output fields are written within `.jld2` files and visualizations in `.pdf` files within the `results` folder.
