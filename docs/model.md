# Model description

The Biogeochemical Reaction Network Simulator (BRNS) represents biogeochemical and reactive transport dynamics in porous media by solving a set of coupled, one-dimensional mass conservation equations for dissolved and solid species along a depth-resolved domain. Transport of solutes and particles is described through a combination of molecular diffusion, advection (e.g., sediment burial or groundwater flow), bioturbation, and bioirrigation, while local concentration changes due to biogeochemical transformations are represented by a user-defined network of kinetic and equilibrium reactions.

## Governing equation

The stand-alone solver represents one-dimensional, multi-component reactive transport along the depth coordinate $x$. For each transported component $k$, the transient code advances the effective concentration form

$$
\frac{\partial C_k}{\partial t} = -v_k\frac{\partial C_k}{\partial x} + \frac{\partial}{\partial x}\left(D_k\frac{\partial C_k}{\partial x}\right) + R_k(\mathbf{C},x,t).
$$

Here, $C_k$ is the component concentration, $v_k$ is its effective advective velocity, $D_k$ is its total dispersion coefficient, and $R_k$ is the net production rate defined by the generated reaction network. The model distinguishes dissolved and solid components. Dissolved components use $v_k=v_d=w+q/(\phi A)$; solid components use $v_k=v_s=w$. The transport step uses the spatially variable $v_k$ and $D_k$ specified at grid faces. Components marked as non-transported are excluded from the transport step.

For dissolved components, the implemented coefficient is

$$
D_k = D_{\mathrm{mol},k} + \alpha_L\left|\frac{q}{\phi A}\right| + D_b,
$$

whereas solid components use $D_k=D_b$. Molecular diffusion is specified directly as an input parameter for each dissolved component.

| Symbol | Fortran quantity | Definition |
| --- | --- | --- |
| $C_k$ | `sp(k,j)` | Concentration of component $k$ at grid point $j$. |
| $R_k$ | generated `residual.f` / `rates.f` | Net component rate from the reaction network; it may depend nonlinearly on all concentrations and kinetic parameters. |
| $v_d$ | `vd(j)` | Dissolved-component velocity, $w+q/(\phi A)$. |
| $v_s$ | `vs(j)` | Solid-component velocity, $w$. |
| $D_k$ | `disp(k,j)` | Total component dispersion coefficient. |
| $D_{\mathrm{mol},k}$ | `dsol_0(k)` | Molecular diffusion coefficient of a dissolved component. |
| $\alpha_L$ | `aL` | Longitudinal dispersivity; defines mechanical dispersion. |
| $D_b$ | `Db0` / `db` | Bioturbation coefficient, optionally depth dependent. |
| $w$ | `w0` / `w` | Solid-matrix advection velocity, optionally depth dependent. |
| $q$ | `q0` | Volumetric fluid-flow parameter. |
| $\phi$ | `por(j)` | Porosity; used to convert $q$ to pore-water velocity and in flux boundary conditions. |
| $A$ | `area(j)` | Cross-sectional area; used with porosity to convert $q$ to pore-water velocity. |
| $\Delta t$ | `delt` | Time-step length. |

## Initial and boundary conditions

The initial condition is prescribed independently for each component,

$$
C_k(x,0)=C_{k,0}(x).
$$

The code supports a spatially uniform initial value, a profile read from `initialconc.txt`, and model-specific input profiles. Concentrations below $10^{-20}$ are raised to this value during initialization.

At both domain boundaries, each component may be assigned one of the following boundary conditions:

| Condition | Mathematical form |
| --- | --- |
| Fixed concentration | $C_k=C_{k,b}$ |
| Prescribed concentration gradient | $\partial C_k/\partial x=g_{k,b}$ |
| Prescribed total advective-diffusive flux | $J_k=J_{k,b}$ |

For a fixed-concentration boundary, the specified concentration is retained during both the transport and reaction calculations.

## Numerical solution

Transport and reaction are coupled by sequential operator splitting. Each time step first calculates transport for every component and then calculates the reaction network separately at every concentration location. The local nonlinear reaction system is solved by Newton iteration with relaxation. Each Newton correction is obtained by LU decomposition of the reaction Jacobian. The initial time-step length is read from the model input and is subsequently adapted from the temporal curvature of a selected master component.

Transport is discretized implicitly with finite differences. Each equation couples a concentration point to its two neighbours. This produces a tridiagonal linear system, which is solved by the Thomas algorithm for tridiagonal systems.

The one-dimensional grid may be non-uniform. Concentrations are stored at the odd-numbered grid locations. Even-numbered locations lie between adjacent concentration locations and define the interfaces. The two end concentration locations form the model boundaries. Interface positions are used to calculate the local distances between concentration locations.
