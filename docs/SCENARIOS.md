# Simulation Scenarios and Benchmarks Guide

`wake` includes five built-in aerodynamic and hydrodynamic simulation scenarios that benchmark the ISPC computational kernels against canonical analytical and experimental fluid mechanics data.

---

## 1. Scenario Catalog Summary

| Scenario Identifier | Dimension | Lattice Model | Collision Operator | Key Physics |
|:---|:---:|:---:|:---:|:---|
| `karman` | 2D | D2Q9 | TRT / MRT | Unsteady vortex shedding, lift oscillation, FFT Strouhal frequency |
| `airfoil` | 2D | D2Q9 | TRT / MRT | NACA 0012 dynamic stall, leading-edge vortex rollup, drag polar |
| `cavity` | 3D | D3Q19 | MRT | Cubic lid-driven circulation, secondary corner vortices, Ghia profile |
| `sphere` | 3D | D3Q19 | Smagorinsky LES | Subgrid turbulence, hairpin vortex loops, Q-criterion iso-surfaces |
| `thermal` | 2D | D2Q9 + ADRE | TRT + Boussinesq | Rayleigh-Benard buoyancy rolls, plume dynamics, Nusselt scaling |

---

## 2. Scenario 1: Karman Vortex Street (`karman`)

### 2.1 Problem Specification
A circular cylinder of diameter $D = 20$ lattice units is placed symmetrically within a channel of dimensions $NX = 400, NY = 100$ at coordinate $(X_c, Y_c) = (80, 50)$.

```
      y = NY-1  ┌─────────────────────────────────────────────────────────┐  No-Slip Wall
                │                                                         │
Inflow:         │ ──►                                                     │  Outflow:
Zou-He Velocity │ ──►               ● Cylinder (D=20)                     │  Zou-He Pressure
u_x = 0.08      │ ──►               (x=80, y=50)                          │  rho = 1.0
                │                                                         │
      y = 0     └─────────────────────────────────────────────────────────┘  No-Slip Wall
                x = 0                                                     x = NX-1
```

* **Inflow:** Zou-He parabolic velocity profile with peak centerline velocity $U_{max} = 0.08$.
* **Outflow:** Zou-He constant pressure boundary ($\rho = 1.0$).
* **Top/Bottom:** No-slip halfway bounce-back walls.
* **Reynolds Number:**
  $$Re = \frac{U_{\infty} D}{\nu} = 100$$
  At $Re = 100$, laminar boundary layers separate alternately from the top and bottom of the cylinder, producing a self-sustaining periodic Karman vortex street.

### 2.2 Aerodynamic Telemetry
The Momentum Exchange Method computes instantaneous drag coefficient $C_D(t)$ and lift coefficient $C_L(t)$ every time step. An automated radix-2 FFT extracts the dominant shedding frequency $f_{shed}$ to evaluate the non-dimensional Strouhal number:

$$St = \frac{f_{shed} D}{U_{\infty}} \approx 0.165 \pm 0.005$$

### 2.3 Running Scenario 1
```bash
./build/bin/wake run --scenario karman --steps 4000 --output output/karman
python3 scripts/plot_telemetry.py --csv output/karman/telemetry.csv --out output/karman/forces.png
python3 scripts/verify_strouhal.py --csv output/karman/telemetry.csv --re 100
```

---

## 3. Scenario 2: NACA 0012 Airfoil Stall (`airfoil`)

### 3.1 Problem Specification
A symmetrical NACA 0012 airfoil with chord length $C = 60$ lattice units is immersed in a flow domain at an angle of attack $\alpha \in [0^{\circ}, 25^{\circ}]$.

The airfoil profile is generated analytically from the 4-digit thickness distribution:

$$y_t(x) = 5 t C \left[ 0.2969 \sqrt{\frac{x}{C}} - 0.1260 \left(\frac{x}{C}\right) - 0.3516 \left(\frac{x}{C}\right)^2 + 0.2843 \left(\frac{x}{C}\right)^3 - 0.1015 \left(\frac{x}{C}\right)^4 \right]$$

with $t = 0.12$. Coordinate clustering uses a non-linear cosine discretization to pack points tightly around the high-curvature leading edge:

$$\frac{x_k}{C} = \frac{1}{2} \left[ 1 - \cos\left(\frac{k \pi}{N_p}\right) \right], \quad k = 0, \dots, N_p$$

### 3.2 Physics of Dynamic Stall
* At low angles ($\alpha < 8^{\circ}$), flow remains attached with linear lift slope $dC_L / d\alpha \approx 2\pi$.
* As $\alpha$ exceeds the critical stall angle ($\alpha \approx 14^{\circ}$), adverse pressure gradients trigger boundary layer separation near the leading edge.
* A coherent leading-edge vortex (LEV) rolls up, traverses the suction surface, and detaches, causing a sudden loss of lift and a sharp rise in form drag.

### 3.3 Running Scenario 2
```bash
./build/bin/wake run --scenario airfoil --steps 3000 --output output/airfoil
```

---

## 4. Scenario 3: 3D Cubic Lid-Driven Cavity (`cavity`)

### 4.1 Problem Specification
A cubic enclosure of size $N \times N \times N$ (e.g. $64 \times 64 \times 64$) contains a fluid initially at rest. Five faces (bottom, left, right, front, back) are rigid no-slip walls. The top lid ($y = NY - 1$) translates horizontally with constant tangential velocity $\mathbf{u} = (U_0, 0, 0)$ with $U_0 = 0.1$.

### 4.2 Benchmark Validation
This setup benchmarks the D3Q19 Multiple-Relaxation-Time (MRT) collision operator against the classical numerical solutions of Ghia et al. (1982).
* **Primary Vortex:** Centered in the upper-right quadrant of the cavity.
* **Secondary Corner Eddies:** Bottom-left and bottom-right corners exhibit counter-rotating secondary recirculation zones.
* **Profile Extraction:** The simulation automatically extracts and outputs normalized centerline velocity profiles $u_x(y)$ at $x = NX/2, z = NZ/2$ and $u_y(x)$ at $y = NY/2, z = NZ/2$.

### 4.3 Running Scenario 3
```bash
./build/bin/wake run --scenario cavity --steps 2500 --output output/cavity
```

---

## 5. Scenario 4: 3D Sphere Wake with Smagorinsky LES (`sphere`)

### 5.1 Problem Specification
A sphere of radius $R = 12$ ($D = 24$) is placed in a 3D channel of dimensions $160 \times 64 \times 64$. The flow operates at Reynolds number $Re = 300$ to $1000$.

### 5.2 Smagorinsky Subgrid Modeling and Q-Criterion
At $Re > 300$, the wake loses planar symmetry, shedding complex hairpin-shaped vortex loops. To resolve subgrid turbulence dissipation without artificial damping, the Smagorinsky LES kernel dynamically evaluates local filtered strain rates:

$$\tau_{eff} = \tau_0 + \tau_t$$

Coherent vortex cores are identified using Hunt's second invariant of the velocity gradient tensor ($Q$-criterion):

$$Q = \frac{1}{2} \left( \|\boldsymbol{\Omega}\|^2 - \|\mathbf{S}\|^2 \right) > 0$$

where $\boldsymbol{\Omega} = \frac{1}{2}(\nabla \mathbf{u} - \nabla \mathbf{u}^T)$ is the vorticity tensor and $\mathbf{S} = \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T)$ is the strain rate tensor. Regions with $Q > 0$ delineate regions where rotation dominates strain, isolating vortex tubes cleanly.

### 5.3 Running Scenario 4
```bash
./build/bin/wake run --scenario sphere --steps 2000 --output output/sphere
python3 scripts/plot_paraview_pipeline.py --vti output/sphere/wake_002000.vti
```

---

## 6. Scenario 5: Rayleigh-Benard Natural Thermal Convection (`thermal`)

### 6.1 Problem Specification
A 2D rectangular cavity of aspect ratio $L/H = 2$ is heated from below and cooled from above:
* Bottom plate ($y = 0$): Fixed hot temperature $T_{hot} = 1.0$.
* Top plate ($y = NY - 1$): Fixed cold temperature $T_{cold} = 0.0$.
* Vertical side walls ($x = 0, x = NX - 1$): Adiabatic ($\partial T / \partial x = 0$) no-slip walls.

Buoyancy coupling is enforced via Boussinesq body forcing:

$$\mathbf{F}_{buoyancy} = \rho \beta (T - T_0) g \hat{\mathbf{y}}$$

### 6.2 Convective Instability and Nusselt Scaling
Below the critical Rayleigh number $Ra_c \approx 1708$, heat transfer is purely conductive (linear temperature gradient, $Nu = 1.0$). Above $Ra_c$, buoyancy overcomes viscous resistance, establishing stable counter-rotating convection rolls.

The engine evaluates the space-averaged Nusselt number at the hot boundary:

$$Nu = \frac{H}{\Delta T} \left. \frac{\partial \langle T \rangle}{\partial y} \right|_{y=0}$$

### 6.3 Running Scenario 5
```bash
./build/bin/wake run --scenario thermal --steps 3000 --output output/thermal
```
