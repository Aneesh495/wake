# Mathematical and Numerical Theory Guide

This document presents the theoretical foundations, multiscale asymptotic analysis, collision operators, turbulence subgrid models, and boundary mechanics implemented in `wake`.

---

## 1. Kinetic Foundations: From Boltzmann to Lattice

The continuum Boltzmann equation describes the statistical evolution of a single-particle phase-space distribution function $f(\mathbf{x}, \boldsymbol{\xi}, t)$:

$$\frac{\partial f}{\partial t} + \boldsymbol{\xi} \cdot \nabla_{\mathbf{x}} f + \frac{\mathbf{F}}{\rho} \cdot \nabla_{\boldsymbol{\xi}} f = \Omega(f)$$

where $\boldsymbol{\xi}$ is the microscopic particle velocity, $\mathbf{F}$ is an external force density, and $\Omega(f)$ is the collision operator modeling intermolecular momentum exchange.

In the Lattice Boltzmann Method (LBM), velocity space is discretized into a finite set of $Q$ discrete velocity vectors $\{\mathbf{c}_i\}_{i=0}^{Q-1}$ that match the spatial lattice symmetries. Integrating along characteristics over time step $\Delta t$ yields the discrete lattice Boltzmann equation:

$$f_i(\mathbf{x} + \mathbf{c}_i \Delta t, t + \Delta t) - f_i(\mathbf{x}, t) = \Omega_i(\mathbf{x}, t) + S_i(\mathbf{x}, t) \Delta t$$

where:
* $f_i(\mathbf{x}, t)$ is the probability density of particles moving with lattice velocity $\mathbf{c}_i$ at position $\mathbf{x}$ and time $t$.
* $\Omega_i$ is the discrete collision operator.
* $S_i$ is the discrete forcing term accounting for external body forces.

The lattice speed of sound $c_s$ is determined by the lattice geometry:

$$c_s = \frac{1}{\sqrt{3}} \frac{\Delta x}{\Delta t}$$

---

## 2. Maxwell-Boltzmann Equilibrium via Hermite Expansions

The continuous Maxwell-Boltzmann equilibrium distribution in a thermal equilibrium state at temperature $T_0$ is:

$$f^{(0)}(\mathbf{x}, \boldsymbol{\xi}, t) = \frac{\rho}{(2\pi c_s^2)^{D/2}} \exp\left(-\frac{|\boldsymbol{\xi} - \mathbf{u}|^2}{2 c_s^2}\right)$$

In `wake`, the discrete equilibrium distribution $f_i^{eq}$ is derived by expanding the continuous Maxwellian in a truncated series of orthogonal Hermite polynomials $\mathcal{H}^{(n)}(\boldsymbol{\xi} / c_s)$ using Gauss-Hermite quadrature:

$$f_i^{eq}(\rho, \mathbf{u}) = w_i \rho \left[ 1 + \frac{\mathbf{c}_i \cdot \mathbf{u}}{c_s^2} + \frac{(\mathbf{c}_i \cdot \mathbf{u})^2 - c_s^2 |\mathbf{u}|^2}{2 c_s^4} + \mathcal{O}(|\mathbf{u}|^3) \right]$$

where $w_i$ are lattice quadrature weights satisfying the moment conservation relations:

$$\sum_{i=0}^{Q-1} w_i = 1, \quad \sum_{i=0}^{Q-1} w_i c_{i\alpha} c_{i\beta} = c_s^2 \delta_{\alpha\beta}, \quad \sum_{i=0}^{Q-1} w_i c_{i\alpha} c_{i\beta} c_{i\gamma} c_{i\delta} = c_s^4 (\delta_{\alpha\beta}\delta_{\gamma\delta} + \delta_{\alpha\gamma}\delta_{\beta\delta} + \delta_{\alpha\delta}\delta_{\beta\gamma})$$

### 2.1 Quadrature Weights

* **D2Q9 Lattice:**
  * Center ($i = 0$): $w_0 = 4/9$
  * Axis-aligned ($i = 1, 2, 3, 4$): $w_i = 1/9$
  * Diagonal ($i = 5, 6, 7, 8$): $w_i = 1/36$

* **D3Q19 Lattice:**
  * Center ($i = 0$): $w_0 = 1/3$
  * Axis-aligned ($i = 1 \dots 6$): $w_i = 1/18$
  * Face diagonals ($i = 7 \dots 18$): $w_i = 1/36$

---

## 3. Multiscale Chapman-Enskog Asymptotic Expansion

The recovery of the Navier-Stokes equations from the discrete kinetic equation is established via multiscale Chapman-Enskog expansion. Introducing a small Knudsen number parameter $\epsilon \ll 1$:

$$f_i = f_i^{(0)} + \epsilon f_i^{(1)} + \epsilon^2 f_i^{(2)} + \mathcal{O}(\epsilon^3)$$

Differential operators are expanded across multi-scale temporal and spatial regimes:

$$\frac{\partial}{\partial t} = \epsilon \frac{\partial}{\partial t_1} + \epsilon^2 \frac{\partial}{\partial t_2}, \quad \nabla = \epsilon \nabla_1$$

Taylor expanding the left hand side of the discrete streaming operator:

$$f_i(\mathbf{x} + \mathbf{c}_i \Delta t, t + \Delta t) - f_i(\mathbf{x}, t) = \sum_{n=1}^{\infty} \frac{\Delta t^n}{n!} D_i^n f_i, \quad D_i \equiv \frac{\partial}{\partial t} + \mathbf{c}_i \cdot \nabla$$

### 3.1 Order $\epsilon^1$ (Euler Equations)

Matching terms at order $\epsilon$:

$$D_{1i} f_i^{(0)} = -\frac{1}{\tau \Delta t} f_i^{(1)}$$

Summing over all velocities $i$ yields the continuity equation at time scale $t_1$:

$$\frac{\partial \rho}{\partial t_1} + \nabla_1 \cdot (\rho \mathbf{u}) = 0$$

Multiplying by $\mathbf{c}_i$ and summing yields the inviscid Euler momentum equation:

$$\frac{\partial (\rho \mathbf{u})}{\partial t_1} + \nabla_1 \cdot \boldsymbol{\Pi}^{(0)} = 0, \quad \boldsymbol{\Pi}^{(0)} = p \mathbf{I} + \rho \mathbf{u} \mathbf{u}$$

where the pressure obeys the ideal isothermal equation of state:

$$p = c_s^2 \rho$$

### 3.2 Order $\epsilon^2$ (Navier-Stokes Dissipation)

Matching terms at order $\epsilon^2$:

$$\frac{\partial f_i^{(0)}}{\partial t_2} + D_{1i} f_i^{(1)} + \frac{\Delta t}{2} D_{1i}^2 f_i^{(0)} = -\frac{1}{\tau \Delta t} f_i^{(2)}$$

Evaluating momentum moments yields the viscous stress contribution:

$$\frac{\partial (\rho \mathbf{u})}{\partial t_2} + \nabla_1 \cdot \left[ \left(1 - \frac{1}{2\tau}\right) \boldsymbol{\Pi}^{(1)} \right] = 0$$

Recombining time scales $t = \epsilon t_1 + \epsilon^2 t_2$ yields the incompressible Navier-Stokes momentum equation:

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \mathbf{u}) = -\nabla p + \nabla \cdot \left[ \rho \nu \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right) \right] + \mathbf{F}$$

The kinematic shear viscosity $\nu$ is directly tied to the relaxation time $\tau$:

$$\nu = c_s^2 \left( \tau - \frac{1}{2} \right) \Delta t$$

---

## 4. Collision Operators

### 4.1 Bhatnagar-Gross-Krook (BGK)

The single-relaxation-time collision operator relaxes all non-equilibrium moments at the identical rate:

$$\Omega_i^{BGK} = -\frac{1}{\tau} \left( f_i - f_i^{eq} \right)$$

While computationally efficient, BGK encounters unphysical wall slip and numerical instability when $\tau \to 0.5$ (very low viscosity).

### 4.2 Two-Relaxation-Time (TRT)

The TRT operator decomposes distribution functions into symmetric ($+$) and anti-symmetric ($-$) components with respect to velocity inversion ($\mathbf{c}_{\bar{i}} = -\mathbf{c}_i$):

$$f_i^+ = \frac{f_i + f_{\bar{i}}}{2}, \quad f_i^- = \frac{f_i - f_{\bar{i}}}{2}$$

The collision step relaxes symmetric moments at rate $\omega_s = 1/\tau_s$ and anti-symmetric moments at rate $\omega_a = 1/\tau_a$:

$$f_i^{coll} = f_i - \omega_s \left( f_i^+ - f_i^{eq,+} \right) - \omega_a \left( f_i^- - f_i^{eq,-} \right)$$

* $\tau_s$ controls the physical kinematic shear viscosity: $\nu = c_s^2 (\tau_s - 0.5) \Delta t$.
* $\tau_a$ controls acoustic and numerical damping.

Setting the TRT "magic parameter":

$$\Lambda \equiv \left(\tau_s - \frac{1}{2}\right)\left(\tau_a - \frac{1}{2}\right) = \frac{1}{4}$$

eliminates numerical slip at bounce-back walls and locks the effective boundary location exactly halfway between grid points for arbitrary viscosity.

### 4.3 Multiple-Relaxation-Time (MRT)

The MRT model transforms discrete distributions into physical moment space using an orthogonal transformation matrix $\mathbf{M}$:

$$\mathbf{m} = \mathbf{M} \mathbf{f}, \quad \mathbf{f} = \mathbf{M}^{-1} \mathbf{m}$$

In moment space, collision is executed using a diagonal relaxation matrix $\mathbf{S}$:

$$\mathbf{m}^{coll} = \mathbf{m} - \mathbf{S} (\mathbf{m} - \mathbf{m}^{eq})$$

For D2Q9, the 9 orthogonal moments are:
1. Density $\rho$ ($s_0 = 0$)
2. Energy $e$ ($s_1 = 1.63$)
3. Energy squared $\epsilon$ ($s_2 = 1.14$)
4. Momentum $j_x$ ($s_3 = 0$)
5. Energy flux $q_x$ ($s_4 = 1.92$)
6. Momentum $j_y$ ($s_5 = 0$)
7. Energy flux $q_y$ ($s_6 = 1.92$)
8. Normal stress $p_{xx}$ ($s_7 = 1/\tau$)
9. Shear stress $p_{xy}$ ($s_8 = 1/\tau$)

Decoupling the relaxation of ghost and energy flux modes from the stress modes provides unconditional linear stability across high Reynolds number regimes.

---

## 5. Body Forcing: Guo Scheme

Applying an external body force density $\mathbf{F}$ in LBM requires correcting the discrete velocity moments to avoid discrete lattice artifacts:

$$\rho \mathbf{u} = \sum_{i=0}^{Q-1} \mathbf{c}_i f_i + \frac{\Delta t}{2} \mathbf{F}$$

The discrete source term $S_i$ introduced by Guo et al. is:

$$S_i = \left( 1 - \frac{1}{2\tau} \right) w_i \left[ \frac{\mathbf{c}_i - \mathbf{u}}{c_s^2} + \frac{(\mathbf{c}_i \cdot \mathbf{u}) \mathbf{c}_i}{c_s^4} \right] \cdot \mathbf{F}$$

This formulation guarantees that the recovered Navier-Stokes momentum equation matches the continuous continuum equations to second-order spatial and temporal accuracy with zero spurious velocity dependencies.

---

## 6. Turbulence: Smagorinsky Large Eddy Simulation (LES)

In turbulent flows, large energetic eddies are resolved directly on the computational lattice, while subgrid-scale (SGS) dissipation is modeled through turbulent eddy viscosity $\nu_t$.

The filtered non-equilibrium stress tensor is evaluated directly from local non-equilibrium populations without finite difference gradients:

$$\Pi_{\alpha\beta}^{(neq)} = \sum_{i=0}^{Q-1} c_{i\alpha} c_{i\beta} \left( f_i - f_i^{eq} \right)$$

The filtered strain rate tensor $S_{\alpha\beta}$ is:

$$S_{\alpha\beta} = -\frac{1}{2 \rho c_s^2 \tau_{eff}} \Pi_{\alpha\beta}^{(neq)}$$

Using the Smagorinsky closure with Smagorinsky constant $C_s \approx 0.12$:

$$\nu_t = (C_s \Delta x)^2 |\bar{S}|, \quad |\bar{S}| = \sqrt{2 S_{\alpha\beta} S_{\alpha\beta}}$$

The effective relaxation time $\tau_{eff} = \tau_0 + \tau_t$ is determined by solving the quadratic equation:

$$\tau_{eff} = \tau_0 + \frac{1}{2} \left[ \sqrt{\tau_0^2 + \frac{2 \sqrt{2} (C_s \Delta x)^2}{\rho c_s^4} |\Pi^{(neq)}|} - \tau_0 \right]$$

This subgrid model dynamically provides local numerical stability in regions of intense shear and vortex rollup without contaminating laminar boundary layers.

---

## 7. Curved Boundaries: Bouzidi-Firdaouss-Lallemand (BFL)

For solid obstacles with curved or slanted boundaries, standard stair-step bounce-back introduces first-order geometric errors. The BFL method interpolates incoming distributions using the exact geometric distance fraction $q$:

$$q = \frac{|\mathbf{x}_f - \mathbf{x}_w|}{|\mathbf{x}_f - \mathbf{x}_b|}, \quad q \in [0, 1)$$

where $\mathbf{x}_f$ is the fluid node, $\mathbf{x}_b$ is the solid node, and $\mathbf{x}_w$ is the physical wall intersection.

```
 Fluid Node (xf)          Wall (xw)         Solid Node (xb)
      ○ ───────────────────── ■ ──────────────────── ●
      │ <───── q Δx ────────> │ <── (1-q) Δx ─────> │
```

Depending on $q$, the reflected distribution $f_{\bar{i}}(\mathbf{x}_f, t + \Delta t)$ is interpolated:

* **For $q < 1/2$:**
  $$f_{\bar{i}}(\mathbf{x}_f, t + \Delta t) = 2q f_i(\mathbf{x}_f, t) + (1 - 2q) f_i(\mathbf{x}_f - \mathbf{c}_i \Delta t, t)$$

* **For $q \ge 1/2$:**
  $$f_{\bar{i}}(\mathbf{x}_f, t + \Delta t) = \frac{1}{2q} f_i(\mathbf{x}_f, t) + \frac{2q - 1}{2q} f_{\bar{i}}(\mathbf{x}_f, t)$$

The BFL scheme maintains strict second-order spatial convergence across arbitrary geometries.

---

## 8. Immersed Boundary Method (IBM)

For moving or compliant boundaries, the Immersed Boundary Method decouples the Eulerian Cartesian fluid grid from Lagrangian boundary markers $\mathbf{X}_k$.

### 8.1 Regularized Dirac Delta Function

Information transfer between Eulerian coordinates $\mathbf{x}$ and Lagrangian coordinates $\mathbf{X}$ is mediated by Peskin's 4-point smoothed kernel:

$$\delta_h(\mathbf{x} - \mathbf{X}) = \frac{1}{h^D} \phi\left(\frac{x - X}{h}\right) \phi\left(\frac{y - Y}{h}\right) \phi\left(\frac{z - Z}{h}\right)$$

where the 1D kernel $\phi(r)$ satisfies:

$$\phi(r) = \begin{cases}
\frac{1}{8}(3 - 2|r| + \sqrt{1 + 4|r| - 4r^2}), & 0 \le |r| \le 1 \\
\frac{1}{8}(5 - 2|r| - \sqrt{-7 + 12|r| - 4r^2}), & 1 < |r| \le 2 \\
0, & |r| > 2
\end{cases}$$

### 8.2 Velocity Interpolation and Force Spreading

1. **Velocity Interpolation:** The fluid velocity at marker position $\mathbf{X}_k$ is interpolated:
   $$\mathbf{U}_k = \sum_{\mathbf{x}} \mathbf{u}(\mathbf{x}) \delta_h(\mathbf{x} - \mathbf{X}_k) \Delta x^D$$

2. **Lagrangian Boundary Force:** Restoring forces enforce no-slip:
   $$\mathbf{F}_k = \frac{\mathbf{U}_k^{target} - \mathbf{U}_k}{\Delta t}$$

3. **Eulerian Force Spreading:** The Lagrangian forces are spread to adjacent Eulerian fluid nodes:
   $$\mathbf{f}(\mathbf{x}) = \sum_{k} \mathbf{F}_k \delta_h(\mathbf{x} - \mathbf{X}_k) \Delta s_k$$

---

## 9. Thermal Convection: Boussinesq Coupling

Thermal flow is governed by the coupled Navier-Stokes and Advection-Diffusion Reaction Equation (ADRE) for temperature $T$:

$$\frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T = \alpha \nabla^2 T$$

where $\alpha$ is thermal diffusivity.

### 9.1 Boussinesq Approximation

Density variations are assumed negligible except in the gravity body force term:

$$\mathbf{F}_{buoyancy} = -\rho \beta (T - T_0) \mathbf{g}$$

where $\beta$ is the volumetric thermal expansion coefficient and $T_0$ is the reference temperature.

### 9.2 Non-Dimensional Numbers

* **Rayleigh Number ($Ra$):** Ratio of buoyancy forces to viscous and thermal dissipation:
  $$Ra = \frac{g \beta (T_{hot} - T_{cold}) H^3}{\nu \alpha}$$

* **Prandtl Number ($Pr$):** Ratio of momentum diffusivity to thermal diffusivity:
  $$Pr = \frac{\nu}{\alpha}$$

* **Nusselt Number ($Nu$):** Ratio of total heat transfer to conductive heat transfer:
  $$Nu = \frac{H}{\Delta T} \left. \frac{\partial \langle T \rangle}{\partial y} \right|_{wall}$$

Natural convection rolls spontaneously form above the critical Rayleigh number $Ra_c \approx 1708$ for plates with rigid boundaries.
