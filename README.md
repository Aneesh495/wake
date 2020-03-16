# wake

A high-throughput computational aerodynamics and Lattice Boltzmann hydrodynamics engine vectorized in ISPC with a C++ host runtime, multi-relaxation-time collision, Smagorinsky LES turbulence, immersed boundary method, and real-time wake visualization.

## Key Subsystems

- Core ISPC SPMD Kernels: SIMD-accelerated D2Q9, D3Q19, and D3Q27 streaming and collision.
- Collision Dynamics: BGK, Two-Relaxation-Time (TRT), and Multiple-Relaxation-Time (MRT) operators.
- Turbulence Modeling: Smagorinsky and Wall-Adapting Local Eddy-viscosity (WALE) Large Eddy Simulation.
- Fluid-Structure Interaction: Diffuse-interface Immersed Boundary Method (IBM) for arbitrary solid meshes.
- Aerodynamic Telemetry: Momentum Exchange Method for instantaneous lift, drag, and Strouhal vortex shedding analysis.
- Live Real-Time Rendering: ANSI TrueColor terminal stream and VTK binary volumetric export for ParaView.
