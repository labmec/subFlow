## subFlow library

subFlow is a lightweight C++ library for analyzing multiphase flow in porous media with three phases: water, oil and gas. It couples a locally conservative Darcy solver and a saturation transport solver to simulate phase displacement and transport in heterogeneous reservoirs.

Key features:
- Multiphase formulation for water–oil–gas systems.
- Darcy flow solver (compressible and incompressible options) using a locally conservative finite element formulation with H(div)–L2 approximation pairs for total flux and pressure.
- Saturation transport solver based on a Finite Volume scheme (robust upwinding and conservation), with gravitational segregation using an Implicit Hybrid Upwind (IHU) strategy.
- Coupling via a Sequential Fully Implicit (SFI) method to handle strong nonlinear coupling between flow and transport.
- Flexible input/output: supports parameterized .json inputs and mesh files (.geo, .msh) used by the repository.
- Designed for extension: modular solver components and clear interfaces for physics, discretization and coupling strategies.

Typical workflow:
1. Build mesh and provide rock/fluids data (input/).
2. Solve the Darcy problem for pressure and total flux (H(div)–L2 FEM).
3. Advance saturations with the Finite Volume transport solver.
4. Iterate SFI steps until convergence for each time step.

Use subFlow when you need a conservative, modular framework for reservoir-scale multiphase simulations with emphasis on accurate flux representation and robust transport coupling.

## Repository Structure

The repository is organized as follows:

```
wann/
├── docs/
├── input/
├── mathematica/
├── src/
└── targets/
```

### Directory Details

- **`docs/`**: Documentation files, including guides and references for using the repository.
- **`input/`**: Input files (.json) and meshe files (.geo and .msh) for running the C++ simulations.
- **`mathematica/`**: Contains Mathematica notebooks and scripts for symbolic computations related to unidimensional multiphase flow simulation.
- **`src/`**: Core C++ source code for the repository.
- **`targets/`**: Where the executables with the main functions are located.

Contributions and suggestions for improvement are encouraged!