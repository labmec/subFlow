# subFlow

**subFlow** is a C++ library designed for the simulation of multiphase flows in porous media. It provides robust and efficient tools for researchers and engineers working on subsurface flow problems, such as groundwater hydrology, oil reservoir engineering, and environmental modeling.

## Features

- **Fluid Flow Analysis**  
    Utilizes a locally conservative mixed Finite Element Method (mixed-FEM) formulation to accurately solve fluid flow in porous media.

- **Saturation Transport**  
    Implements a Finite Volume (FV) scheme for the simulation of saturation transport, ensuring mass conservation and stability.

- **Coupled Analyses**  
    Supports fully coupled simulations using a Sequential Fully Implicit (SFI) method, enabling robust and efficient multiphase flow modeling.

## Key Capabilities

- Accurate simulation of multiphase flow in heterogeneous porous media
- Modular and extensible C++ codebase
- Designed for high performance and scalability

## Getting Started

To use subFlow in your project, clone the repository and follow the build instructions in the [INSTALL.md](INSTALL.md) file.

```bash
git clone https://github.com/labmec/subFlow.git
cd subFlow
# Follow build instructions
```

## Documentation

Comprehensive documentation and example problems are available in the `docs/` directory.

## License

subFlow is released under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

For questions, bug reports, or contributions, please open an issue or contact the maintainers.

