# DDM_release: Discrete Damage Mechanics (DDM6)

Discrete Damage Mechanics (DDM6) software distribution repository. This package provides active source code, build scripts, libraries, and examples for damage prediction and material property identification for laminated composites.

## Repository Structure

- **`src/`**: Core Fortran source files for the active DDM6 implementation.
- **`lib/`**: Compiled modules and library assets.
- **`examples/`**: Verification and application examples.
- **`legacy/`**: Historical and archival baseline files (such as the 2017.obj library/plugin for Abaqus) preserved for backward compatibility and reference. 

## Active Distribution
The active codebase is built around **DDM6-gfortran**. Future releases for other compiler environments (such as Intel Fortran and Abaqus user subroutines) will be organized into dedicated paths or distinct release tags to prevent version confusion.

---

## License & Attribution
This project includes original work by Ever Barbero. Copyright (c) 2026 Ever Barbero.

Redistributions of this software, in source or binary form, must retain the above attribution in any documentation, about dialog, or other materials provided with the distribution.

Public, academic, or scientific use of this software should credit the original author:
* **Ever Barbero** 
* ORCID [0000-0001-9121-6277](https://github.com/everbarbero)