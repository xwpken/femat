# femat

FEMAT is a finite element analysis project based on code compiled during my master's studies.
It is used to implement:

- 2D linear elasticity
- 3D linear elasticity
- geometric nonlinearity in Total Lagrangian (TL) format

The core source code is MATLAB-based, and the project also provides a Python interface for integration workflows.

## Project description

This repository consolidates research and engineering code developed for graduate-level finite element method studies.
The goal is to provide a practical and reproducible implementation of linear-elastic and geometrically nonlinear solid mechanics analyses.

## Installation

### MATLAB

1. Clone this repository:
   ```bash
   git clone https://github.com/xwpken/femat.git
   ```
2. Open MATLAB and add the repository directory to your path.
3. Run your FEMAT scripts/functions from the MATLAB environment.

### Python interface

If you want to drive FEMAT workflows from Python, install Python and required dependencies in your environment, then configure access to MATLAB-based routines from your Python entry points provided in this repository.

## Citation

If this project helps your research or teaching, please cite it using the GitHub-repository BibTeX format:

```bibtex
@misc{femat_github,
  author       = {xwpken},
  title        = {femat},
  year         = {2026},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{https://github.com/xwpken/femat}}
}
```
