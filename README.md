# 2D Incompressible CFD Simulator

This is a full-featured 2D incompressible flow simulator using the vorticity-streamfunction formulation, written in Python and Fortran. It supports a GUI-based parameter setup, flow around a circular obstacle, multiple solver choices, and full field visualizations. The simulation is powered by a Fortran backend and a Python frontend that handles user input and result visualization.

## Features

- Vorticity-streamfunction formulation of 2D incompressible Navier-Stokes
- Built-in GUI for parameter selection (domain, resolution, viscosity, solvers)
- Support for optional cylindrical obstruction
- Time integration using multiple Runge-Kutta methods (Euler, Heun, RK4, RK6, RK8)
- Poisson solver choices: Jacobi, Gauss-Seidel, SOR
- Field visualizations:
  - Velocity magnitude and direction
  - Streamfunction contours
  - Vorticity distribution
  - Relative pressure field
- Scalar plots over time:
  - Lift and drag coefficients (Cl, Cd)
  - Pressure drop
  - Enstrophy, circulation, kinetic energy
  - Mass flow rates

## Folder Structure

```
.
├── main.py              # Main GUI and simulation interface
├── cfd_visuals.py       # Visualization logic (heatmaps, quivers, streamlines)
├── parameters.txt       # Auto-generated file storing simulation inputs
├── theme.json           # GUI color/style configuration
├── cfd_backend.f90      # Fortran source file (compile manually)
├── cfd.exe              # Compiled binary (do not upload)
├── requirements.txt     # Python package dependencies
├── .gitignore           # Ignore binary/graph files
├── README.md            # This file
```

## Getting Started

### 1. Compile the Fortran Backend

Ensure you have `gfortran` installed. Compile the backend like this:

```bash
gfortran -O2 -o cfd.exe cfd_backend.f90
```


### 2. Install Python Dependencies

Install using:

```bash
pip install -r requirements.txt
```

### 3. Run the Simulator

```bash
python main.py
```

Use the GUI to set domain size, inlet velocity, time steps, solver types, and whether to include a cylinder. After simulation, visualizations will be automatically saved.

### 4. Output Visuals

Generated visuals will be saved to the `graphs/` folder:
- Streamfunction and velocity fields
- Vorticity heatmaps
- Pressure contour plots
- Quiver plots
- Scalar vs time graphs (e.g., Cd, Cl, pressure drop)

## Requirements

Listed in `requirements.txt`:
- `numpy`
- `matplotlib`
- `customtkinter`

## Version

**v0.9.0** — Stable public release. Backend Fortran solver, full GUI, and visualizations complete. Minor improvements in solver modularity and advanced analysis (e.g., streamline density control) are planned for future releases.

## License

MIT License
