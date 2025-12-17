# Numerical-FEM Project

## Overview
This project is a 2D FEM solver for the heat conduction equation, designed to solve numerical problems using various computational methods. The project includes implementations for the stationary, instationary and non linear heat equation with the corresponding numerical solvers, Lagrange interpolation, and more. It is structured to provide modular and reusable components for numerical analysis and visualization.

This project specifically includes a 2D FEM solver for the heat conduction equation.

## Features
- **Finite Element Methods (FEM):**
  - Stationary and instationary FEM solvers.
  - Mesh generation and visualization.
- **Numerical Solvers:**
  - Gaussian elimination.
  - Gradient descent and conjugate gradient methods.
- **Interpolation:**
  - Lagrange interpolation in 1D and 2D.
- **Integration:**
  - Gaussian quadrature.
  - Trapezoidal and midpoint rules.
- **Ordinary Differential Equations (ODEs):**
  - Forward Euler, Backward Euler, and more.
- **Visualization:**
  - Mesh plotting and solution visualization.

## File Structure
- `FEMstationar.py`: Stationary FEM solver.
- `FEMinstationar.py`: Instationary FEM solver.
- `FEMstatmeshgen.py`: Mesh generation for FEM.
- `integration.py`: Numerical integration methods.
- `lagrange.py`: Lagrange interpolation in 1D.
- `lagrange2D.py`: Lagrange interpolation in 2D.
- `linear_solver.py`: Linear solvers.
- `nlinear.py`: Nonlinear FEM solver.
- `ode.py`: Solvers for ordinary differential equations.
- `triagplot.py`: Visualization of triangular and quadrilateral meshes.

## How to Run
1. Clone the repository or download the project files.
2. Ensure you have Python installed (version 3.8 or higher).
3. Install the required libraries:
   ```bash
   pip install numpy matplotlib scipy
   ```
4. Run the desired script. For example, to run the stationary FEM solver:
   ```bash
   python FEMstationar.py
   ```

## Dependencies
- Python 3.8+
- NumPy
- Matplotlib
- SciPy

## Examples
### Running the Stationary FEM Solver
```bash
python FEMstationar.py
```
This will compute the solution for a stationary FEM problem and visualize the results.

### Visualizing a Mesh
```bash
python FEMstatmeshgen.py
```
This will use the mesh generation script to solve the stationary FEM problem.

### Modifying Geometry and Constants
The geometry used in this project is a plate with one corner cut off in the form of a circle. You can modify the dimensions of the plate and the radius of the circular cut-off directly in the code. For example:
- **Plate Dimensions:** Adjust the width and height variables in the correlating scripts (i.e. `FEMstatmeshgen.py`).
- **Circular Cut-Off:** Change the radius parameter to modify the size of the circular cut-off.
- **Physical Constants:** Update values like `thermal_conductivity` in the FEM solver scripts to reflect different material properties.
-**Mesh Generation and Geometry Changes:** Change the Mesh used for the FEM solver or vary the geometry entirely.



