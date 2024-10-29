### Adaptive Mesh Refinement in a MacCormack’s Method Navier-Stokes Solver

#### University of California San Diego
#### Course: 290C - Spring 2024
#### Authors: Eden Shokrgozar, Anthony Calderon-Schuler

---

#### Overview
This project adapts MacCormack's method in a Navier-Stokes solver for enhanced computational efficiency through adaptive mesh refinement. It focuses on accurately capturing flow characteristics over a supersonic flat plate, emphasizing areas with intense activity like shock waves and boundary layers.

#### Features
- **Adaptive Mesh Refinement:** Enhances resolution in critical flow regions.
- **MacCormack Method:** Utilized for its robust second-order accuracy.
- **Dynamic Grid Adjustment:** Incorporates changes in the grid structure to optimize computational resources.

#### Usage
- **Main Script:** Run `FinalProject.m` for executing the solver with adaptive mesh refinement.
- **Parameter Adjustment:** Modify simulation parameters through `UserParameters.m` for customized runs.

#### File Structure
- **/functions:** Contains auxiliary functions for the simulation.
- **/Results:** Stores output files, including animations and plots showcasing the refinement process.
- **`FinalProject.m`:** Primary script file.
- **`UserParameters.m`:** Configuration file for setting simulation parameters.

#### How to Run
1. **Set Parameters:** Open `UserParameters.m` and adjust the parameters like Mach number, CFL condition, and visualization preferences.
2. **Execute Simulation:** Load and run `FinalProject.m` to start the simulation. Adjust `iteration` and `refining_interval` within the script for different refinement dynamics.
3. **View Results:** Check the `/Results` directory for output visualizations and data.

