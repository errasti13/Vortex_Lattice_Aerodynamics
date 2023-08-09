# Vortex-Lattice Method for Incompressible Aerodynamics Analysis

## Overview

This repository contains MATLAB code that implements the Vortex-Lattice method for analyzing incompressible steady-state aerodynamics of a given wing geometry. The code calculates and visualizes aerodynamic coefficients, lift and drag distributions, and other relevant parameters for a specified range of angles of attack (alpha) or sweep angles (psi).

The code is part of a collaborative project undertaken by Grupo 13: Alonso Lucas, Sara; Errasti Odriozola, Jon; Sarabia Vargas, Alejandro; and Terreros Sanchez, Carlos. Serves as an educational tool to understand and demonstrate the application of the Vortex-Lattice method.

## Features

- Calculation of lift and drag coefficients for various angles of attack or sweep angles.
- Visualization of lift and drag distributions along the wing span.
- Customizable parameters for wing geometry, flow conditions, and analysis range.
- Utilizes a custom function for vortex panel integration (LeyBiotSavart3D).

## Usage

1. Clone this repository to your local machine.
2. Open the MATLAB script `MainScript.m`.
3. Configure the parameters in the script's setup section for your analysis.
4. Run the script to perform aerodynamic calculations and generate results.
5. View and interpret the generated plots to analyze aerodynamic behavior.

## Dependencies

- MATLAB (R201X or later)
- LeyBiotSavart3D function (provided in the repository)

## License

This project is licensed under the [MIT License](LICENSE).

## Contributing

Contributions to the project are welcome. If you find any issues or have suggestions for improvements, please feel free to open an issue or submit a pull request.
