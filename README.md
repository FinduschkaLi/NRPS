# The Network Reduced Power System (NRPS) and derived models

This repository contain files for running the NRPS and related models in MATLAB, for generating the admittance matrix, performing the network reduction and for visualizing the graph of the network.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Files](#files)
- [References](#references)

## Features

The Network Reduced Power System (NRPS) Model is a powerful tool for modeling electric grids of (virtual) synchronous machines. This repository contains files for running the Nonlinear Kuramoto, the NRPS model, the Friction Enhanced Power System (FEPS) model, the Reactive Active Power Power System (RAPS) model.

## Getting started

To use the NRPS Model MATLAB code, follow these steps:

1. Clone this repository: `git clone https://github.com/FinduschkaLi/NRPS.git`
2. Open MATLAB and navigate to the project directory.
3. Run the file `example_NRPS.m` to see an example of how to use the models.

## Files

### Parameterization

- **reduceNetwork.m:** File used for generating the reduced admittance matrix YRED. Input is an array of tie-line impedances and an array of loads aswell as the number of generators.
- **generateNRPSParameters.m:** File used for generating the matrices A, PHI and the vector P for the NRPS model. Input is the admittance matrix YRED.
- **generateRAPSParameters.m:** File used to calculate the admittance matrix between the output / filters of the VSMs and the matrix to obtain v from e.

### Power calculations

- **getP.m:** File used to calculate the output power of the VSMs / generators
- **getVPQ.m:** File used to calculate the voltage and powers of the VSMs / generators for RAPS and NRPS models

### Execution

- **RungeK.m:** File used for solving a differential equation for a time vector t

### Linearization

- **linearizeFEPS.m:** This file calculates a linear model for the FEPS model

### Model files

- **NK.m:** Nonlinear Kuramoto model
- **NRPS.m:** Regular Network Reduced Power System model.
- **FEPS.m:** Friction enhanced power system model (NRPS with Virtual Friction)
- **RAPS.m:** Reactive Active Power Power System model (third order model with voltage dynamics, an extended version of the NRPS model)

### Example files

- **example_NRPS.m:** Example file for running the NRPS model and linearizing it
- **example_RAPS.m:** Example file for running the RAPS model

### Visualization

- **drawNetwork.m:** File used for drawing the network graph

## References

The following papers give an explanation of the models

- NRPS model: [1] Weiss, George, Florian Dörfler, and Yoash Levron. "A stability theorem for networks containing synchronous generators." Systems & Control Letters 134 (2019): 104561.
- FEPS model: [2] Reissner, Florian, Hang Yin, and George Weiss. "A stability result for network reduced power systems using virtual friction and inertia." IEEE Transactions on Smart Grid 13.3 (2021): 1668-1678.
- RAPS model: [3] Reissner, Florian, and George Weiss. "The NRPS model and its extensions for modeling grids with multiple
  (virtual) synchronous machines."8th IEEE Workshop on the Electronic Grid (EGRID). IEEE, 2023 (accepted).
