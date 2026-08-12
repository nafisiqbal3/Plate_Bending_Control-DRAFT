# Feedforward-Feedback Control Methodology for Low-Volume Precision Plate Bending [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository accompanies the paper **"Feedforward-Feedback Control Methodology for Low-Volume Precision Plate Bending"**, authored by Partho Kundu, Nafis Iqbal, Zongze Li, Balark Tiwari, and Robert G. Landers (*Integrating Materials and Manufacturing Innovation*).

## Overview
The repository contains the experimental data, MATLAB implementation of the process simulation model, and real-time C++ code implementing the image-processing and process-control algorithms presented in the paper.

The experimental setup used in this study comprises a servo motor-driven three-axis motion-control system, with feedbacks coming from the servo encoders, a visual camera and a force sensor. The system is operated by a multithreaded, real-time C++ application. * This repository does **not** contain the complete code that runs the motion and sensory system in real-time, as that part is hardware specific. Only the image processing code and the process controlling code are part of this repository.*

The study proposes a feedforward-feedback control methodology for efficiently and accurately bending a metal plate. The proposed process model is first used to simulate the bending process, and the corresponding MATLAB implementation is provided in the `src/simulation` directory.

The experimental system is then used to validate the proposed model and control methodology. The real-time C++ image-processing code in the `src/image_processing` directory processes raw camera images to determine the bending angle of the plate. The calculated bending angle is subsequently provided as feedback to the process-control algorithm.

The real-time C++ implementation of the process-control algorithm is provided in the `src/process_controller` directory. The controller uses the measured bending angle and force-sensor feedback to identify the relevant process parameters and calculate the feedforward control action based on the proposed model. A feedback control action is also implemented to compensate for modeling and process uncertainties and achieve the desired bending angle with high precision.

The `data` directory contains the relevant data that have been used to generate the figures in the paper. 

Each directory contains its own *README* file that further explains the usage of the code and data.

## Repository Structure

The repository is organized as follows:

```text
.
├── data/                   # Dataset for Figures 9-22 in the Paper
├── src/                    
│   ├── simulation/         # Process model and simulation (MATLAB)
│   ├── image_processing/   # Real-time image processing algorithm (C++)
│   └── process_controller/ # Real-time process control algorithm (C++)
├── CITATION.cff
├── LICENSE                 
└── README.md              
