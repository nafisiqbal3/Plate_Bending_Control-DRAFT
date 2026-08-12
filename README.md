# Feedforward-Feedback Control Methodology for Low-Volume Precision Plate Bending [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository accompanies the paper **"Feedforward-Feedback Control Methodology for Low-Volume Precision Plate Bending"**, authored by Partho Kundu, Nafis Iqbal, Zongze Li, Balark Tiwari, and Robert G. Landers (*Integrating Materials and Manufacturing Innovation*).

## Overview
The repository contains the experimental data, MATLAB implementation of the process simulation model, and real-time C++ code implementing the image-processing and process-control algorithms presented in the paper. The repository is organized as follows:

```text
.
├── data/                   # Dataset for Figures 9-22 in the Paper
├── src/                    
│   ├── simulation/         # Process model and simulation (MATLAB)
│   └── real_time_code/     # Real-time image processing and process control algorithm (C++)
├── CITATION.cff
├── LICENSE                 
└── README.md      
```

The study proposes a feedforward-feedback control methodology for efficiently and accurately bending a metal plate. The proposed process model is first used to simulate the bending process, and the corresponding MATLAB implementation is provided in the `src/simulation` directory. There is also a `README_simulation.md` file that explains how to use that code.

The process model and control strategy proposed in the paper is validated by experiments. The experimental setup used in this study comprises a servo motor-driven three-axis motion-control system, with feedbacks coming from the servo encoders, a visual camera and a force sensor. The system is operated by a multithreaded, real-time C++ application. *This repository does **not** contain the complete code that runs the motion and sensory system in real-time, as that would make the repository hardware specific.* Only the image processing code and the process controlling code are included in `src/real_time_code/` directory. There is an accompanying `README_real_time.md` file that explains how those real-time C++ codes can be used.

The `data/` directory contains the relevant experimental data that have been used to generate the Figure 9-22 in the paper. Each figure has its own subfolder that contains `.csv` files with relevant data from the sensors or simulation results. The accompanying `README.md` files can be instructive.

## Prerequisites
The simulation code only requires MATLAB and accompanying data to be run.

However, the real-time codes are written in C++ and they are only part of a larger multi-threaded real-time program. The code can not be used independently of a closed-loop control system that has a real-time position controlled ram press, a force sensor on the ram, and a visual camera that feeds the control algorithm with real-time data. The program that runs the system should be real-time constrained, better be multi-threaded. We have used a Linux Ubuntu 22.04-based system with `PREEMPT-RT` patch. The additional C++ libraries are:
1. OpenCV: For image processing.
2.


## Citation
If you use this code, dataset, or methodology in your research, please cite our paper:

```text
@article{kundu2026feedforward,
  title={Feedforward-Feedback Control Methodology for Low-Volume Precision Plate Bending},
  author={Kundu, Partho and Iqbal, Nafis and Li, Zongze and Tiwari, Balark and Landers, Robert G.},
  journal={Integrating Materials and Manufacturing Innovation},
  year={2026},
  doi={10.5281/zenodo.XXXXXXX}
}
```
        
