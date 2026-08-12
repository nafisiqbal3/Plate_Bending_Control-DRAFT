# Introduction

This repository is associated with the paper **"Feedforward-Feedback Control Methodology for Low-Volume Precision Plate Bending"**, authored by Partho Kundu, Nafis Iqbal, Zongze Li, Balark Tiwari and Robert G. Landers.

The repository contains the experimental data, MATLAB implementation of the process simulation model, and real-time C++ code implementing the image-processing and process-control algorithms presented in the paper.

A real-time, multithreaded C++ program is used to operate the experimental setup, which consists of a servo motor-driven three-axis motion-control system and receives feedback from the servo motion encoders, a camera, and a force sensor.

The study proposes a feedforward-feedback control methodology for efficiently and accurately bending a metal plate. The proposed process model is first used to simulate the bending process, and the corresponding MATLAB implementation is provided in the `src/simulation` directory.

The experimental system is then used to validate the proposed model and control methodology. The real-time C++ image-processing code in the `src/image_processing` directory processes raw camera images to determine the bending angle of the plate. The calculated bending angle is subsequently provided as feedback to the process-control algorithm.

The real-time C++ implementation of the process-control algorithm is provided in the `src/process_controller` directory. The controller uses the measured bending angle and force-sensor feedback to identify the relevant process parameters and calculate the feedforward control action based on the proposed model. A feedback control action is also implemented to compensate for modeling and process uncertainties and achieve the desired bending angle with high precision.
