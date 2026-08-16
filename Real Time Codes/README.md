# Real Time Codes

Real-time image-processing and process-control C++ code for the experimental plate-bending/unbending system. These are the vision and control-logic building blocks that plug into a larger multithreaded, real-time program driving the motion-control system, force sensor, and camera — this code is **not** a standalone, runnable application.

## Contents

| File | Description |
|---|---|
| `image_processing_func.cpp` | `image_processor()` — takes a raw BGR camera frame, thresholds a white-plate ROI to find its edge, fits a line to the detected points with OpenCV, and returns the plate's bending angle (relative to a baseline captured on loop 15). |
| `process_controller_func.cpp` | `process_controller()` — the feedforward-feedback state machine that drives the process. Sequences the ram through search, yield detection, spring-back/parameter identification, feedforward stroke calculation, and closed-loop angle correction, issuing motion commands (`ActuatorHandoff`) at each step. |
| `Geometric_S1_and_theta_1_calculation.h` | `springback_model::Geometric_S1_and_theta_1_calculation()` — geometry model for the **bending** process: computes the wrap length/arc-length factor `S1` and loaded bend angle `theta_1` from ram displacement and tooling geometry, handling both the pre-contact and post-contact (die-wrap) phases. |
| `Displacement_from_Unloaded_Angle.h` | `springback_model::Displacement_from_Unloaded_Angle()` — forward process model for **bending**, relating ram displacement to force, moment, curvature, and spring-back across the elastic/plastic/die-conforming regions. Used to build the `x`-vs-`theta_1`/spring-back curves, the integral-controller gain schedule, and to solve for the ram displacement (`x_required`) needed to hit a desired unloaded bend angle. |
| `Geometric_S1_and_theta_1_calculation_unbending.h` | `springback_model::Geometric_S1_and_theta_1_calculation_unbending()` — mirror-image geometry model for the **unbending** process: starts from a plate pre-bent at a negative initial angle and computes `S1`/`theta_1` as it's straightened, including a reverse-phase solve to map the initial offset angle back to an equivalent ram displacement. |
| `Displacement_from_Unloaded_Angle_Unbending.h` | `springback_model::Displacement_from_Unloaded_Angle_Unbending()` — forward process model for **unbending**, the counterpart to the bending displacement model above, parameterized by an `initial_angle_offset` and used the same way to compute the feedforward stroke and gain schedule for a target unbent angle. |

## Bending vs. unbending

There are two parallel sets of process-model headers:

- **Bending** — `Geometric_S1_and_theta_1_calculation.h` + `Displacement_from_Unloaded_Angle.h` — for a plate starting flat (or at a positive angle) and being bent further.
- **Unbending** — `Geometric_S1_and_theta_1_calculation_unbending.h` + `Displacement_from_Unloaded_Angle_Unbending.h` — for a plate starting pre-bent at a negative angle and being straightened back out.

`process_controller_func.cpp` includes both sets — use whichever one matches your experiment; they aren't meant to run simultaneously.

## How the pieces fit together

1. **Vision** (`image_processing_func.cpp`) continuously estimates the current plate bend angle from the camera feed.
2. **Process controller** (`process_controller_func.cpp`) consumes that angle plus encoder position and force-sensor readings, and steps through its state machine:
   - **States 1–2** — drive the ram forward and detect the material yield point from a deviation in the force–displacement slope.
   - **State 3** — use the geometric header to back out the geometric factor `S1` and estimate yield stress (`sigma_0`) at the detected yield point.
   - **State 4** — retract, measure spring-back, estimate effective modulus (`E'`), and call the displacement header to compute the feedforward ram stroke and gain schedule needed to reach the target angle.
   - **States 5–6** — execute the feedforward stroke and detect tool detachment.
   - **State 7** — run iterative closed-loop (integral) correction using the vision-measured angle until the target angle is reached within tolerance.
   - **States 8–9** — retract the ram to home and complete the cycle.

## Prerequisites

- C++20 (or later) compiler
- [OpenCV](https://opencv.org/) (for `image_processing_func.cpp`)
- No external dependencies for the four `springback_model` headers beyond the standard library

This code assumes it is wired into a real-time control loop with live encoder, force, and camera feedback — constants marked `//Change accordingly` in `process_controller_func.cpp` (tolerances, tool geometry, timing) should be set to match your specific press, tooling, and material.

## Reference

See the top-level [README](../README.md) and the associated paper for the full experimental setup and derivation of the process models.
