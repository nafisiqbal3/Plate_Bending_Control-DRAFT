#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>


#include "Geometric_S1_and_theta_1_calculation.h"
#include "Displacement_from_Unloaded_Angle.h"

using namespace springback_model;

// Command type returned to the low-level motor hardware
enum class MotionCommandType {
    NONE,           // No action required
    RELATIVE_MOVE,  // Dispatch relative displacement command
    STOP_IMMEDIATE  // Interrupt ongoing motion immediately
};

// Output structure containing hardware handoff commands
struct ActuatorHandoff {
    MotionCommandType command_type = MotionCommandType::NONE;
    double mm_increment_z = 0.0;   // Displacement command for Z-axis (mm)
    double feed_rate = 0.0;        // Target velocity / feed rate (mm/s)
};

// Persistent state structure maintained across control loop iterations
struct ProcessControllerState {
    int state = 1;
    int phase = 0;
    
    // Yield detection tracking
    int state2_sample_count = 0;
    double state2_anchor_z = 0.0;
    double state2_anchor_f = 0.0;
    bool state2_anchor_captured = false;
    double state2_initial_slope = 0.0;
    bool state2_initial_slope_captured = false;
    double yield_force_detected = 0.0;
    double yield_z_detected = 0.0;
    bool yield_point_detected = false;
    
    // Model & Identification outputs
    float theta1 = 0.0;
    float theta2 = 0.0;
    float S1 = 0.0;
    float sigma_0 = 0.0;
    float E_prime = 0.0;
    float F_current = 0.0;
    float Z = 0.0;
    float total_punch = 0.0;
    float scaled_total_punch = 0.0;
    float required_punch = 0.0;
    std::vector<double> Ki_array = std::vector<double>(100, 0.0);
    std::vector<double> x_range = std::vector<double>(100, 0.0);

    // Iterative integral controller tracking
    bool initialized = false;
    float starting_force = 0.0;
    float starting_z = 0.0;
    float prev_command = 0.0;
    float new_command = 0.0;
    float prev_angle = 0.0;
    float total_punch_distance = 0.0;

    ProcessControllerState() {
        for (size_t i = 0; i < x_range.size(); ++i) {
            x_range[i] = 50.0 * static_cast<double>(i) / static_cast<double>(x_range.size() - 1);
        }
    }
};

/**
 * @brief Core feedforward-feedback process controller state machine.
 * 
 * @param[in]  bendAngle             Current real-time plate bending angle (deg) from vision.
 * @param[in]  enc_z                 Current Z-axis encoder displacement (mm).
 * @param[in]  ref_z                 Current Z-axis reference position (mm).
 * @param[in]  F_mag                 Current measured force magnitude (N).
 * @param[in]  command_from_pc_req   Flag indicating low-level motion system is ready for new command.
 * @param[in,out] ctrlState          Persistent state machine variables across iterations.
 * @param[in]  target_angle          Desired target bending angle (deg). Default = 30.0 deg.
 * @return ActuatorHandoff           Structure containing displacement, feed rate, and override commands.
 */
ActuatorHandoff process_controller(double bendAngle,
                                         double enc_z,
                                         double ref_z,
                                         double F_mag,
                                         bool command_from_pc_req,
                                         ProcessControllerState& ctrlState,
                                         double target_angle = 30.0)
{
    ActuatorHandoff handoff;

    // Constants & process parameters
    const float v = 5.0f;
    const float angle_tolerance = 0.1f;
    const float force_tolerance = 1.5f;
    const float t = 1.6f;
    const float w = 25.4f * 2.0f;
    const float distance = 52.0f;
    const float model_correction_factor = 0.95f;
    
    const int state2_initial_skip_N = 10;
    const int state2_slope_step_N = 8;
    const double state2_slope_dev_threshold = 0.015; // 1.5% slope deviation for yield

    switch (ctrlState.state) {

        // --- STATE 1: Initiate Forward Search Move ---
        case 1: {
            if (command_from_pc_req) {
                handoff.command_type = MotionCommandType::RELATIVE_MOVE;
                handoff.mm_increment_z = 80.0;
                handoff.feed_rate = 5.0;
                handoff.override_motion = false;

                ctrlState.state = 2;
            }
            break;
        }

        // --- STATE 2: Yield Detection via Force-Displacement Slope Deviation ---
        case 2: {
            ctrlState.F_current = F_mag;
            ctrlState.Z = enc_z;
            ctrlState.state2_sample_count++;

            if (ctrlState.state2_sample_count < state2_initial_skip_N) {
                break;
            }

            if (!ctrlState.state2_anchor_captured) {
                ctrlState.state2_anchor_z = ctrlState.Z;
                ctrlState.state2_anchor_f = ctrlState.F_current;
                ctrlState.state2_anchor_captured = true;
                break;
            }

            int samples_after_anchor = ctrlState.state2_sample_count - state2_initial_skip_N;
            if (samples_after_anchor % state2_slope_step_N != 0) {
                break;
            }

            double delta_z = ctrlState.Z - ctrlState.state2_anchor_z;
            double delta_f = ctrlState.F_current - ctrlState.state2_anchor_f;

            if (std::abs(delta_z) < 1.0e-6) {
                break;
            }

            double current_slope = delta_f / delta_z;

            if (!ctrlState.state2_initial_slope_captured) {
                if (std::abs(current_slope) > 1.0e-9) {
                    ctrlState.state2_initial_slope = current_slope;
                    ctrlState.state2_initial_slope_captured = true;
                }
                break;
            }

            double slope_deviation = std::abs(current_slope - ctrlState.state2_initial_slope) /
                                      std::abs(ctrlState.state2_initial_slope);

            if (slope_deviation >= state2_slope_dev_threshold) {
                ctrlState.yield_point_detected = true;
                ctrlState.yield_force_detected = ctrlState.F_current;
                ctrlState.yield_z_detected = ctrlState.Z;

                // Stop motion immediately
                handoff.command_type = MotionCommandType::STOP_IMMEDIATE;
                handoff.mm_increment_z = 0.0;
                handoff.feed_rate = 3.0;
                handoff.interrupt_motion = true;
                handoff.override_motion = true;

                ctrlState.state = 3;
            }
            break;
        }

        // --- STATE 3: Compute Yield Stress (sigma_0) and Retract Ram ---
        case 3: {
            if (command_from_pc_req) {
                if (ctrlState.yield_point_detected) {
                    ctrlState.F_current = ctrlState.yield_force_detected;
                    ctrlState.Z = ctrlState.yield_z_detected;
                } else {
                    ctrlState.F_current = F_mag;
                    ctrlState.Z = enc_z;
                }
                ctrlState.theta1 = bendAngle;

                // Analytical calculation of Geometric S1 factor
                auto out = springback_model::Geometric_S1_and_theta_1_calculation(
                    ctrlState.Z, 6.35, 25.4, t, distance, 30.0, 0.0, model_correction_factor
                );
                ctrlState.S1 = out.S1;
                ctrlState.sigma_0 = (6.0f * ctrlState.F_current * ctrlState.S1) / (w * t * t);

                handoff.command_type = MotionCommandType::RELATIVE_MOVE;
                handoff.mm_increment_z = -ref_z;
                handoff.feed_rate = 5.0;
                handoff.override_motion = false;

                ctrlState.state = 4;
            }
            break;
        }

        // --- STATE 4: Compute Effective Modulus (E') and Feedforward Required Stroke ---
        case 4: {
            if (command_from_pc_req) {
                ctrlState.theta2 = bendAngle;
                float spring_back_angle = (ctrlState.theta1 - ctrlState.theta2) * 3.1416f / 180.0f;

                ctrlState.E_prime = (6.0f * ctrlState.F_current * ctrlState.S1 * ctrlState.S1) / 
                                    (spring_back_angle * t * t * t * w);

                auto result = Displacement_from_Unloaded_Angle(
                    ctrlState.E_prime, 0.3, ctrlState.sigma_0, 950.0, 0.1, 6.35, 25.4, 
                    t, w, distance, 30.0, ctrlState.x_range, 1.0, 2.5, 1.5, 
                    model_correction_factor, 0.0, 0.9 * target_angle
                );

                ctrlState.total_punch = result.x_required;
                ctrlState.Ki_array = result.integral_controller_gain;
                ctrlState.scaled_total_punch = 1.0f * ctrlState.total_punch;
                ctrlState.required_punch = ctrlState.scaled_total_punch - enc_z;

                handoff.command_type = MotionCommandType::RELATIVE_MOVE;
                handoff.mm_increment_z = ctrlState.required_punch;
                handoff.feed_rate = v;

                ctrlState.state = 5;
            }
            break;
        }

        // --- STATE 5: Retract after Feedforward Execution ---
        case 5: {
            if (command_from_pc_req) {
                handoff.command_type = MotionCommandType::RELATIVE_MOVE;
                handoff.mm_increment_z = -ctrlState.required_punch;
                handoff.feed_rate = 5.0;

                ctrlState.state = 6;
            }
            break;
        }

        // --- STATE 6: Detect Tool Detachment ---
        case 6: {
            ctrlState.F_current = F_mag;
            if (ctrlState.F_current < force_tolerance) {
                handoff.command_type = MotionCommandType::STOP_IMMEDIATE;
                handoff.mm_increment_z = 0.0;
                handoff.feed_rate = 3.0;
                handoff.interrupt_motion = true;
                handoff.override_motion = true;

                ctrlState.state = 7;
            }
            break;
        }

        // --- STATE 7: Iterative Integral Closed-Loop Feedback Control ---
        case 7: {
            double current_angle = bendAngle;
            double angle_error = target_angle - current_angle;

            if (ctrlState.phase == 0 && command_from_pc_req) {
                if (angle_error < angle_tolerance) {
                    ctrlState.state = 8;
                } else {
                    float detachment_distance = enc_z;
                    handoff.override_motion = false;

                    double gain = 0.75;
                    auto it = std::lower_bound(ctrlState.x_range.begin(), ctrlState.x_range.end(), ctrlState.total_punch_distance);
                    size_t idx = std::distance(ctrlState.x_range.begin(), it);
                    gain = (idx == 0) ? ctrlState.Ki_array.front() :
                           (idx == ctrlState.x_range.size()) ? ctrlState.Ki_array.back() :
                           ctrlState.Ki_array[idx - 1] + (ctrlState.Ki_array[idx] - ctrlState.Ki_array[idx - 1]) *
                           (ctrlState.total_punch_distance - ctrlState.x_range[idx - 1]) /
                           (ctrlState.x_range[idx] - ctrlState.x_range[idx - 1]);

                    if (!ctrlState.initialized) {
                        ctrlState.starting_force = F_mag;
                        ctrlState.starting_z = enc_z;
                        ctrlState.prev_command = ctrlState.scaled_total_punch - enc_z;

                        auto it_init = std::lower_bound(ctrlState.x_range.begin(), ctrlState.x_range.end(), ctrlState.scaled_total_punch);
                        size_t idx_init = std::distance(ctrlState.x_range.begin(), it_init);
                        gain = (idx_init == 0) ? ctrlState.Ki_array.front() :
                               (idx_init == ctrlState.x_range.size()) ? ctrlState.Ki_array.back() :
                               ctrlState.Ki_array[idx_init - 1] + (ctrlState.Ki_array[idx_init] - ctrlState.Ki_array[idx_init - 1]) *
                               (ctrlState.scaled_total_punch - ctrlState.x_range[idx_init - 1]) /
                               (ctrlState.x_range[idx_init] - ctrlState.x_range[idx_init - 1]);

                        ctrlState.initialized = true;
                    }

                    ctrlState.new_command = ctrlState.prev_command + gain * angle_error;
                    ctrlState.prev_command = ctrlState.new_command;
                    ctrlState.new_command -= (detachment_distance - ctrlState.starting_z);
                    ctrlState.prev_angle = current_angle;

                    ctrlState.total_punch_distance = ctrlState.new_command + enc_z;

                    handoff.command_type = MotionCommandType::RELATIVE_MOVE;
                    handoff.mm_increment_z = ctrlState.new_command;
                    handoff.feed_rate = v;

                    ctrlState.phase = 1;
                }
            } else if (ctrlState.phase == 1 && command_from_pc_req) {
                handoff.command_type = MotionCommandType::RELATIVE_MOVE;
                handoff.mm_increment_z = -ctrlState.new_command;
                handoff.feed_rate = 5.0;

                ctrlState.phase = 2;
            } else if (ctrlState.phase == 2) {
                ctrlState.F_current = F_mag;
                if (ctrlState.F_current < force_tolerance) {
                    handoff.command_type = MotionCommandType::STOP_IMMEDIATE;
                    handoff.mm_increment_z = 0.0;
                    handoff.feed_rate = 5.0;
                    handoff.interrupt_motion = true;
                    handoff.override_motion = true;

                    ctrlState.phase = 0;
                }
            }
            break;
        }

        // --- STATE 8: Process Complete - Retreat Ram to Home Position ---
        case 8: {
            if (command_from_pc_req) {
                handoff.command_type = MotionCommandType::RELATIVE_MOVE;
                handoff.mm_increment_z = -enc_z;
                handoff.feed_rate = 15.0;

                ctrlState.state = 9;
            }
            break;
        }

        // --- STATE 9: Final State ---
        case 9: {
            if (command_from_pc_req) {
                ctrlState.state = 10; // Complete
            }
            break;
        }

        default:
            break;
    }

    return handoff;
}
