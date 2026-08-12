#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace springback_model {

struct GeometricS1Theta1Result {
    double S1 = std::numeric_limits<double>::quiet_NaN();
    double theta_1 = std::numeric_limits<double>::quiet_NaN();
    double Wrapped_arc_length = std::numeric_limits<double>::quiet_NaN();

    int phase_S1 = 0;       // 1 = phase 1, 2 = phase 2
    int phase_theta_1 = 0;  // 1 = phase 1, 2 = phase 2
};

namespace geom_detail {

constexpr double PI = 3.14159265358979323846;

inline double deg2rad(double deg)
{
    return deg * PI / 180.0;
}

inline double rad2deg(double rad)
{
    return rad * 180.0 / PI;
}

inline double clamp_val(double v, double lo, double hi)
{
    return std::max(lo, std::min(v, hi));
}

inline double local_theta1_phase1_geom(double xx, double S, double Delta_y)
{
    const double D = std::sqrt((S - xx) * (S - xx) + Delta_y * Delta_y);
    const double delta = std::atan2(Delta_y, std::max(S - xx, 1e-12));
    const double lambda = std::acos(clamp_val(S / D, -1.0, 1.0));

    return rad2deg(delta - lambda);
}

inline double solve_theta_v_phase2_geom(
    double xx,
    double S,
    double Delta_y,
    double Rt,
    double t,
    double alpha_rad)
{
    auto residual = [&](double th) -> double
    {
        return th - std::atan2(
            Delta_y - Rt * std::sin(alpha_rad) - (t / 2.0) * std::cos(th),
            -S + xx + Rt * std::cos(alpha_rad) + (t / 2.0) * std::sin(th)
        );
    };

    const int n_guess = 200;
    const double lo = 1e-4;
    const double hi = PI / 2.0 - 1e-4;

    double theta_v = std::numeric_limits<double>::quiet_NaN();
    double best_res = std::numeric_limits<double>::infinity();

    for (int i = 0; i < n_guess; ++i)
    {
        const double th = lo + (hi - lo) * static_cast<double>(i) /
                                   static_cast<double>(n_guess - 1);

        const double res = std::abs(residual(th));

        if (std::isfinite(res) && res < best_res)
        {
            best_res = res;
            theta_v = th;
        }
    }

    if (!std::isfinite(theta_v))
    {
        theta_v = 0.0;
    }

    return theta_v;
}

inline double local_theta1_phase2_geom(
    double xx,
    double S,
    double Delta_y,
    double Rt,
    double Rd,
    double t,
    double alpha_rad)
{
    const double theta_v = solve_theta_v_phase2_geom(
        xx, S, Delta_y, Rt, t, alpha_rad
    );

    double denom_sin = std::sin(theta_v);

    if (std::abs(denom_sin) < 1e-9)
    {
        denom_sin = (denom_sin >= 0.0 ? 1.0 : -1.0) * 1e-9;
    }

    const double Rn =
        (Delta_y - Rt * std::sin(alpha_rad) - (t / 2.0) * std::cos(theta_v))
        / denom_sin;

    const double ratio_asin = clamp_val((Rd + t / 2.0) / Rn, -1.0, 1.0);
    const double theta_w = std::asin(ratio_asin);

    return rad2deg((PI / 2.0) - theta_v + theta_w);
}

} // namespace geom_detail


inline GeometricS1Theta1Result Geometric_S1_and_theta_1_calculation(
    double x,
    double Rd,
    double Rt,
    double t,
    double d,
    double alpha_deg,
    double extended_length,
    double model_correction_factor)
{
    GeometricS1Theta1Result result;

    const double S = Rd + Rt + t;
    const double alpha_rad = geom_detail::deg2rad(alpha_deg);

    // Use original d for S1.
    const double d_S1 = d;

    // Use corrected d for theta_1.
    const double d_theta = d * model_correction_factor;

    const double Delta_y_S1 = Rd + d_S1;
    const double Delta_y_theta = Rd + d_theta;

    // -----------------------------
    // S1 calculation
    // -----------------------------
    const double theta1_phase1_S1 =
        geom_detail::local_theta1_phase1_geom(x, S, Delta_y_S1);

    if (theta1_phase1_S1 <= alpha_deg)
    {
        result.phase_S1 = 1;

        const double D =
            std::sqrt((S - x) * (S - x) + Delta_y_S1 * Delta_y_S1);

        const double delta =
            std::atan2(Delta_y_S1, std::max(S - x, 1e-12));

        const double lambda =
            std::acos(geom_detail::clamp_val(S / D, -1.0, 1.0));

        const double straight_len =
            std::sqrt(std::max(D * D - S * S, 0.0));

        const double wrap_len =
            (Rd + t / 2.0) * (delta - lambda);

        result.S1 = straight_len + wrap_len + extended_length;
        result.Wrapped_arc_length = wrap_len;
    }
    else
    {
        result.phase_S1 = 2;

        const double theta_v_S1 =
            geom_detail::solve_theta_v_phase2_geom(
                x, S, Delta_y_S1, Rt, t, alpha_rad
            );

        double denom_sin = std::sin(theta_v_S1);

        if (std::abs(denom_sin) < 1e-9)
        {
            denom_sin = (denom_sin >= 0.0 ? 1.0 : -1.0) * 1e-9;
        }

        const double Rn =
            (Delta_y_S1 - Rt * std::sin(alpha_rad)
             - (t / 2.0) * std::cos(theta_v_S1))
            / denom_sin;

        const double straight_len =
            std::sqrt(
                std::max(
                    Rn * Rn - (Rd + t / 2.0) * (Rd + t / 2.0),
                    0.0
                )
            );

        const double theta_w =
            std::asin(
                geom_detail::clamp_val((Rd + t / 2.0) / Rn, -1.0, 1.0)
            );

        const double theta1_tmp =
            (geom_detail::PI / 2.0) - theta_v_S1 + theta_w;

        const double wrap_len =
            (Rd + t / 2.0) * theta1_tmp;

        result.S1 = straight_len + wrap_len + extended_length;
        result.Wrapped_arc_length = wrap_len;
    }

    // -----------------------------
    // theta_1 calculation
    // -----------------------------
    const double theta1_phase1_theta =
        geom_detail::local_theta1_phase1_geom(x, S, Delta_y_theta);

    if (theta1_phase1_theta <= alpha_deg)
    {
        result.phase_theta_1 = 1;
        result.theta_1 = theta1_phase1_theta;
    }
    else
    {
        result.phase_theta_1 = 2;
        result.theta_1 =
            geom_detail::local_theta1_phase2_geom(
                x, S, Delta_y_theta, Rt, Rd, t, alpha_rad
            );
    }

    return result;
}

} // namespace springback_model