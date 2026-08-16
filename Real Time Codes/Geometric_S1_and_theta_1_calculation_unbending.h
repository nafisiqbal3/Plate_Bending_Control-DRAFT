#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace springback_model {

struct GeometricS1Theta1Result {
    double S1 = std::numeric_limits<double>::quiet_NaN();
    double theta_1 = std::numeric_limits<double>::quiet_NaN();
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

inline double safe_sin_denom(double v)
{
    if (std::abs(v) < 1e-9) {
        return (v >= 0.0 ? 1.0 : -1.0) * 1e-9;
    }
    return v;
}

inline double reverse_phase1_x_from_theta(
    double theta_rad_positive,
    double Rd,
    double Rt,
    double t,
    double d)
{
    const double Rt_n = Rt + t / 2.0;

    return (Rt_n * std::sin(theta_rad_positive) + d + Rd) * std::tan(theta_rad_positive)
           + Rt_n * std::cos(theta_rad_positive)
           - Rt_n;
}

inline double reverse_phase1_S1_from_theta(
    double theta_rad_positive,
    double Rd,
    double Rt,
    double t,
    double d)
{
    const double Rt_n = Rt + t / 2.0;
    const double denom = std::cos(theta_rad_positive);

    if (std::abs(denom) < 1e-12) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return (Rt_n * std::sin(theta_rad_positive) + d + Rd) / denom;
}

inline double solve_reverse_theta_from_x(
    double x_target_positive,
    double theta_offset_positive_rad,
    double Rd,
    double Rt,
    double t,
    double d)
{
    auto fun = [&](double th) -> double
    {
        return reverse_phase1_x_from_theta(th, Rd, Rt, t, d) - x_target_positive;
    };

    double th_lo = 0.0;
    double th_hi = theta_offset_positive_rad;

    if (th_lo > th_hi) {
        std::swap(th_lo, th_hi);
    }

    const double f_lo = fun(th_lo);
    const double f_hi = fun(th_hi);

    if (std::isfinite(f_lo) && std::isfinite(f_hi) && f_lo * f_hi <= 0.0) {
        double a = th_lo;
        double b = th_hi;
        double fa = f_lo;

        for (int iter = 0; iter < 100; ++iter) {
            const double mid = 0.5 * (a + b);
            const double fm = fun(mid);

            if (!std::isfinite(fm)) {
                break;
            }

            if (std::abs(fm) < 1e-12 || std::abs(b - a) < 1e-12) {
                return mid;
            }

            if (fa * fm <= 0.0) {
                b = mid;
            } else {
                a = mid;
                fa = fm;
            }
        }

        return 0.5 * (a + b);
    }

    const int n_scan = 1000;
    double best_theta = th_lo;
    double best_res = std::numeric_limits<double>::infinity();

    for (int i = 0; i < n_scan; ++i) {
        const double th = th_lo + (th_hi - th_lo) * static_cast<double>(i) /
                                      static_cast<double>(n_scan - 1);

        const double res = std::abs(fun(th));

        if (std::isfinite(res) && res < best_res) {
            best_res = res;
            best_theta = th;
        }
    }

    return best_theta;
}

inline double local_theta1_phase1(
    double xx,
    double S,
    double Delta_y)
{
    const double D = std::sqrt((S - xx) * (S - xx) + Delta_y * Delta_y);
    const double delta = std::atan2(Delta_y, std::max(S - xx, 1e-12));
    const double lambda = std::acos(clamp_val(S / D, -1.0, 1.0));

    return rad2deg(delta - lambda);
}

inline double local_S1_phase1(
    double xx,
    double S,
    double Delta_y,
    double Rd,
    double t)
{
    const double D = std::sqrt((S - xx) * (S - xx) + Delta_y * Delta_y);
    const double delta = std::atan2(Delta_y, std::max(S - xx, 1e-12));
    const double lambda = std::acos(clamp_val(S / D, -1.0, 1.0));

    const double straight_len = std::sqrt(std::max(D * D - S * S, 0.0));
    const double wrap_len = (Rd + t / 2.0) * (delta - lambda);

    return straight_len + wrap_len;
}

inline double solve_theta_v_phase2(
    double xx,
    double S,
    double Delta_y,
    double Rt,
    double t,
    double alpha_rad)
{
    auto fun = [&](double th) -> double
    {
        return th - std::atan2(
            Delta_y - Rt * std::sin(alpha_rad) - (t / 2.0) * std::cos(th),
            -S + xx + Rt * std::cos(alpha_rad) + (t / 2.0) * std::sin(th)
        );
    };

    const int n_guess = 200;
    const double lo = 1e-6;
    const double hi = PI - 1e-6;

    double theta_v = std::numeric_limits<double>::quiet_NaN();
    double best_res = std::numeric_limits<double>::infinity();

    for (int i = 0; i < n_guess; ++i) {
        const double th = lo + (hi - lo) * static_cast<double>(i) /
                                   static_cast<double>(n_guess - 1);

        const double res = std::abs(fun(th));

        if (std::isfinite(res) && res < best_res) {
            best_res = res;
            theta_v = th;
        }
    }

    if (!std::isfinite(theta_v)) {
        theta_v = PI / 4.0;
    }

    return clamp_val(theta_v, 1e-6, PI - 1e-6);
}

inline double local_theta1_phase2(
    double xx,
    double S,
    double Delta_y,
    double Rt,
    double Rd,
    double t,
    double alpha_rad)
{
    const double theta_v = solve_theta_v_phase2(xx, S, Delta_y, Rt, t, alpha_rad);
    const double denom_sin = safe_sin_denom(std::sin(theta_v));

    const double L_center =
        (Delta_y - Rt * std::sin(alpha_rad) - (t / 2.0) * std::cos(theta_v))
        / denom_sin;

    const double ratio_theta_w =
        clamp_val((Rd + t / 2.0) / L_center, -1.0, 1.0);

    const double theta_w = std::asin(ratio_theta_w);

    return rad2deg((PI / 2.0) - theta_v + theta_w);
}

inline double local_S1_phase2(
    double xx,
    double S,
    double Delta_y,
    double Rt,
    double Rd,
    double t,
    double alpha_rad)
{
    const double theta_v = solve_theta_v_phase2(xx, S, Delta_y, Rt, t, alpha_rad);
    const double denom_sin = safe_sin_denom(std::sin(theta_v));

    const double L_center =
        (Delta_y - Rt * std::sin(alpha_rad) - (t / 2.0) * std::cos(theta_v))
        / denom_sin;

    const double straight_len =
        std::sqrt(std::max(L_center * L_center - (Rd + t / 2.0) * (Rd + t / 2.0), 0.0));

    const double ratio_theta_w =
        clamp_val((Rd + t / 2.0) / L_center, -1.0, 1.0);

    const double theta_w = std::asin(ratio_theta_w);

    const double theta_tmp = (PI / 2.0) - theta_v + theta_w;
    const double wrap_len = (Rd + t / 2.0) * theta_tmp;

    return straight_len + wrap_len;
}

inline double find_phase23_x_alpha(
    double S,
    double Delta_y_theta,
    double alpha_deg)
{
    const double theta_raw_0 = local_theta1_phase1(0.0, S, Delta_y_theta);

    auto fun = [&](double xx) -> double
    {
        return local_theta1_phase1(xx, S, Delta_y_theta) - theta_raw_0 - alpha_deg;
    };

    double x_lo = 0.0;
    double x_hi = 1.0;

    while (std::isfinite(fun(x_hi)) && fun(x_hi) < 0.0 && x_hi < 1000.0) {
        x_hi = x_hi * 1.5 + 1.0;
    }

    const double f_lo = fun(x_lo);
    const double f_hi = fun(x_hi);

    if (std::isfinite(f_lo) && std::isfinite(f_hi) && f_lo * f_hi <= 0.0) {
        double a = x_lo;
        double b = x_hi;
        double fa = f_lo;

        for (int iter = 0; iter < 100; ++iter) {
            const double mid = 0.5 * (a + b);
            const double fm = fun(mid);

            if (!std::isfinite(fm)) {
                break;
            }

            if (std::abs(fm) < 1e-10 || std::abs(b - a) < 1e-10) {
                return mid;
            }

            if (fa * fm <= 0.0) {
                b = mid;
            } else {
                a = mid;
                fa = fm;
            }
        }

        return 0.5 * (a + b);
    }

    const int n_scan = 2000;
    double best_x = x_lo;
    double best_res = std::numeric_limits<double>::infinity();

    for (int i = 0; i < n_scan; ++i) {
        const double xx = x_lo + (x_hi - x_lo) * static_cast<double>(i) /
                                      static_cast<double>(n_scan - 1);

        const double res = std::abs(fun(xx));

        if (std::isfinite(res) && res < best_res) {
            best_res = res;
            best_x = xx;
        }
    }

    return best_x;
}

}


inline GeometricS1Theta1Result Geometric_S1_and_theta_1_calculation_unbending(
    double x,
    double Rd,
    double Rt,
    double t,
    double d,
    double alpha_deg,
    double extended_length,
    double curvature_factor,
    double theta_initial_deg)
{
    GeometricS1Theta1Result result;

    if (x < 0.0) {
        return result;
    }

    if (theta_initial_deg >= 0.0) {
        return result;
    }

    const double S = Rd + Rt + t;
    const double Rn = Rd + t / 2.0;

    const double alpha_rad = geom_detail::deg2rad(alpha_deg);

    const double d_theta = d * curvature_factor;
    const double d_S1 = d;

    const double Delta_y_theta = Rd + d_theta;
    const double Delta_y_S1 = Rd + d_S1;

    const double theta_offset_abs_deg = std::abs(theta_initial_deg);
    const double theta_offset_abs_rad = geom_detail::deg2rad(theta_offset_abs_deg);

    const double x_offset_positive =
        geom_detail::reverse_phase1_x_from_theta(theta_offset_abs_rad, Rd, Rt, t, d_theta);

    const double x_offset_geom = -x_offset_positive;

    const double x_geom = x + x_offset_geom;

    const double theta_shift = theta_offset_abs_deg;

    const double S1_reverse_end =
        geom_detail::reverse_phase1_S1_from_theta(0.0, Rd, Rt, t, d_S1);

    const double theta_raw_phase2_0 =
        geom_detail::local_theta1_phase1(0.0, S, Delta_y_theta);

    const double S1_raw_phase2_0 =
        geom_detail::local_S1_phase1(0.0, S, Delta_y_S1, Rd, t);

    const double S1_phase2_offset =
        S1_reverse_end - S1_raw_phase2_0;

    const double x_alpha_geom =
        geom_detail::find_phase23_x_alpha(S, Delta_y_theta, alpha_deg);

    const double theta_raw_phase3_alpha =
        geom_detail::local_theta1_phase2(
            x_alpha_geom, S, Delta_y_theta, Rt, Rd, t, alpha_rad
        );

    const double theta_phase3_offset =
        alpha_deg - theta_raw_phase3_alpha;

    const double S1_phase2_alpha_raw =
        geom_detail::local_S1_phase1(
            x_alpha_geom, S, Delta_y_S1, Rd, t
        );

    const double S1_phase2_alpha =
        S1_phase2_alpha_raw + S1_phase2_offset;

    const double S1_phase3_alpha_raw =
        geom_detail::local_S1_phase2(
            x_alpha_geom, S, Delta_y_S1, Rt, Rd, t, alpha_rad
        );

    const double S1_phase3_offset =
        S1_phase2_alpha - S1_phase3_alpha_raw;

    if (x_geom <= 0.0) {
        const double x_phase1_positive = -x_geom;

        const double theta_geom_rad =
            geom_detail::solve_reverse_theta_from_x(
                x_phase1_positive,
                theta_offset_abs_rad,
                Rd,
                Rt,
                t,
                d_theta
            );

        const double theta_geom_deg =
            geom_detail::rad2deg(theta_geom_rad);

        result.theta_1 =
            theta_offset_abs_deg - theta_geom_deg;

        result.S1 =
            geom_detail::reverse_phase1_S1_from_theta(
                theta_geom_rad,
                Rd,
                Rt,
                t,
                d_S1
            )
            + extended_length;

        return result;
    }

    if (x_geom <= x_alpha_geom) {
        const double theta_raw_phase2_x =
            geom_detail::local_theta1_phase1(x_geom, S, Delta_y_theta);

        const double theta_phase2_data =
            theta_raw_phase2_x - theta_raw_phase2_0;

        result.theta_1 =
            theta_phase2_data + theta_shift;

        const double S1_phase2_data =
            geom_detail::local_S1_phase1(x_geom, S, Delta_y_S1, Rd, t);

        result.S1 =
            S1_phase2_data + S1_phase2_offset + extended_length;

        return result;
    }

    const double theta_phase3_data =
        geom_detail::local_theta1_phase2(
            x_geom,
            S,
            Delta_y_theta,
            Rt,
            Rd,
            t,
            alpha_rad
        );

    result.theta_1 =
        theta_phase3_data + theta_phase3_offset + theta_shift;

    const double S1_phase3_data =
        geom_detail::local_S1_phase2(
            x_geom,
            S,
            Delta_y_S1,
            Rt,
            Rd,
            t,
            alpha_rad
        );

    result.S1 =
        S1_phase3_data + S1_phase3_offset + extended_length;

    return result;
}

}