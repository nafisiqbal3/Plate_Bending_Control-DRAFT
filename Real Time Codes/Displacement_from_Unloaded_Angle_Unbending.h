#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace springback_model {

struct UnbendingModelData {
    std::vector<double> x;
    std::vector<double> theta_1;
    std::vector<double> theta_2;
    std::vector<double> Springback;
    std::vector<double> Clamped_Moment;
    std::vector<double> Clamped_Curvature;
    std::vector<double> Force;
    std::vector<double> Trial_Force;
    std::vector<double> S1;
    std::vector<double> Sy;
    std::vector<double> Wrapped_arc_length;
    std::vector<double> delta_elastic_local;
    std::vector<double> delta_R2_local;
    std::vector<double> delta_plastic_local;
    std::vector<double> Ki_array;
    std::vector<int> phase_flag;
};

struct DisplacementUnbendingResult {
    double x_required = std::numeric_limits<double>::quiet_NaN();
    double theta_1_at_x = std::numeric_limits<double>::quiet_NaN();
    double springback_at_x = std::numeric_limits<double>::quiet_NaN();

    std::vector<double> Ki_array;
    UnbendingModelData model_data;
};

namespace unbending_detail {

constexpr double PI = 3.14159265358979323846;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

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

inline bool finite(double v)
{
    return std::isfinite(v);
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
    double d_theta)
{
    const double Rt_n = Rt + t / 2.0;

    return (Rt_n * std::sin(theta_rad_positive) + d_theta + Rd)
               * std::tan(theta_rad_positive)
           + Rt_n * std::cos(theta_rad_positive)
           - Rt_n;
}

inline double reverse_phase1_S1_from_theta(
    double theta_rad_positive,
    double Rd,
    double Rt,
    double t,
    double d_S1)
{
    const double Rt_n = Rt + t / 2.0;
    const double denom = std::cos(theta_rad_positive);

    if (std::abs(denom) < 1e-12) {
        return NaN;
    }

    return (Rt_n * std::sin(theta_rad_positive) + d_S1 + Rd) / denom;
}

inline double solve_reverse_theta_from_x(
    double x_target_positive,
    double theta_offset_positive_rad,
    double Rd,
    double Rt,
    double t,
    double d_theta)
{
    auto fun = [&](double th) -> double {
        return reverse_phase1_x_from_theta(th, Rd, Rt, t, d_theta)
               - x_target_positive;
    };

    double th_lo = 0.0;
    double th_hi = theta_offset_positive_rad;

    if (th_lo > th_hi) {
        std::swap(th_lo, th_hi);
    }

    const double f_lo = fun(th_lo);
    const double f_hi = fun(th_hi);

    if (finite(f_lo) && finite(f_hi) && f_lo * f_hi <= 0.0) {
        double a = th_lo;
        double b = th_hi;
        double fa = f_lo;

        for (int iter = 0; iter < 100; ++iter) {
            const double mid = 0.5 * (a + b);
            const double fm = fun(mid);

            if (!finite(fm)) {
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
        const double th = th_lo + (th_hi - th_lo) * static_cast<double>(i)
                                      / static_cast<double>(n_scan - 1);
        const double res = std::abs(fun(th));

        if (finite(res) && res < best_res) {
            best_res = res;
            best_theta = th;
        }
    }

    return best_theta;
}

inline double local_theta1_phase1(double xx, double S, double Delta_y)
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
    double t,
    double* wrapped_arc_length = nullptr)
{
    const double Rn = Rd + t / 2.0;

    const double D = std::sqrt((S - xx) * (S - xx) + Delta_y * Delta_y);
    const double delta = std::atan2(Delta_y, std::max(S - xx, 1e-12));
    const double lambda = std::acos(clamp_val(S / D, -1.0, 1.0));

    const double phi = delta - lambda;
    const double straight_len = std::sqrt(std::max(D * D - S * S, 0.0));
    const double wrap_len = Rn * phi;

    if (wrapped_arc_length) {
        *wrapped_arc_length = wrap_len;
    }

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
    auto fun = [&](double th) -> double {
        return th - std::atan2(
            Delta_y - Rt * std::sin(alpha_rad) - (t / 2.0) * std::cos(th),
            -S + xx + Rt * std::cos(alpha_rad) + (t / 2.0) * std::sin(th)
        );
    };

    const int n_guess = 200;
    const double lo = 1e-6;
    const double hi = PI - 1e-6;

    double theta_v = PI / 4.0;
    double best_res = std::numeric_limits<double>::infinity();

    for (int i = 0; i < n_guess; ++i) {
        const double th = lo + (hi - lo) * static_cast<double>(i)
                                   / static_cast<double>(n_guess - 1);

        const double res = std::abs(fun(th));

        if (finite(res) && res < best_res) {
            best_res = res;
            theta_v = th;
        }
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
    const double Rn = Rd + t / 2.0;

    const double theta_v = solve_theta_v_phase2(xx, S, Delta_y, Rt, t, alpha_rad);
    const double denom_sin = safe_sin_denom(std::sin(theta_v));

    const double L_center =
        (Delta_y - Rt * std::sin(alpha_rad) - (t / 2.0) * std::cos(theta_v))
        / denom_sin;

    const double ratio_theta_w = clamp_val(Rn / L_center, -1.0, 1.0);
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
    double alpha_rad,
    double* wrapped_arc_length = nullptr)
{
    const double Rn = Rd + t / 2.0;

    const double theta_v = solve_theta_v_phase2(xx, S, Delta_y, Rt, t, alpha_rad);
    const double denom_sin = safe_sin_denom(std::sin(theta_v));

    const double L_center =
        (Delta_y - Rt * std::sin(alpha_rad) - (t / 2.0) * std::cos(theta_v))
        / denom_sin;

    const double ratio_theta_w = clamp_val(Rn / L_center, -1.0, 1.0);
    const double theta_w = std::asin(ratio_theta_w);
    const double theta_tmp = (PI / 2.0) - theta_v + theta_w;

    const double straight_len = std::sqrt(std::max(L_center * L_center - Rn * Rn, 0.0));
    const double wrap_len = Rn * theta_tmp;

    if (wrapped_arc_length) {
        *wrapped_arc_length = wrap_len;
    }

    return straight_len + wrap_len;
}

inline double find_phase23_x_alpha(
    double S,
    double Delta_y_theta,
    double alpha_deg)
{
    const double theta_raw_0 = local_theta1_phase1(0.0, S, Delta_y_theta);

    auto fun = [&](double xx) -> double {
        return local_theta1_phase1(xx, S, Delta_y_theta)
               - theta_raw_0
               - alpha_deg;
    };

    double x_lo = 0.0;
    double x_hi = 1.0;

    while (finite(fun(x_hi)) && fun(x_hi) < 0.0 && x_hi < 1000.0) {
        x_hi = x_hi * 1.5 + 1.0;
    }

    const double f_lo = fun(x_lo);
    const double f_hi = fun(x_hi);

    if (finite(f_lo) && finite(f_hi) && f_lo * f_hi <= 0.0) {
        double a = x_lo;
        double b = x_hi;
        double fa = f_lo;

        for (int iter = 0; iter < 100; ++iter) {
            const double mid = 0.5 * (a + b);
            const double fm = fun(mid);

            if (!finite(fm)) {
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
        const double xx = x_lo + (x_hi - x_lo) * static_cast<double>(i)
                                      / static_cast<double>(n_scan - 1);

        const double res = std::abs(fun(xx));

        if (finite(res) && res < best_res) {
            best_res = res;
            best_x = xx;
        }
    }

    return best_x;
}

inline double local_Kspan_3region_scalar(
    double m,
    double Ky,
    double K2,
    double Kd,
    double My,
    double M2,
    double Md,
    double Ep,
    double t,
    double A,
    double B,
    double C)
{
    m = std::max(m, 0.0);
    const double ms = M2 / My;

    if (m <= 1.0) {
        return Ky * m;
    }

    if (m <= ms) {
        const double inside = A - B * m;
        if (inside <= 0.0) {
            return Kd;
        }

        return Ky * (1.0 / std::sqrt(inside) - C);
    }

    const double M = m * My;
    double Kval = K2 + (M - M2) / (Ep * t * t * t / 12.0);

    if (Kval > Kd) {
        Kval = Kd;
    }

    return Kval;
}

inline double curvature_from_m_3region(
    double m,
    double Ky,
    double K2,
    double Kd,
    double My,
    double M2,
    double Md,
    double Ep,
    double t,
    double A,
    double B,
    double C)
{
    return local_Kspan_3region_scalar(m, Ky, K2, Kd, My, M2, Md, Ep, t, A, B, C);
}

inline double delta_of_P_local(
    double P,
    double L,
    double Ky,
    double My,
    double K2,
    double M2,
    double Kd,
    double Md,
    double Ep,
    double t,
    double A,
    double B,
    double C,
    int Ns)
{
    if (L <= 0.0 || Ns < 2) {
        return NaN;
    }

    const double P1 = My / L;
    const double ds = L / static_cast<double>(Ns - 1);

    double integral = 0.0;

    for (int i = 0; i < Ns; ++i) {
        const double s = ds * static_cast<double>(i);
        const double wgt = L - s;

        double m = (P / P1) * (1.0 - s / L);
        m = std::min(m, Md / My);

        const double K = local_Kspan_3region_scalar(m, Ky, K2, Kd, My, M2, Md, Ep, t, A, B, C);

        const double y = wgt * K;
        const double trap_weight = (i == 0 || i == Ns - 1) ? 0.5 : 1.0;

        integral += trap_weight * y;
    }

    return integral * ds;
}

inline void solveABC_region2(double ms, double kh, double& A, double& B, double& C)
{
    auto g = [&](double u) -> double {
        const double inside = 1.0 - 2.0 * u * (ms - 1.0);
        if (inside <= 0.0 || u <= 0.0) {
            return NaN;
        }

        return 1.0 / (u * std::sqrt(inside)) - (kh + 1.0 / u - 1.0);
    };

    if (ms <= 1.0 || kh <= 1.0) {
        A = B = C = NaN;
        return;
    }

    const double u_max = 1.0 / (2.0 * (ms - 1.0));
    double u_lo = 1e-10;
    double g_lo = g(u_lo);

    if (!finite(g_lo)) {
        u_lo = 1e-8;
        g_lo = g(u_lo);
    }

    double u_hi = NaN;
    bool found = false;

    for (int p = 3; p <= 16; ++p) {
        const double eps = std::pow(10.0, -p);
        const double u_try = u_max * (1.0 - eps);
        const double g_try = g(u_try);

        if (finite(g_try) && finite(g_lo) && g_lo * g_try <= 0.0) {
            u_hi = u_try;
            found = true;
            break;
        }
    }

    if (!found) {
        const int n_scan = 5000;
        double prev_u = u_lo;
        double prev_g = g_lo;

        for (int i = 1; i < n_scan; ++i) {
            const double frac = static_cast<double>(i) / static_cast<double>(n_scan - 1);
            const double u_try = u_lo + (u_max * (1.0 - 1e-12) - u_lo) * frac;
            const double g_try = g(u_try);

            if (finite(prev_g) && finite(g_try) && prev_g * g_try <= 0.0) {
                u_lo = prev_u;
                u_hi = u_try;
                found = true;
                break;
            }

            prev_u = u_try;
            prev_g = g_try;
        }
    }

    if (!found) {
        A = B = C = NaN;
        return;
    }

    double a = u_lo;
    double b = u_hi;
    double fa = g(a);

    for (int iter = 0; iter < 100; ++iter) {
        const double mid = 0.5 * (a + b);
        const double fm = g(mid);

        if (!finite(fm)) {
            break;
        }

        if (std::abs(fm) < 1e-12 || std::abs(b - a) < 1e-12) {
            const double u = mid;
            B = 2.0 * u * u * u;
            C = 1.0 / u - 1.0;
            A = u * u + 2.0 * u * u * u;
            return;
        }

        if (fa * fm <= 0.0) {
            b = mid;
        } else {
            a = mid;
            fa = fm;
        }
    }

    const double u = 0.5 * (a + b);
    B = 2.0 * u * u * u;
    C = 1.0 / u - 1.0;
    A = u * u + 2.0 * u * u * u;
}

inline void get_transition_displacements_local(
    double L,
    double Ky,
    double My,
    double K2,
    double M2,
    double Kd,
    double Md,
    double Ep,
    double t,
    double A,
    double B,
    double C,
    double& delta_R2,
    double& delta_sat)
{
    const int Ns = 1000;

    const double P2 = M2 / L;
    const double Ps = Md / L;

    delta_R2 = delta_of_P_local(P2, L, Ky, My, K2, M2, Kd, Md, Ep, t, A, B, C, Ns);
    delta_sat = delta_of_P_local(Ps, L, Ky, My, K2, M2, Kd, Md, Ep, t, A, B, C, Ns);
}

inline double solve_P_for_displacement(
    double x,
    double L,
    double Ps,
    double Ky,
    double My,
    double K2,
    double M2,
    double Kd,
    double Md,
    double Ep,
    double t,
    double A,
    double B,
    double C)
{
    const int Ns = 300;

    auto Froot = [&](double P) -> double {
        return delta_of_P_local(P, L, Ky, My, K2, M2, Kd, Md, Ep, t, A, B, C, Ns) - x;
    };

    const double f0 = Froot(0.0);
    const double f1 = Froot(Ps);

    if (finite(f0) && finite(f1) && f0 * f1 <= 0.0) {
        double a = 0.0;
        double b = Ps;
        double fa = f0;

        for (int iter = 0; iter < 80; ++iter) {
            const double mid = 0.5 * (a + b);
            const double fm = Froot(mid);

            if (!finite(fm)) {
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

    const int n_grid = 3000;
    double best_P = 0.0;
    double best_res = std::numeric_limits<double>::infinity();

    for (int i = 0; i < n_grid; ++i) {
        const double P = Ps * static_cast<double>(i) / static_cast<double>(n_grid - 1);
        const double res = std::abs(Froot(P));

        if (finite(res) && res < best_res) {
            best_res = res;
            best_P = P;
        }
    }

    return best_P;
}

inline double linear_interpolate(
    double x0,
    double y0,
    double x1,
    double y1,
    double xq)
{
    if (std::abs(x1 - x0) < 1e-12) {
        return y0;
    }

    const double a = (xq - x0) / (x1 - x0);
    return y0 + a * (y1 - y0);
}

} // namespace unbending_detail


inline DisplacementUnbendingResult Displacement_from_Unloaded_Angle_Unbending(
    double E,
    double nu,
    double sigma_0,
    double Ep,
    double f,
    double Rd,
    double Rt,
    double t,
    double w,
    double d,
    double alpha,
    const std::vector<double>& x_range,
    double curvature_factor,
    double ep_curvature_coeff,
    double ep_moment_coeff,
    double extended_length,
    double initial_angle_offset,
    double Desired_angle)
{
    using namespace unbending_detail;

    DisplacementUnbendingResult out;

    const int nx = static_cast<int>(x_range.size());

    if (nx == 0) {
        return out;
    }

    if (initial_angle_offset >= 0.0) {
        throw std::invalid_argument("initial_angle_offset must be negative for unbending.");
    }

    const double S = Rd + Rt + t;
    const double Rn = Rd + t / 2.0;

    const double alpha_rad = deg2rad(alpha);

    const double d_theta = d * curvature_factor;
    const double d_S1 = d;

    const double Delta_y_theta = Rd + d_theta;
    const double Delta_y_S1 = Rd + d_S1;

    const double theta_offset_abs_deg = std::abs(initial_angle_offset);
    const double theta_offset_abs_rad = deg2rad(theta_offset_abs_deg);

    const double x_offset_positive =
        reverse_phase1_x_from_theta(theta_offset_abs_rad, Rd, Rt, t, d_theta);

    const double x_offset_geom = -x_offset_positive;

    const double theta_shift = theta_offset_abs_deg;

    const double S1_reverse_end =
        reverse_phase1_S1_from_theta(0.0, Rd, Rt, t, d_S1);

    const double theta_raw_phase2_0 =
        local_theta1_phase1(0.0, S, Delta_y_theta);

    const double S1_raw_phase2_0 =
        local_S1_phase1(0.0, S, Delta_y_S1, Rd, t);

    const double S1_phase2_offset =
        S1_reverse_end - S1_raw_phase2_0;

    const double x_alpha_geom =
        find_phase23_x_alpha(S, Delta_y_theta, alpha);

    const double theta_raw_phase3_alpha =
        local_theta1_phase2(x_alpha_geom, S, Delta_y_theta, Rt, Rd, t, alpha_rad);

    const double theta_phase3_offset =
        alpha - theta_raw_phase3_alpha;

    const double S1_phase2_alpha_raw =
        local_S1_phase1(x_alpha_geom, S, Delta_y_S1, Rd, t);

    const double S1_phase2_alpha =
        S1_phase2_alpha_raw + S1_phase2_offset;

    const double S1_phase3_alpha_raw =
        local_S1_phase2(x_alpha_geom, S, Delta_y_S1, Rt, Rd, t, alpha_rad);

    const double S1_phase3_offset =
        S1_phase2_alpha - S1_phase3_alpha_raw;

    out.model_data.x = x_range;
    out.model_data.theta_1.assign(nx, NaN);
    out.model_data.theta_2.assign(nx, NaN);
    out.model_data.Springback.assign(nx, NaN);
    out.model_data.Clamped_Moment.assign(nx, NaN);
    out.model_data.Clamped_Curvature.assign(nx, NaN);
    out.model_data.Force.assign(nx, NaN);
    out.model_data.Trial_Force.assign(nx, NaN);
    out.model_data.S1.assign(nx, NaN);
    out.model_data.Sy.assign(nx, 0.0);
    out.model_data.Wrapped_arc_length.assign(nx, NaN);
    out.model_data.delta_elastic_local.assign(nx, NaN);
    out.model_data.delta_R2_local.assign(nx, NaN);
    out.model_data.delta_plastic_local.assign(nx, NaN);
    out.model_data.Ki_array.assign(nx, NaN);
    out.model_data.phase_flag.assign(nx, 0);

    for (int k = 0; k < nx; ++k) {
        const double x_input = x_range[k];

        if (!finite(x_input) || x_input < 0.0) {
            continue;
        }

        const double x_geom = x_input + x_offset_geom;

        if (x_geom <= 0.0) {
            out.model_data.phase_flag[k] = 1;

            const double x_phase1_positive = -x_geom;

            const double theta_geom_rad =
                solve_reverse_theta_from_x(
                    x_phase1_positive,
                    theta_offset_abs_rad,
                    Rd,
                    Rt,
                    t,
                    d_theta
                );

            const double theta_geom_deg = rad2deg(theta_geom_rad);

            out.model_data.theta_1[k] = theta_offset_abs_deg - theta_geom_deg;

            out.model_data.S1[k] =
                reverse_phase1_S1_from_theta(theta_geom_rad, Rd, Rt, t, d_S1)
                + extended_length;

            out.model_data.Wrapped_arc_length[k] = 0.0;
        }
        else if (x_geom <= x_alpha_geom) {
            out.model_data.phase_flag[k] = 2;

            const double theta_raw_phase2_x =
                local_theta1_phase1(x_geom, S, Delta_y_theta);

            const double theta_phase2_data =
                theta_raw_phase2_x - theta_raw_phase2_0;

            out.model_data.theta_1[k] =
                theta_phase2_data + theta_shift;

            double wrap_len = 0.0;

            const double S1_phase2_data =
                local_S1_phase1(x_geom, S, Delta_y_S1, Rd, t, &wrap_len);

            out.model_data.S1[k] =
                S1_phase2_data + S1_phase2_offset + extended_length;

            out.model_data.Wrapped_arc_length[k] = wrap_len;
        }
        else {
            out.model_data.phase_flag[k] = 3;

            const double theta_phase3_data =
                local_theta1_phase2(x_geom, S, Delta_y_theta, Rt, Rd, t, alpha_rad);

            out.model_data.theta_1[k] =
                theta_phase3_data + theta_phase3_offset + theta_shift;

            double wrap_len = 0.0;

            const double S1_phase3_data =
                local_S1_phase2(x_geom, S, Delta_y_S1, Rt, Rd, t, alpha_rad, &wrap_len);

            out.model_data.S1[k] =
                S1_phase3_data + S1_phase3_offset + extended_length;

            out.model_data.Wrapped_arc_length[k] = wrap_len;
        }
    }

    const double E_prime = E / (1.0 - nu * nu);
    const double I = t * t * t / 12.0;

    const double Ky = (2.0 * sigma_0) / (E * t);
    const double My = (sigma_0 * t * t) / 6.0;

    const double Kp = ep_curvature_coeff * Ky;
    const double Mp = ep_moment_coeff * My;

    const double Kd = (1.0 / (Rd + t / 2.0));
    const double Md = (Ep * t * t * t / 12.0) * (Kd - Kp) + Mp;

    double A_r2 = NaN;
    double B_r2 = NaN;
    double C_r2 = NaN;

    solveABC_region2(Mp / My, Kp / Ky, A_r2, B_r2, C_r2);

    const double Ld = d_theta + Rd + extended_length;

    for (int k = 0; k < nx; ++k) {
        const double xk = x_range[k];
        const double Lk = out.model_data.S1[k];
        const double theta1 = out.model_data.theta_1[k];

        if (!finite(xk) || !finite(Lk) || !finite(theta1) || Lk <= 1e-9 || theta1 < 0.0) {
            continue;
        }

        out.model_data.delta_elastic_local[k] = My * Lk * Lk / (3.0 * E * I);

        double delta_R2_k = NaN;
        double delta_sat_k = NaN;

        get_transition_displacements_local(
            Lk,
            Ky,
            My,
            Kp,
            Mp,
            Kd,
            Md,
            Ep,
            t,
            A_r2,
            B_r2,
            C_r2,
            delta_R2_k,
            delta_sat_k
        );

        out.model_data.delta_R2_local[k] = delta_R2_k;
        out.model_data.delta_plastic_local[k] = delta_sat_k;

        const double delta_e_k = out.model_data.delta_elastic_local[k];

        if (xk <= delta_e_k) {
            out.model_data.Clamped_Moment[k] =
                3.0 * E * I * xk / (Lk * Lk);

            out.model_data.Clamped_Curvature[k] =
                out.model_data.Clamped_Moment[k] / (E * I);
        }
        else {
            if (xk >= delta_sat_k) {
                out.model_data.Clamped_Moment[k] = Md;
                out.model_data.Clamped_Curvature[k] = Kd;
            }
            else {
                const double Ps = Md / Lk;

                const double P_sol =
                    solve_P_for_displacement(
                        xk,
                        Lk,
                        Ps,
                        Ky,
                        My,
                        Kp,
                        Mp,
                        Kd,
                        Md,
                        Ep,
                        t,
                        A_r2,
                        B_r2,
                        C_r2
                    );

                out.model_data.Clamped_Moment[k] =
                    Lk * clamp_val(P_sol, 0.0, Ps);

                const double m0 = out.model_data.Clamped_Moment[k] / My;

                double K0 =
                    curvature_from_m_3region(
                        m0,
                        Ky,
                        Kp,
                        Kd,
                        My,
                        Mp,
                        Md,
                        Ep,
                        t,
                        A_r2,
                        B_r2,
                        C_r2
                    );

                if (K0 >= Kd || out.model_data.Clamped_Moment[k] >= Md) {
                    K0 = Kd;
                    out.model_data.Clamped_Moment[k] = Md;
                }

                out.model_data.Clamped_Curvature[k] = K0;
            }
        }

        const double M = out.model_data.Clamped_Moment[k];

        out.model_data.Springback[k] =
            rad2deg((6.0 * M * Lk) / (E_prime * t * t * t));

        out.model_data.theta_2[k] =
            theta1 - out.model_data.Springback[k];

        const double theta_rad = deg2rad(theta1);

        const double denom =
            Ld *
            (
                (1.0 - (Rt / Ld) * std::sin(theta_rad))
                +
                ((std::tan(theta_rad) - f) / (1.0 + f * std::tan(theta_rad)))
                *
                (
                    (xk / Ld)
                    - (Rt / Ld) * (1.0 - std::cos(theta_rad))
                )
            );

        if (std::abs(denom) >= 1e-12) {
            out.model_data.Force[k] = (M * w) / denom;
        }

        if (std::abs(Lk) >= 1e-12) {
            out.model_data.Trial_Force[k] = (M / Lk) * w;
        }

        if (M > 1e-12) {
            out.model_data.Sy[k] = Lk * (1.0 - My / M);
            if (out.model_data.Sy[k] < 0.0) {
                out.model_data.Sy[k] = 0.0;
            }
        }
    }

    std::vector<double> system_gain(nx, NaN);
    out.Ki_array.assign(nx, NaN);

    if (nx >= 2) {
        for (int k = 0; k < nx; ++k) {
            int i0 = std::max(0, k - 1);
            int i1 = std::min(nx - 1, k + 1);

            if (k == 0) {
                i0 = 0;
                i1 = 1;
            }
            else if (k == nx - 1) {
                i0 = nx - 2;
                i1 = nx - 1;
            }

            const double dx = x_range[i1] - x_range[i0];
            const double dtheta = out.model_data.theta_2[i1] - out.model_data.theta_2[i0];

            if (finite(dx) && finite(dtheta) && std::abs(dx) > 1e-12) {
                system_gain[k] = dtheta / dx;

                if (std::abs(system_gain[k]) > 1e-12) {
                    out.Ki_array[k] = 1.0 / system_gain[k];
                }
            }
        }
    }

    out.model_data.Ki_array = out.Ki_array;

    bool found_target = false;

    for (int k = 0; k < nx - 1; ++k) {
        const double y0 = out.model_data.theta_2[k];
        const double y1 = out.model_data.theta_2[k + 1];

        if (!finite(y0) || !finite(y1)) {
            continue;
        }

        const double f0 = y0 - Desired_angle;
        const double f1 = y1 - Desired_angle;

        if (f0 == 0.0) {
            out.x_required = x_range[k];
            out.theta_1_at_x = out.model_data.theta_1[k];
            out.springback_at_x = out.model_data.Springback[k];
            found_target = true;
            break;
        }

        if (f0 * f1 <= 0.0) {
            out.x_required =
                linear_interpolate(y0, x_range[k], y1, x_range[k + 1], Desired_angle);

            out.theta_1_at_x =
                linear_interpolate(x_range[k], out.model_data.theta_1[k],
                                   x_range[k + 1], out.model_data.theta_1[k + 1],
                                   out.x_required);

            out.springback_at_x =
                linear_interpolate(x_range[k], out.model_data.Springback[k],
                                   x_range[k + 1], out.model_data.Springback[k + 1],
                                   out.x_required);

            found_target = true;
            break;
        }
    }

    if (!found_target) {
        double best_err = std::numeric_limits<double>::infinity();
        int best_idx = -1;

        for (int k = 0; k < nx; ++k) {
            const double y = out.model_data.theta_2[k];

            if (!finite(y)) {
                continue;
            }

            const double err = std::abs(y - Desired_angle);

            if (err < best_err) {
                best_err = err;
                best_idx = k;
            }
        }

        if (best_idx >= 0) {
            out.x_required = x_range[best_idx];
            out.theta_1_at_x = out.model_data.theta_1[best_idx];
            out.springback_at_x = out.model_data.Springback[best_idx];
        }
    }

    return out;
}

} // namespace springback_model