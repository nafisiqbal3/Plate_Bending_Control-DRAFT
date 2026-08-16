#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace springback_model {

struct ModelData {
    std::vector<double> x;
    std::vector<double> theta_1;
    std::vector<double> theta_2;
    std::vector<double> Springback;
    std::vector<double> Clamped_Moment;
    std::vector<double> Clamped_Curvature;
    std::vector<double> Force;
    std::vector<double> S1;
    std::vector<double> integral_controller_gain;
    double delta_elastic_final = std::numeric_limits<double>::quiet_NaN();
    double delta_R2_final = std::numeric_limits<double>::quiet_NaN();
    double My = std::numeric_limits<double>::quiet_NaN();
    double Ky = std::numeric_limits<double>::quiet_NaN();
    double Md = std::numeric_limits<double>::quiet_NaN();
    double Kd = std::numeric_limits<double>::quiet_NaN();
};

struct DisplacementFromUnloadedAngleResult {
    double x_required = std::numeric_limits<double>::quiet_NaN();
    double theta_1_at_x = std::numeric_limits<double>::quiet_NaN();
    double springback_at_x = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> integral_controller_gain;
    ModelData model_data;
};

namespace disp_detail {

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

inline bool is_finite(double v)
{
    return std::isfinite(v);
}

inline double trapz(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    if (x.size() != y.size() || x.size() < 2)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum = 0.0;

    for (size_t i = 1; i < x.size(); ++i)
    {
        sum += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    }

    return sum;
}

inline double linear_interp(
    double xq,
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    if (x.size() != y.size() || x.empty())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (x.size() == 1)
    {
        return y.front();
    }

    if (xq <= x.front())
    {
        return y.front();
    }

    if (xq >= x.back())
    {
        return y.back();
    }

    auto it = std::lower_bound(x.begin(), x.end(), xq);
    const size_t i1 = static_cast<size_t>(it - x.begin());
    const size_t i0 = i1 - 1;

    const double x0 = x[i0];
    const double x1 = x[i1];
    const double y0 = y[i0];
    const double y1 = y[i1];

    if (std::abs(x1 - x0) < 1e-14)
    {
        return y0;
    }

    const double q = (xq - x0) / (x1 - x0);
    return y0 + q * (y1 - y0);
}

inline double nearest_x_for_target(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double target)
{
    double best_err = std::numeric_limits<double>::infinity();
    double best_x = std::numeric_limits<double>::quiet_NaN();

    for (size_t i = 0; i < x.size(); ++i)
    {
        if (!is_finite(y[i]))
        {
            continue;
        }

        const double err = std::abs(y[i] - target);

        if (err < best_err)
        {
            best_err = err;
            best_x = x[i];
        }
    }

    return best_x;
}

inline double invert_y_to_x_linear(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double target)
{
    for (size_t i = 1; i < x.size(); ++i)
    {
        if (!is_finite(y[i - 1]) || !is_finite(y[i]))
        {
            continue;
        }

        const double f0 = y[i - 1] - target;
        const double f1 = y[i] - target;

        if (std::abs(f0) < 1e-14)
        {
            return x[i - 1];
        }

        if (std::abs(f1) < 1e-14)
        {
            return x[i];
        }

        if ((f0 < 0.0 && f1 > 0.0) || (f0 > 0.0 && f1 < 0.0))
        {
            const double denom = y[i] - y[i - 1];

            if (std::abs(denom) < 1e-14)
            {
                return x[i - 1];
            }

            const double q = (target - y[i - 1]) / denom;
            return x[i - 1] + q * (x[i] - x[i - 1]);
        }
    }

    return nearest_x_for_target(x, y, target);
}

inline std::vector<double> calculate_integral_controller_gain(
    const std::vector<double>& x,
    const std::vector<double>& theta_2)
{
    const size_t n = x.size();

    std::vector<double> gain(
        n,
        std::numeric_limits<double>::quiet_NaN()
    );

    if (n != theta_2.size() || n < 2)
    {
        return gain;
    }

    for (size_t i = 0; i < n; ++i)
    {
        size_t i0 = 0;
        size_t i1 = 0;

        if (i == 0)
        {
            i0 = 0;
            i1 = 1;
        }
        else if (i == n - 1)
        {
            i0 = n - 2;
            i1 = n - 1;
        }
        else
        {
            i0 = i - 1;
            i1 = i + 1;
        }

        if (!is_finite(x[i0]) || !is_finite(x[i1]) ||
            !is_finite(theta_2[i0]) || !is_finite(theta_2[i1]))
        {
            gain[i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        const double delta_x = x[i1] - x[i0];
        const double delta_theta2 = theta_2[i1] - theta_2[i0];

        if (std::abs(delta_theta2) < 1e-14)
        {
            gain[i] = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            gain[i] = delta_x / delta_theta2;
        }
    }

    return gain;
}

inline double bracketed_bisection(
    const std::function<double(double)>& f,
    double a,
    double b,
    int max_iter = 120,
    double tol = 1e-10)
{
    double fa = f(a);
    double fb = f(b);

    if (!is_finite(fa) || !is_finite(fb))
    {
        throw std::runtime_error("Non-finite function value in bracketed_bisection.");
    }

    if (fa == 0.0)
    {
        return a;
    }

    if (fb == 0.0)
    {
        return b;
    }

    if (fa * fb > 0.0)
    {
        throw std::runtime_error("Root not bracketed in bracketed_bisection.");
    }

    double left = a;
    double right = b;
    double fleft = fa;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        double mid = 0.5 * (left + right);
        double fm = f(mid);

        if (!is_finite(fm))
        {
            mid = std::nextafter(mid, right);
            fm = f(mid);
        }

        if (std::abs(fm) < tol || std::abs(right - left) < tol)
        {
            return mid;
        }

        if (fleft * fm <= 0.0)
        {
            right = mid;
        }
        else
        {
            left = mid;
            fleft = fm;
        }
    }

    return 0.5 * (left + right);
}

inline double local_theta1_phase1_disp(
    double xx,
    double S,
    double Delta_y)
{
    const double D = std::sqrt((S - xx) * (S - xx) + Delta_y * Delta_y);
    const double delta = std::atan2(Delta_y, std::max(S - xx, 1e-12));
    const double lambda = std::acos(clamp_val(S / D, -1.0, 1.0));

    return rad2deg(delta - lambda);
}

inline double solve_theta_v_phase2_disp(
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

    const int n_guess = 7;
    const double lo = 1e-4;
    const double hi = PI / 2.0 - 1e-4;

    std::vector<double> guess_list;
    guess_list.reserve(n_guess);

    for (int i = 0; i < n_guess; ++i)
    {
        guess_list.push_back(
            lo + (hi - lo) * static_cast<double>(i) /
                     static_cast<double>(n_guess - 1)
        );
    }

    double theta_v = std::numeric_limits<double>::quiet_NaN();
    double best_res = std::numeric_limits<double>::infinity();

    for (int j = 0; j < static_cast<int>(guess_list.size()) - 1; ++j)
    {
        const double a = guess_list[j];
        const double b = guess_list[j + 1];

        const double fa = fun(a);
        const double fb = fun(b);

        if (!is_finite(fa) || !is_finite(fb))
        {
            continue;
        }

        try
        {
            double th_try = std::numeric_limits<double>::quiet_NaN();

            if (fa == 0.0)
            {
                th_try = a;
            }
            else if (fb == 0.0)
            {
                th_try = b;
            }
            else if (fa * fb < 0.0)
            {
                th_try = bracketed_bisection(fun, a, b, 100, 1e-11);
            }
            else
            {
                continue;
            }

            const double res = std::abs(fun(th_try));

            if (is_finite(th_try) && res < best_res)
            {
                theta_v = th_try;
                best_res = res;
            }
        }
        catch (...)
        {
        }
    }

    if (!is_finite(theta_v))
    {
        for (double g : guess_list)
        {
            const double res = std::abs(fun(g));

            if (is_finite(res) && res < best_res)
            {
                theta_v = g;
                best_res = res;
            }
        }
    }

    if (!is_finite(theta_v))
    {
        theta_v = 0.0;
    }

    return theta_v;
}

inline double local_theta1_phase2_disp(
    double xx,
    double S,
    double Delta_y,
    double Rt,
    double Rd,
    double t,
    double alpha_rad)
{
    const double theta_v =
        solve_theta_v_phase2_disp(xx, S, Delta_y, Rt, t, alpha_rad);

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

inline double local_safe_g(
    const std::function<double(double)>& g,
    double u)
{
    try
    {
        const double val = g(u);
        return is_finite(val) ? val : std::numeric_limits<double>::quiet_NaN();
    }
    catch (...)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

inline void solveABC_region2(
    double ms,
    double kh,
    double& A,
    double& B,
    double& C)
{
    auto g = [&](double u) -> double
    {
        return 1.0 / (u * std::sqrt(1.0 - 2.0 * u * (ms - 1.0)))
               - (kh + 1.0 / u - 1.0);
    };

    const double u_max = 1.0 / (2.0 * (ms - 1.0));

    double u_lo = 1e-10;
    double g_lo = g(u_lo);

    if (!is_finite(g_lo))
    {
        u_lo = 1e-8;
        g_lo = g(u_lo);
    }

    bool found = false;
    double u_hi = std::numeric_limits<double>::quiet_NaN();

    for (int p = 3; p <= 16; ++p)
    {
        const double u_try = u_max * (1.0 - std::pow(10.0, -p));
        const double g_try = g(u_try);

        if (is_finite(g_try)
            && ((g_try > 0.0 && g_lo < 0.0)
                || (g_try < 0.0 && g_lo > 0.0)))
        {
            u_hi = u_try;
            found = true;
            break;
        }
    }

    if (!found)
    {
        const int N = 5000;

        std::vector<double> us(N);
        std::vector<double> gv(N);

        for (int i = 0; i < N; ++i)
        {
            const double expo =
                -16.0 + (-3.0 + 16.0) * static_cast<double>(i) /
                            static_cast<double>(N - 1);

            const double val = std::pow(10.0, expo);

            us[i] = u_max * (1.0 - val);
            gv[i] = local_safe_g(g, us[i]);
        }

        int j_found = -1;

        for (int j = 0; j < N - 1; ++j)
        {
            if (!is_finite(gv[j]) || !is_finite(gv[j + 1]))
            {
                continue;
            }

            if (gv[j] == 0.0 || gv[j] * gv[j + 1] < 0.0)
            {
                j_found = j;
                break;
            }
        }

        if (j_found < 0)
        {
            throw std::runtime_error("Could not bracket u in solveABC_region2.");
        }

        u_lo = us[j_found];
        u_hi = us[j_found + 1];
    }

    const double u = bracketed_bisection(g, u_lo, u_hi, 150, 1e-12);

    B = 2.0 * u * u * u;
    C = 1.0 / u - 1.0;
    A = u * u + 2.0 * u * u * u;
}

inline double curvature_from_m_3region(
    double m,
    double Ky,
    double K2,
    double Kd,
    double My,
    double M2,
    double Ep,
    double t,
    double A,
    double B,
    double C)
{
    m = std::max(m, 0.0);

    const double ms = M2 / My;

    if (m <= 1.0)
    {
        return Ky * m;
    }
    else if (m <= ms)
    {
        return Ky * (1.0 / std::sqrt(A - B * m) - C);
    }
    else
    {
        const double M = m * My;
        double Kval = K2 + (M - M2) / (Ep * std::pow(t, 3) / 12.0);

        if (Kval > Kd)
        {
            Kval = Kd;
        }

        return Kval;
    }
}

inline std::vector<double> local_Kspan_3region(
    const std::vector<double>& mvec,
    double Ky,
    double K2,
    double Kd,
    double My,
    double M2,
    double Ep,
    double t,
    double A,
    double B,
    double C)
{
    std::vector<double> Kspan(mvec.size(), 0.0);

    const double ms = M2 / My;

    for (size_t i = 0; i < mvec.size(); ++i)
    {
        const double m = mvec[i];

        if (m <= 1.0)
        {
            Kspan[i] = Ky * m;
        }
        else if (m <= ms)
        {
            Kspan[i] = Ky * (1.0 / std::sqrt(A - B * m) - C);
        }
        else
        {
            const double M = m * My;

            Kspan[i] =
                K2 + (M - M2) / (Ep * std::pow(t, 3) / 12.0);

            Kspan[i] = std::min(Kspan[i], Kd);
        }
    }

    return Kspan;
}

inline std::pair<double, double> get_transition_displacements_local(
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
    double C)
{
    const int Ns = 3000;

    std::vector<double> s(Ns);
    std::vector<double> wgt(Ns);

    for (int i = 0; i < Ns; ++i)
    {
        s[i] = L * static_cast<double>(i) / static_cast<double>(Ns - 1);
        wgt[i] = L - s[i];
    }

    const double P1 = My / L;
    const double P2 = M2 / L;
    const double Ps = Md / L;

    auto delta_of_P = [&](double P) -> double
    {
        std::vector<double> mvec(Ns);

        for (int i = 0; i < Ns; ++i)
        {
            mvec[i] =
                std::min((P / P1) * (1.0 - s[i] / L), Md / My);
        }

        const std::vector<double> Kspan =
            local_Kspan_3region(
                mvec, Ky, K2, Kd, My, M2, Ep, t, A, B, C
            );

        std::vector<double> integrand(Ns);

        for (int i = 0; i < Ns; ++i)
        {
            integrand[i] = wgt[i] * Kspan[i];
        }

        return trapz(s, integrand);
    };

    const double delta_R2 = delta_of_P(P2);
    const double delta_sat = delta_of_P(Ps);

    return {delta_R2, delta_sat};
}

}


inline DisplacementFromUnloadedAngleResult Displacement_from_Unloaded_Angle(
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
    double Desired_angle)
{
    if (x_range.empty())
    {
        throw std::invalid_argument("x_range must not be empty.");
    }

    DisplacementFromUnloadedAngleResult out;
    ModelData& model_data = out.model_data;

    const size_t nx = x_range.size();
    const double alpha_rad = disp_detail::deg2rad(alpha);

    const double S = Rd + Rt + t;

    const double d_S1 = d;
    const double d_theta = d * curvature_factor;

    const double Delta_y_S1 = Rd + d_S1;
    const double Delta_y_theta = Rd + d_theta;

    const double E_prime = E / (1.0 - nu * nu);
    const double I = std::pow(t, 3) / 12.0;

    const double Ky = (2.0 * sigma_0) / (E * t);
    const double My = (sigma_0 * t * t) / 6.0;

    const double Kp = ep_curvature_coeff * Ky;
    const double Mp = ep_moment_coeff * My;

    const double Kd = (1.0 / (Rd + t / 2.0));
    const double Md =
        (Ep * std::pow(t, 3) / 12.0) * (Kd - Kp) + Mp;

    std::vector<double> theta_1(nx, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> theta_2(nx, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> Springback(nx, 0.0);

    std::vector<double> Wrapped_arc_length(nx, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> S1(nx, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> S1_only(nx, std::numeric_limits<double>::quiet_NaN());

    std::vector<double> Clamped_Moment(nx, 0.0);
    std::vector<double> Clamped_Curvature(nx, 0.0);
    std::vector<double> Force(nx, 0.0);

    std::vector<double> delta_elastic_local(nx, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> delta_R2_local(nx, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> delta_plastic_local(nx, std::numeric_limits<double>::quiet_NaN());

    const double Ld = d_theta + Rd + extended_length;

    std::vector<double> theta_phase1_S1_vals(nx);

    for (size_t i = 0; i < nx; ++i)
    {
        theta_phase1_S1_vals[i] =
            disp_detail::local_theta1_phase1_disp(
                x_range[i], S, Delta_y_S1
            );
    }

    double x_alpha_S1 = std::numeric_limits<double>::quiet_NaN();

    const double theta1_x0_S1 = theta_phase1_S1_vals.front();
    const double theta1_xend_S1 = theta_phase1_S1_vals.back();

    if (disp_detail::is_finite(theta1_x0_S1)
        && disp_detail::is_finite(theta1_xend_S1)
        && theta1_x0_S1 <= alpha
        && theta1_xend_S1 >= alpha)
    {
        auto froot = [&](double xx) -> double
        {
            return disp_detail::local_theta1_phase1_disp(
                       xx, S, Delta_y_S1
                   )
                   - alpha;
        };

        try
        {
            x_alpha_S1 =
                disp_detail::bracketed_bisection(
                    froot,
                    x_range.front(),
                    x_range.back(),
                    120,
                    1e-10
                );
        }
        catch (...)
        {
            x_alpha_S1 =
                disp_detail::nearest_x_for_target(
                    x_range, theta_phase1_S1_vals, alpha
                );
        }
    }
    else
    {
        x_alpha_S1 =
            disp_detail::nearest_x_for_target(
                x_range, theta_phase1_S1_vals, alpha
            );
    }

    std::vector<double> theta_phase1_theta_vals(nx);

    for (size_t i = 0; i < nx; ++i)
    {
        theta_phase1_theta_vals[i] =
            disp_detail::local_theta1_phase1_disp(
                x_range[i], S, Delta_y_theta
            );
    }

    double x_alpha_theta = std::numeric_limits<double>::quiet_NaN();

    const double theta1_x0_theta = theta_phase1_theta_vals.front();
    const double theta1_xend_theta = theta_phase1_theta_vals.back();

    if (disp_detail::is_finite(theta1_x0_theta)
        && disp_detail::is_finite(theta1_xend_theta)
        && theta1_x0_theta <= alpha
        && theta1_xend_theta >= alpha)
    {
        auto froot = [&](double xx) -> double
        {
            return disp_detail::local_theta1_phase1_disp(
                       xx, S, Delta_y_theta
                   )
                   - alpha;
        };

        try
        {
            x_alpha_theta =
                disp_detail::bracketed_bisection(
                    froot,
                    x_range.front(),
                    x_range.back(),
                    120,
                    1e-10
                );
        }
        catch (...)
        {
            x_alpha_theta =
                disp_detail::nearest_x_for_target(
                    x_range, theta_phase1_theta_vals, alpha
                );
        }
    }
    else
    {
        x_alpha_theta =
            disp_detail::nearest_x_for_target(
                x_range, theta_phase1_theta_vals, alpha
            );
    }

    // S1 calculation

    for (size_t k = 0; k < nx; ++k)
    {
        const double xk = x_range[k];

        if (xk <= x_alpha_S1)
        {
            const double D =
                std::sqrt((S - xk) * (S - xk) + Delta_y_S1 * Delta_y_S1);

            const double delta =
                std::atan2(Delta_y_S1, std::max(S - xk, 1e-12));

            const double lambda =
                std::acos(disp_detail::clamp_val(S / D, -1.0, 1.0));

            const double straight_len =
                std::sqrt(std::max(D * D - S * S, 0.0));

            const double wrap_len =
                (Rd + t / 2.0) * (delta - lambda);

            S1_only[k] = straight_len + wrap_len;
            Wrapped_arc_length[k] = wrap_len;
        }
        else
        {
            const double theta_v_S1 =
                disp_detail::solve_theta_v_phase2_disp(
                    xk, S, Delta_y_S1, Rt, t, alpha_rad
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
                        Rn * Rn - std::pow(Rd + t / 2.0, 2),
                        0.0
                    )
                );

            const double theta_w =
                std::asin(
                    disp_detail::clamp_val(
                        (Rd + t / 2.0) / Rn,
                        -1.0,
                        1.0
                    )
                );

            const double theta1_tmp =
                (disp_detail::PI / 2.0) - theta_v_S1 + theta_w;

            const double wrap_len =
                (Rd + t / 2.0) * theta1_tmp;

            S1_only[k] = straight_len + wrap_len;
            Wrapped_arc_length[k] = wrap_len;
        }

        S1[k] = S1_only[k] + extended_length;
    }

    // ----------------------------------------------------
    // theta_1 calculation
    // ----------------------------------------------------
    for (size_t k = 0; k < nx; ++k)
    {
        const double xk = x_range[k];

        if (xk <= x_alpha_theta)
        {
            theta_1[k] =
                disp_detail::local_theta1_phase1_disp(
                    xk, S, Delta_y_theta
                );
        }
        else
        {
            theta_1[k] =
                disp_detail::local_theta1_phase2_disp(
                    xk, S, Delta_y_theta, Rt, Rd, t, alpha_rad
                );
        }
    }

    std::vector<bool> valid(nx, false);

    int k_max = -1;
    int k_ref = -1;

    for (size_t k = 0; k < nx; ++k)
    {
        valid[k] =
            disp_detail::is_finite(S1[k])
            && disp_detail::is_finite(theta_1[k])
            && S1[k] > 1e-9
            && theta_1[k] >= 0.0;

        if (valid[k] && k_ref < 0)
        {
            k_ref = static_cast<int>(k);
        }

        if (valid[k])
        {
            k_max = static_cast<int>(k);
        }
    }

    if (k_max < 0)
    {
        throw std::runtime_error(
            "No valid theta_1/S1 computed. Check geometry and phase transition."
        );
    }

    // ----------------------------------------------------
    // Region-2 nonlinear coefficients
    // ----------------------------------------------------
    double A_r2 = 0.0;
    double B_r2 = 0.0;
    double C_r2 = 0.0;

    disp_detail::solveABC_region2(
        Mp / My,
        Kp / Ky,
        A_r2,
        B_r2,
        C_r2
    );

    // ----------------------------------------------------
    // Local transition displacements
    // ----------------------------------------------------
    for (int k = 0; k <= k_max; ++k)
    {
        const size_t kk = static_cast<size_t>(k);
        const double Lk = S1[kk];

        delta_elastic_local[kk] =
            My * Lk * Lk / (3.0 * E * I);

        const auto transition =
            disp_detail::get_transition_displacements_local(
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
                C_r2
            );

        delta_R2_local[kk] = transition.first;
        delta_plastic_local[kk] = transition.second;
    }

    const double delta_elastic_final =
        delta_elastic_local[static_cast<size_t>(k_max)];

    const double delta_R2_final =
        delta_R2_local[static_cast<size_t>(k_max)];

    const int Nx_int = 300;

    for (int k = 0; k <= k_max; ++k)
    {
        const size_t kk = static_cast<size_t>(k);

        if (!valid[kk])
        {
            continue;
        }

        const double xk = x_range[kk];
        const double Lk = S1[kk];

        const double delta_e_k = delta_elastic_local[kk];
        const double delta_sat_k = delta_plastic_local[kk];

        std::vector<double> s_int(Nx_int);
        std::vector<double> wgt(Nx_int);

        for (int i = 0; i < Nx_int; ++i)
        {
            s_int[i] = Lk * static_cast<double>(i) /
                           static_cast<double>(Nx_int - 1);

            wgt[i] = Lk - s_int[i];
        }

        const double P1 = My / Lk;
        const double Ps = Md / Lk;

        auto delta_of_P = [&](double P) -> double
        {
            std::vector<double> mvec(Nx_int);

            for (int i = 0; i < Nx_int; ++i)
            {
                mvec[i] =
                    std::min((P / P1) * (1.0 - s_int[i] / Lk), Md / My);
            }

            const std::vector<double> Kspan =
                disp_detail::local_Kspan_3region(
                    mvec,
                    Ky,
                    Kp,
                    Kd,
                    My,
                    Mp,
                    Ep,
                    t,
                    A_r2,
                    B_r2,
                    C_r2
                );

            std::vector<double> integrand(Nx_int);

            for (int i = 0; i < Nx_int; ++i)
            {
                integrand[i] = wgt[i] * Kspan[i];
            }

            return disp_detail::trapz(s_int, integrand);
        };

        if (xk <= delta_e_k)
        {
            Clamped_Moment[kk] = 3.0 * E * I * xk / (Lk * Lk);
            Clamped_Curvature[kk] = Clamped_Moment[kk] / (E * I);
        }
        else
        {
            if (xk >= delta_sat_k)
            {
                Clamped_Moment[kk] = Md;
                Clamped_Curvature[kk] = Kd;
            }
            else
            {
                auto Froot = [&](double P) -> double
                {
                    return delta_of_P(P) - xk;
                };

                double P_sol = std::numeric_limits<double>::quiet_NaN();

                try
                {
                    const double f0 = Froot(0.0);
                    const double f1 = Froot(Ps);

                    if (disp_detail::is_finite(f0)
                        && disp_detail::is_finite(f1)
                        && (f0 == 0.0 || f1 == 0.0 || f0 * f1 <= 0.0))
                    {
                        P_sol =
                            disp_detail::bracketed_bisection(
                                Froot,
                                0.0,
                                Ps,
                                120,
                                1e-10
                            );
                    }
                }
                catch (...)
                {
                }

                if (!disp_detail::is_finite(P_sol))
                {
                    const int Ngrid = 3000;
                    double best_err = std::numeric_limits<double>::infinity();

                    for (int i = 0; i < Ngrid; ++i)
                    {
                        const double Pg =
                            Ps * static_cast<double>(i) /
                            static_cast<double>(Ngrid - 1);

                        const double dg = delta_of_P(Pg);
                        const double err = std::abs(dg - xk);

                        if (disp_detail::is_finite(err) && err < best_err)
                        {
                            best_err = err;
                            P_sol = Pg;
                        }
                    }
                }

                P_sol = disp_detail::clamp_val(P_sol, 0.0, Ps);

                Clamped_Moment[kk] = Lk * P_sol;

                const double m0 = Clamped_Moment[kk] / My;

                double K0 =
                    disp_detail::curvature_from_m_3region(
                        m0,
                        Ky,
                        Kp,
                        Kd,
                        My,
                        Mp,
                        Ep,
                        t,
                        A_r2,
                        B_r2,
                        C_r2
                    );

                if (K0 >= Kd || Clamped_Moment[kk] >= Md)
                {
                    K0 = Kd;
                    Clamped_Moment[kk] = Md;
                }

                Clamped_Curvature[kk] = K0;
            }
        }

        Springback[kk] =
            disp_detail::rad2deg(
                (6.0 * Clamped_Moment[kk] * Lk)
                / (E_prime * std::pow(t, 3))
            );

        theta_2[kk] = theta_1[kk] - Springback[kk];

        const double theta1_rad = disp_detail::deg2rad(theta_1[kk]);

        const double denom =
            Ld
            * (
                (1.0 - (Rt / Ld) * std::sin(theta1_rad))
                + (
                    (std::tan(theta1_rad) - f)
                    / (1.0 + f * std::tan(theta1_rad))
                  )
                  * (
                      (xk / Ld)
                      - (Rt / Ld) * (1.0 - std::cos(theta1_rad))
                    )
              );

        if (std::abs(denom) < 1e-12)
        {
            Force[kk] = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            Force[kk] = (Clamped_Moment[kk] * w) / denom;
        }
    }

    // ----------------------------------------------------
    // New output:
    // Integral controller gain
    //
    // System gain:
    //      G_system = delta_theta2 / delta_x
    //
    // Integral controller gain:
    //      K_integral = 1 / G_system = delta_x / delta_theta2
    //
    // The array size is the same as x_range.
    // Endpoint points use forward/backward difference.
    // Interior points use central difference.
    // ----------------------------------------------------
    std::vector<double> integral_controller_gain =
        disp_detail::calculate_integral_controller_gain(
            x_range,
            theta_2
        );

    // ----------------------------------------------------
    // Invert desired unloaded angle theta_2 to x_required
    // ----------------------------------------------------
    const double x_required =
        disp_detail::invert_y_to_x_linear(
            x_range,
            theta_2,
            Desired_angle
        );

    const double theta_1_at_x =
        disp_detail::linear_interp(
            x_required,
            x_range,
            theta_1
        );

    const double springback_at_x =
        disp_detail::linear_interp(
            x_required,
            x_range,
            Springback
        );

    out.x_required = x_required;
    out.theta_1_at_x = theta_1_at_x;
    out.springback_at_x = springback_at_x;
    out.integral_controller_gain = integral_controller_gain;

    model_data.x = x_range;
    model_data.theta_1 = theta_1;
    model_data.theta_2 = theta_2;
    model_data.Springback = Springback;
    model_data.Clamped_Moment = Clamped_Moment;
    model_data.Clamped_Curvature = Clamped_Curvature;
    model_data.Force = Force;
    model_data.S1 = S1;
    model_data.integral_controller_gain = integral_controller_gain;

    model_data.delta_elastic_final = delta_elastic_final;
    model_data.delta_R2_final = delta_R2_final;
    model_data.My = My;
    model_data.Ky = Ky;
    model_data.Md = Md;
    model_data.Kd = Kd;

    return out;
}

}