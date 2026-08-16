clc; close all; clear;

curvature_factor = 1.0;         % Compensates for local bending
E = 69000;                      % Elastic modulus (MPa)
nu = 0.0;
sigma_0 = 300;                  % Yield Stress (MPa)
Ep =950;                        % Plastic modulus (MPa)

Rd = 6.35;                      % mm die radius
Rt = 25.4;                      % mm tool radius
t  = 1.6;                       % mm plate thickness
w  = 25.4*2;                    % mm width
d_base  = 52;                   % mm die surface to tool center
S = Rd + Rt + t;
alpha = 30;                     % Semi arc angle (degree)
alpha_rad = deg2rad(alpha);

d_S1    = d_base;
d_theta = d_base * curvature_factor;

Delta_y_S1    = Rd + d_S1;
Delta_y_theta = Rd + d_theta;

nx = 500;                       % Number of displacement points
x  = linspace(0, 50, nx)';      % Punch displacement range

E_prime = E/(1 - nu^2);
I  = t^3/12;

Ky = (2*sigma_0)/(E*t);
My = (sigma_0 * t^2) / 6;

Kp = 2.5 * Ky;
Mp = 1.5 * My;

Kd = (1/(Rd + t/2));
Md = (Ep*t^3/12)*(Kd - Kp) + Mp;

theta_1 = nan(nx,1);
theta_2 = nan(nx,1);
theta1_only = nan(nx,1);
Springback = zeros(nx,1);

S1 = nan(nx,1);

Clamped_Moment = zeros(nx,1);
Clamped_Curvature = zeros(nx,1);
Trial_Force = zeros(nx,1);

delta_elastic_local = nan(nx,1);
delta_R2_local      = nan(nx,1);
delta_plastic_local = nan(nx,1);


theta1_phase1_geom_theta = @(xx) local_theta1_phase1(xx, S, Delta_y_theta);

theta1_x0_theta   = theta1_phase1_geom_theta(x(1));
theta1_xend_theta = theta1_phase1_geom_theta(x(end));

if isfinite(theta1_x0_theta) && isfinite(theta1_xend_theta) && ...
        theta1_x0_theta <= alpha && theta1_xend_theta >= alpha
    try
        x_alpha_common = fzero(@(xx) local_theta1_phase1(xx, S, Delta_y_theta) - alpha, [x(1), x(end)]);
    catch
        [~, idxa] = min(abs(arrayfun(theta1_phase1_geom_theta, x) - alpha));
        x_alpha_common = x(idxa);
    end
else
    [~, idxa] = min(abs(arrayfun(theta1_phase1_geom_theta, x) - alpha));
    x_alpha_common = x(idxa);
end

S1_phase1_at_xa = local_S1_phase1(x_alpha_common, S, Delta_y_S1, Rd, t);
S1_phase2_at_xa = local_S1_phase2(x_alpha_common, S, Delta_y_S1, Rt, Rd, t, alpha_rad);

S1_phase2_offset = S1_phase1_at_xa - S1_phase2_at_xa;

fprintf('\n================ Phase Calculations =================\n');
fprintf('d                                  = %.6f mm\n', d_S1);
fprintf('Phase shift                        = %.6f mm\n', x_alpha_common);
fprintf('Bending angle at phase shift       = %.6f deg\n', local_theta1_phase1(x_alpha_common, S, Delta_y_theta));
fprintf('Arc length at Phase-1 end          = %.6f mm\n', S1_phase1_at_xa);
fprintf('Arc length at Phase-2 start        = %.6f mm\n', S1_phase2_at_xa);

for k = 1:nx
    if x(k) <= x_alpha_common

        D = sqrt((S - x(k))^2 + Delta_y_S1^2);
        delta  = atan2(Delta_y_S1, max(S - x(k), 1e-12));
        lambda = acos(max(min(S / D, 1), -1));

        straight_len = sqrt(max(D^2 - S^2, 0));
        wrap_len = (Rd + t/2) * (delta - lambda);

        S1(k) = straight_len + wrap_len;

    else

        theta_v_S1 = solve_theta_v_phase2(x(k), S, Delta_y_S1, Rt, t, alpha_rad);

        denom_sin = sin(theta_v_S1);
        if abs(denom_sin) < 1e-9
            denom_sin = sign(denom_sin + eps) * 1e-9;
        end

        Rn = (Delta_y_S1 - Rt*sin(alpha_rad) - (t/2)*cos(theta_v_S1)) / denom_sin;

        straight_len = sqrt(max(Rn^2 - (Rd + t/2)^2, 0));

        theta_w = asin(max(min((Rd + t/2)/Rn, 1), -1));
        theta1_tmp = (pi/2) - theta_v_S1 + theta_w;
        wrap_len = (Rd + t/2) * theta1_tmp;
        S1(k) = straight_len + wrap_len + S1_phase2_offset;
    end

end

for k = 1:nx
    if x(k) <= x_alpha_common
        theta1_only(k) = local_theta1_phase1(x(k), S, Delta_y_theta);
    else
        theta1_only(k) = local_theta1_phase2(x(k), S, Delta_y_theta, Rt, Rd, t, alpha_rad);
    end
end

theta_1 = theta1_only;

valid = isfinite(S1) & isfinite(theta_1) & (S1 > 1e-9) & (theta_1 >= 0);
k_max = find(valid, 1, 'last');
if isempty(k_max)
    error('No valid theta_1/S1 computed. Check geometry and phase transition.');
end

k_ref = find(valid, 1, 'first');
S1_ref = S1(k_ref);

[A_r2, B_r2, C_r2] = solveABC_region2(Mp/My, Kp/Ky);

for k = 1:k_max
    Lk = S1(k);

    delta_elastic_local(k) = My * Lk^2 / (3 * E * I);

    [delta_R2_k, delta_sat_k] = get_transition_displacements_local( ...
        Lk, Ky, My, Kp, Mp, Kd, Md, Ep, t, A_r2, B_r2, C_r2);

    delta_R2_local(k) = delta_R2_k;
    delta_plastic_local(k) = delta_sat_k;
end

delta_elastic_final = delta_elastic_local(k_max);
delta_R2_final      = delta_R2_local(k_max);
delta_plastic_final = delta_plastic_local(k_max);

Nx_int = 300;

for k = 1:k_max
    if ~valid(k), continue; end

    Lk = S1(k);

    delta_e_k   = delta_elastic_local(k);
    delta_R2_k  = delta_R2_local(k);
    delta_sat_k = delta_plastic_local(k);

    s_int = linspace(0, Lk, Nx_int);
    wgt   = (Lk - s_int);

    P1 = My / Lk;
    Ps = Md / Lk;

    delta_of_P = @(P) trapz(s_int, wgt .* local_Kspan_3region( ...
        min((P/P1) .* (1 - s_int./Lk), Md/My), ...
        Ky, Kp, Kd, My, Mp, Md, Ep, t, A_r2, B_r2, C_r2));

    if x(k) <= delta_e_k
        Clamped_Moment(k) = 3 * E * I * x(k) / (Lk^2);
        Clamped_Curvature(k) = Clamped_Moment(k) / (E * I);
    else
        if x(k) >= delta_sat_k
            Clamped_Moment(k) = Md;
            Clamped_Curvature(k) = Kd;
        else
            Froot = @(P) delta_of_P(P) - x(k);

            try
                P_sol = fzero(Froot, [0, Ps]);
            catch
                P_grid = linspace(0, Ps, 3000);
                d_grid = arrayfun(delta_of_P, P_grid);
                [~, idxMin] = min(abs(d_grid - x(k)));
                P_sol = P_grid(idxMin);
            end

            P_sol = max(0, min(P_sol, Ps));
            Clamped_Moment(k) = Lk * P_sol;

            m0 = Clamped_Moment(k) / My;
            K0 = curvature_from_m_3region(m0, Ky, Kp, Kd, My, Mp, Md, Ep, t, A_r2, B_r2, C_r2);

            if K0 >= Kd || Clamped_Moment(k) >= Md
                K0 = Kd;
                Clamped_Moment(k) = Md;
            end

            Clamped_Curvature(k) = K0;
        end
    end

    Springback(k) = rad2deg((6 * Clamped_Moment(k) * Lk) / (E_prime * t^3));
    theta_2(k)    = theta_1(k) - Springback(k);

    if abs(Lk) < 1e-12
        Trial_Force(k) = NaN;
    else
        Trial_Force(k) = (Clamped_Moment(k) / Lk) * w;
    end

end

idx_elastic_final = find(Clamped_Moment(1:k_max) >= My, 1, 'first');
idx_R2_curve      = find(Clamped_Moment(1:k_max) >= Mp, 1, 'first');
idx_plastic_final = find(Clamped_Moment(1:k_max) >= Md - 1e-9, 1, 'first');


fprintf('Elastic limit                      = %.6f mm\n', delta_elastic_final);
fprintf('Elastoplastic limit                = %.6f mm\n', delta_R2_final);

theta2_for_gain = theta_2(1:k_max);
x_for_gain = x(1:k_max);

valid_gain = isfinite(theta2_for_gain) & isfinite(x_for_gain);
theta2_for_gain = theta2_for_gain(valid_gain);
x_for_gain = x_for_gain(valid_gain);

system_gain = nan(size(theta2_for_gain));
controller_gain = nan(size(theta2_for_gain));

if numel(x_for_gain) >= 2
    system_gain = gradient(theta2_for_gain, x_for_gain);

    valid_controller_gain = isfinite(system_gain) & abs(system_gain) > 1e-12;
    controller_gain(valid_controller_gain) = 1 ./ system_gain(valid_controller_gain);
end

pStyle.Width       = 3.5; % Figure width (in)
pStyle.Height      = 2.5; % Figure height (in)
pStyle.FontName    = 'Arial';
pStyle.FontSize    = 10; % pt
pStyle.FontWeight  = 'bold';
pStyle.LineWidth   = 1.5;
pStyle.MarkerSize  = 4;
pStyle.BorderWidth = 1;

defl = x(1:k_max);

fig = figure('Color','w'); box on; hold on;
fig.Units = 'inches'; fig.Position = [1, 1, pStyle.Width, pStyle.Height];
plot(defl, S1(1:k_max), 'LineWidth', pStyle.LineWidth, 'MarkerSize', pStyle.MarkerSize, 'DisplayName', 'arc length');
xline(x_alpha_common, 'm', 'LineWidth', pStyle.LineWidth, 'DisplayName', 'phase change');
xlabel('punch displacement (mm)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
ylabel('arc length (mm)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
xlim([0 max(defl)]);
ylim([40 70]);
ax = gca;
ax.FontName = pStyle.FontName; ax.FontSize = pStyle.FontSize; ax.FontWeight = pStyle.FontWeight; ax.LineWidth = pStyle.BorderWidth;
set(findall(fig, 'Type', 'text'), 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
lgd = legend('Location','best','box','off'); lgd.FontName = pStyle.FontName; lgd.FontSize = pStyle.FontSize; lgd.FontWeight = pStyle.FontWeight;
grid off;
hold off;

fig = figure('Color','w'); box on; hold on;
fig.Units = 'inches'; fig.Position = [1, 1, pStyle.Width, pStyle.Height];
plot(defl, theta_1(1:k_max), '-', 'LineWidth', pStyle.LineWidth, 'MarkerSize', pStyle.MarkerSize, 'DisplayName', 'bending angle');
xline(x_alpha_common, 'm', 'LineWidth', pStyle.LineWidth, 'DisplayName', 'phase change');
xlabel('punch displacement (mm)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
ylabel('bending angle (deg)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
xlim([0 max(defl)]);
ax = gca;
ax.FontName = pStyle.FontName; ax.FontSize = pStyle.FontSize; ax.FontWeight = pStyle.FontWeight; ax.LineWidth = pStyle.BorderWidth;
set(findall(fig, 'Type', 'text'), 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
lgd = legend('Location','best','box','off'); lgd.FontName = pStyle.FontName; lgd.FontSize = pStyle.FontSize; lgd.FontWeight = pStyle.FontWeight;
grid off;
hold off;

fig = figure('Color','w'); box on; hold on;
fig.Units = 'inches'; fig.Position = [1, 1, pStyle.Width, pStyle.Height];
plot(defl, theta_2(1:k_max), '--', 'LineWidth', pStyle.LineWidth, 'MarkerSize', pStyle.MarkerSize, 'DisplayName', 'unloaded angle');
xline(delta_elastic_final, ':k', 'LineWidth', pStyle.LineWidth, 'DisplayName', 'elastic limit');
xline(delta_R2_final, ':m', 'LineWidth', pStyle.LineWidth, 'DisplayName', 'elastoplastic limit');
xlabel('punch displacement (mm)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
ylabel('unloaded angle (deg)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
xlim([0 max(defl)]);
ylim([-5 45]);
ax = gca;
ax.FontName = pStyle.FontName; ax.FontSize = pStyle.FontSize; ax.FontWeight = pStyle.FontWeight; ax.LineWidth = pStyle.BorderWidth;
set(findall(fig, 'Type', 'text'), 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
lgd = legend('Location','best','box','off'); lgd.FontName = pStyle.FontName; lgd.FontSize = pStyle.FontSize; lgd.FontWeight = pStyle.FontWeight;
grid off;
hold off;

fig = figure('Color','w'); box on; hold on;
fig.Units = 'inches'; fig.Position = [1, 1, pStyle.Width, pStyle.Height];
plot(defl, Trial_Force(1:k_max), 'LineWidth', pStyle.LineWidth, 'MarkerSize', pStyle.MarkerSize);
xlabel('punch displacement (mm)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
ylabel('bending force (N)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
xlim([0 max(defl)]);
ax = gca;
ax.FontName = pStyle.FontName; ax.FontSize = pStyle.FontSize; ax.FontWeight = pStyle.FontWeight; ax.LineWidth = pStyle.BorderWidth;
set(findall(fig, 'Type', 'text'), 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
grid off;
hold off;

fig = figure('Color','w'); box on; hold on;
fig.Units = 'inches'; fig.Position = [1, 1, pStyle.Width, pStyle.Height];
plot(x_for_gain, system_gain, 'LineWidth', pStyle.LineWidth, 'MarkerSize', pStyle.MarkerSize, 'DisplayName', 'process gain');
xline(delta_elastic_final, ':k', 'LineWidth', pStyle.LineWidth, 'DisplayName', 'elastic limit');
xline(delta_R2_final, ':m', 'LineWidth', pStyle.LineWidth, 'DisplayName', 'elastoplastic limit');
xlabel('punch displacement (mm)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
ylabel('process gain (deg/mm)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
xlim([min(x_for_gain) max(x_for_gain)]);
ylim([min(system_gain) 1.1*max(system_gain)]);
ax = gca;
ax.FontName = pStyle.FontName; ax.FontSize = pStyle.FontSize; ax.FontWeight = pStyle.FontWeight; ax.LineWidth = pStyle.BorderWidth;
set(findall(fig, 'Type', 'text'), 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
lgd = legend('Location','southeast','box','off'); lgd.FontName = pStyle.FontName; lgd.FontSize = pStyle.FontSize; lgd.FontWeight = pStyle.FontWeight;
grid off;
hold off;

fig = figure('Color','w'); box on; hold on;
fig.Units = 'inches'; fig.Position = [1, 1, pStyle.Width, pStyle.Height];
plot(x_for_gain, controller_gain, 'LineWidth', pStyle.LineWidth, 'MarkerSize', pStyle.MarkerSize, 'DisplayName', 'system gain');
xline(delta_elastic_final, ':k', 'LineWidth', pStyle.LineWidth, 'DisplayName', 'elastic limit');
xline(delta_R2_final, ':m', 'LineWidth', pStyle.LineWidth, 'DisplayName', 'elastoplastic limit');
xlabel('punch displacement (mm)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
ylabel('system gain (mm/deg)', 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
xlim([min(x_for_gain) max(x_for_gain)]);
ylim([-10 10]);
ax = gca;
ax.FontName = pStyle.FontName; ax.FontSize = pStyle.FontSize; ax.FontWeight = pStyle.FontWeight; ax.LineWidth = pStyle.BorderWidth;
set(findall(fig, 'Type', 'text'), 'FontName', pStyle.FontName, 'FontSize', pStyle.FontSize, 'FontWeight', pStyle.FontWeight);
lgd = legend('Location','best','box','off'); lgd.FontName = pStyle.FontName; lgd.FontSize = pStyle.FontSize; lgd.FontWeight = pStyle.FontWeight;
grid off;
hold off;

function th1_deg = local_theta1_phase1(xx, S, Delta_y)
    D = sqrt((S - xx).^2 + Delta_y.^2);
    delta  = atan2(Delta_y, max(S - xx, 1e-12));
    lambda = acos(max(min(S ./ D, 1), -1));
    th1_deg = rad2deg(delta - lambda);
end

function S1_val = local_S1_phase1(xx, S, Delta_y, Rd, t)
    D = sqrt((S - xx).^2 + Delta_y.^2);
    delta  = atan2(Delta_y, max(S - xx, 1e-12));
    lambda = acos(max(min(S ./ D, 1), -1));
    straight_len = sqrt(max(D.^2 - S^2, 0));
    wrap_len = (Rd + t/2) * (delta - lambda);
    S1_val = straight_len + wrap_len;
end

function theta_v = solve_theta_v_phase2(xx, S, Delta_y, Rt, t, alpha_rad)
    fun = @(th) th - atan2(Delta_y - Rt*sin(alpha_rad) - (t/2)*cos(th), ...
                           -S + xx + Rt*cos(alpha_rad) + (t/2)*sin(th));

    guess_list = linspace(1e-4, pi/2 - 1e-4, 7);
    theta_v = NaN;
    best_res = inf;

    for j = 1:numel(guess_list)
        try
            th_try = fzero(fun, guess_list(j));
            res = abs(fun(th_try));
            if isfinite(th_try) && isreal(th_try) && res < best_res
                theta_v = th_try;
                best_res = res;
            end
        catch
        end
    end

    if ~isfinite(theta_v)
        theta_v = 0;
    end
end

function S1_val = local_S1_phase2(xx, S, Delta_y, Rt, Rd, t, alpha_rad)
    theta_v = solve_theta_v_phase2(xx, S, Delta_y, Rt, t, alpha_rad);

    denom_sin = sin(theta_v);
    if abs(denom_sin) < 1e-9
        denom_sin = sign(denom_sin + eps) * 1e-9;
    end

    Rn = (Delta_y - Rt*sin(alpha_rad) - (t/2)*cos(theta_v)) / denom_sin;

    straight_len = sqrt(max(Rn^2 - (Rd + t/2)^2, 0));

    theta_w = asin(max(min((Rd + t/2)/Rn, 1), -1));
    theta1_tmp = (pi/2) - theta_v + theta_w;
    wrap_len = (Rd + t/2) * theta1_tmp;

    S1_val = straight_len + wrap_len;
end

function th1_deg = local_theta1_phase2(xx, S, Delta_y, Rt, Rd, t, alpha_rad)
    theta_v = solve_theta_v_phase2(xx, S, Delta_y, Rt, t, alpha_rad);

    denom_sin = sin(theta_v);
    if abs(denom_sin) < 1e-9
        denom_sin = sign(denom_sin + eps) * 1e-9;
    end

    Rn = (Delta_y - Rt*sin(alpha_rad) - (t/2)*cos(theta_v)) / denom_sin;
    ratio_asin = (Rd + t/2) / Rn;
    ratio_asin = max(-1, min(1, ratio_asin));
    theta_w = asin(ratio_asin);

    th1_deg = rad2deg((pi/2) - theta_v + theta_w);
end

function [A,B,C] = solveABC_region2(ms, kh)
    g = @(u) 1./(u.*sqrt(1 - 2*u*(ms-1))) - (kh + 1./u - 1);

    u_max = 1/(2*(ms-1));
    u_lo  = 1e-10;
    g_lo  = g(u_lo);
    if ~isfinite(g_lo)
        u_lo = 1e-8;
        g_lo = g(u_lo);
    end

    found = false;
    u_hi = NaN;

    for p = 3:16
        u_try = u_max * (1 - 10^(-p));
        g_try = g(u_try);
        if isfinite(g_try) && sign(g_try) ~= sign(g_lo)
            u_hi = u_try;
            found = true;
            break
        end
    end

    if ~found
        us = u_max * (1 - logspace(-16, -3, 5000));
        gv = arrayfun(@(uu) localSafeG(g, uu), us);
        sgn = sign(gv);
        j = find(sgn(1:end-1).*sgn(2:end) < 0, 1, 'first');
        if isempty(j)
            error('Could not bracket u.');
        end
        u_lo = us(j);
        u_hi = us(j+1);
    end

    u = fzero(g, [u_lo, u_hi]);

    B = 2*u^3;
    C = 1/u - 1;
    A = u^2 + 2*u^3;
end

function val = localSafeG(g, u)
    try
        val = g(u);
        if ~isfinite(val)
            val = NaN;
        end
    catch
        val = NaN;
    end
end

function Kval = curvature_from_m_3region(m, Ky, K2, Kd, My, M2, ~, Ep, t, A, B, C)
    m = max(m, 0);
    ms = M2 / My;

    if m <= 1
        Kval = Ky * m;
    elseif m <= ms
        Kval = Ky * (1 ./ sqrt(A - B*m) - C);
    else
        M = m * My;
        Kval = K2 + (M - M2) / (Ep * t^3 / 12);
        if Kval > Kd
            Kval = Kd;
        end
    end
end

function Kspan = local_Kspan_3region(mvec, Ky, K2, Kd, My, M2, ~, Ep, t, A, B, C)
    Kspan = zeros(size(mvec));
    ms = M2 / My;

    idx1 = (mvec <= 1);
    idx2 = (mvec > 1) & (mvec <= ms);
    idx3 = (mvec > ms);

    Kspan(idx1) = Ky .* mvec(idx1);
    Kspan(idx2) = Ky .* (1 ./ sqrt(A - B.*mvec(idx2)) - C);

    Mvec = mvec(idx3) .* My;
    Kspan(idx3) = K2 + (Mvec - M2) ./ (Ep * t^3 / 12);
    Kspan(idx3) = min(Kspan(idx3), Kd);
end

function [delta_R2, delta_sat] = get_transition_displacements_local( ...
    L, Ky, My, K2, M2, Kd, Md, Ep, t, A, B, C)

    Ns = 3000;
    s = linspace(0, L, Ns);
    wgt = (L - s);

    P1 = My / L;
    P2 = M2 / L;
    Ps = Md / L;

    delta_of_P = @(P) trapz(s, wgt .* local_Kspan_3region( ...
        min((P/P1) .* (1 - s./L), Md/My), ...
        Ky, K2, Kd, My, M2, Md, Ep, t, A, B, C));

    delta_R2  = delta_of_P(P2);
    delta_sat = delta_of_P(Ps);
end
