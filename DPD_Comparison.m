% DPD_Comparison.m — ILC-DPD vs B-spline DPD vs GMP DPD
%
% ==========================================================================
% Three DPD methods, same memoryless PA, same waveform, head-to-head.
%
% Method 1: ILC-DPD (Iterative Learning Control DPD)
%   - No model, per-sample complex ratio correction
%   - I/Q averaging for noise robustness
%
% Method 2: B-spline DPD
%   - Degree-3 cubic B-splines, 11 break points, K=13 coefficients
%   - WLS additive update with amplitude weighting
%
% Method 3: GMP DPD (memoryless polynomial)
%   - Odd-order polynomial: orders [1, 3, 5, 7, 9], K=5 coefficients
%   - WLS additive update with same weighting
%
% All trained with AWGN noise (SNR=45 dB) for realistic comparison.
% ==========================================================================

clear; clc; close all;
rng(42);

%% 1. Load waveform
load('RefX.mat');
x = RefX;
r = abs(x);
N = length(x);
fprintf('Signal: N=%d  max|x|=%.4f\n', N, max(r));

%% 2. Memoryless PA model (Rapp AM/AM + AM/PM)
G_pa   = 1.0;
A_sat  = 0.85;    p_rapp = 2.5;
phi1   = 22;      phi2   = 18;
r1     = 0.70;    r2     = 0.82;

amam     = @(rv) G_pa*rv ./ (1 + (rv/A_sat).^(2*p_rapp)).^(1/(2*p_rapp));
ampm_rad = @(rv) (phi1*rv.^2./(r1^2+rv.^2) + phi2*rv.^6./(r2^6+rv.^6)) * (pi/180);
pa_apply = @(xv) amam(abs(xv)) ./ max(abs(xv),1e-12) .* xv .* exp(1j*ampm_rad(abs(xv)));

%% 3. AWGN noise
SNR_dB    = 30;
sig_pow   = mean(abs(x).^2);
noise_std = sqrt(sig_pow / 10^(SNR_dB/10) / 2);
add_noise = @(y) y + noise_std * (randn(size(y)) + 1j*randn(size(y)));

%% 4. Psat probe
G_target = 1.0;
r_cap    = 1.8;
r_fine   = linspace(0, r_cap, 100000)';
x_probe  = (r_cap ./ max(r, 1e-12)) .* x;

Psat       = max(abs(pa_apply(x_probe))) * 0.99;
[~, idx]   = min(abs(amam(r_fine) - Psat));
r_drive_max = r_fine(idx);
fprintf('Psat=%.4f  r_drive_max=%.4f\n', Psat, r_drive_max);

% Baseline: PA without DPD
Y_pa_only  = pa_apply(x);
mse_no_dpd = mean(abs(Y_pa_only - G_target*x).^2);
fprintf('No DPD MSE = %.3e\n\n', mse_no_dpd);

%% Common parameters
n_iter  = 20;
w_alpha = 20;
w_pow   = 4;
lam     = 1e-7;

% =====================================================================
%% 5. Method 1: ILC-DPD
% =====================================================================
fprintf('=== Method 1: ILC-DPD ===\n');
mu    = 0.8;
N_avg = 8;

P_iter     = x;
mse_iter   = zeros(1, n_iter+1);
mse_iter(1) = mse_no_dpd;
Psat_it    = Psat;

for k = 1:n_iter
    % I/Q averaged capture
    M_avg = zeros(N, 1);
    for j = 1:N_avg
        M_avg = M_avg + add_noise(pa_apply(P_iter));
    end
    M_avg = M_avg / N_avg;
    Psat_it = max(Psat_it, max(abs(M_avg)) * 0.99);

    % Target
    r_tgt = min(G_target * r, Psat_it);
    T     = r_tgt ./ max(r, 1e-12) .* x;

    % Complex ratio correction
    corr   = T ./ max(abs(M_avg), 1e-12) .* exp(-1j*angle(M_avg));
    P_iter = P_iter .* ((1 - mu) + mu * corr);

    % Clamp
    P_iter = P_iter .* min(1, r_drive_max ./ max(abs(P_iter), 1e-12));

    Y_chk = pa_apply(P_iter);
    mse_iter(k+1) = mean(abs(Y_chk - G_target*x).^2);
    fprintf('  Iter %2d  MSE = %.3e\n', k, mse_iter(k+1));
end
Y_dpd_iter = pa_apply(P_iter);
fprintf('  Final: %.3e  (%.1f dB)\n\n', mse_iter(end), 10*log10(mse_no_dpd/mse_iter(end)));

% =====================================================================
%% 6. Method 2: B-spline DPD
% =====================================================================
fprintf('=== Method 2: B-spline DPD (deg=3, K=13) ===\n');

deg   = 3;
t_brk = (0:0.1:1)';
t_kv  = [repmat(t_brk(1), deg, 1); t_brk; repmat(t_brk(end), deg, 1)];
K_bsp = length(t_kv) - deg - 1;   % = 13

B_real   = buildB(r, t_kv, deg, K_bsp);
Phi_bsp  = B_real .* x;
w_vec    = 1 + w_alpha * r.^w_pow;
sw       = sqrt(w_vec);
Phi_sw   = Phi_bsp .* sw;
AtWA_bsp = Phi_sw' * Phi_sw + lam * eye(K_bsp);

% Identity init
c_bsp      = AtWA_bsp \ (Phi_sw' * (x .* sw));
mse_bsp    = zeros(1, n_iter+1);
mse_bsp(1) = mse_no_dpd;
Psat_bs    = Psat;

for k = 1:n_iter
    xd      = Phi_bsp * c_bsp;
    xd      = xd .* min(1, r_drive_max ./ max(abs(xd), 1e-12));
    y       = add_noise(pa_apply(xd));
    Psat_bs = max(Psat_bs, max(abs(y)) * 0.99);
    r_tgt   = min(G_target * r, Psat_bs);
    x_tgt   = r_tgt ./ max(r, 1e-12) .* x;
    dc      = AtWA_bsp \ (Phi_sw' * ((x_tgt - y) .* sw));
    c_bsp   = c_bsp + dc;

    xd_chk = Phi_bsp * c_bsp;
    xd_chk = xd_chk .* min(1, r_drive_max ./ max(abs(xd_chk), 1e-12));
    Y_chk  = pa_apply(xd_chk);
    mse_bsp(k+1) = mean(abs(Y_chk - G_target*x).^2);
    fprintf('  Iter %2d  MSE = %.3e\n', k, mse_bsp(k+1));
end
xd_bsp   = Phi_bsp * c_bsp;
xd_bsp   = xd_bsp .* min(1, r_drive_max ./ max(abs(xd_bsp), 1e-12));
Y_dpd_bsp = pa_apply(xd_bsp);
fprintf('  Final: %.3e  (%.1f dB)\n\n', mse_bsp(end), 10*log10(mse_no_dpd/mse_bsp(end)));

% =====================================================================
%% 7. Method 3: GMP DPD (memoryless polynomial)
% =====================================================================
fprintf('=== Method 3: GMP DPD (orders [1,3,5,7,9], K=5) ===\n');

gmp_orders = [1 3 5 7 9];
K_gmp      = length(gmp_orders);

Phi_gmp = zeros(N, K_gmp);
for kk = 1:K_gmp
    Phi_gmp(:, kk) = x .* r.^(gmp_orders(kk) - 1);
end
Phi_gmp_sw = Phi_gmp .* sw;
AtWA_gmp   = Phi_gmp_sw' * Phi_gmp_sw + lam * eye(K_gmp);

% Identity init
c_gmp      = AtWA_gmp \ (Phi_gmp_sw' * (x .* sw));
mse_gmp    = zeros(1, n_iter+1);
mse_gmp(1) = mse_no_dpd;
Psat_gm    = Psat;

for k = 1:n_iter
    xd      = Phi_gmp * c_gmp;
    xd      = xd .* min(1, r_drive_max ./ max(abs(xd), 1e-12));
    y       = add_noise(pa_apply(xd));
    Psat_gm = max(Psat_gm, max(abs(y)) * 0.99);
    r_tgt   = min(G_target * r, Psat_gm);
    x_tgt   = r_tgt ./ max(r, 1e-12) .* x;
    dc      = AtWA_gmp \ (Phi_gmp_sw' * ((x_tgt - y) .* sw));
    c_gmp   = c_gmp + dc;

    xd_chk = Phi_gmp * c_gmp;
    xd_chk = xd_chk .* min(1, r_drive_max ./ max(abs(xd_chk), 1e-12));
    Y_chk  = pa_apply(xd_chk);
    mse_gmp(k+1) = mean(abs(Y_chk - G_target*x).^2);
    fprintf('  Iter %2d  MSE = %.3e\n', k, mse_gmp(k+1));
end
xd_gmp    = Phi_gmp * c_gmp;
xd_gmp    = xd_gmp .* min(1, r_drive_max ./ max(abs(xd_gmp), 1e-12));
Y_dpd_gmp = pa_apply(xd_gmp);
fprintf('  Final: %.3e  (%.1f dB)\n\n', mse_gmp(end), 10*log10(mse_no_dpd/mse_gmp(end)));

% =====================================================================
%% 8. Summary Table
% =====================================================================
fprintf('================================================================\n');
fprintf('                    DPD Comparison Summary\n');
fprintf('================================================================\n');
fprintf('PA: Memoryless Rapp (A_sat=%.2f, p=%.1f, phi1=%d, phi2=%d)\n', A_sat, p_rapp, phi1, phi2);
fprintf('Training: %d iters, SNR=%d dB\n\n', n_iter, SNR_dB);
fprintf('  %-25s  %-12s  %-12s  %-10s\n', 'Method', 'MSE', 'Improvement', 'Coefficients');
fprintf('  %-25s  %-12s  %-12s  %-10s\n', '------', '---', '-----------', '------------');
fprintf('  %-25s  %.3e     ---         ---\n', 'No DPD', mse_no_dpd);
fprintf('  %-25s  %.3e     %.1f dB      %d (samples)\n', ...
    'ILC-DPD', mse_iter(end), 10*log10(mse_no_dpd/mse_iter(end)), N);
fprintf('  %-25s  %.3e     %.1f dB      %d\n', ...
    'B-spline (deg=3, K=13)', mse_bsp(end), 10*log10(mse_no_dpd/mse_bsp(end)), K_bsp);
fprintf('  %-25s  %.3e     %.1f dB      %d\n', ...
    'GMP (orders 1-9, K=5)', mse_gmp(end), 10*log10(mse_no_dpd/mse_gmp(end)), K_gmp);
fprintf('================================================================\n');

% =====================================================================
%% 9. Visualization
% =====================================================================

r_ax = linspace(0, 1, 500)';
hardclip_ax = min(G_target * r_ax, Psat);
[r_sorted, si] = sort(r);

% Colors (solid, clearly distinguishable)
c_iter = [0.85 0.1 0.1];   % red
c_bsp  = [0.0  0.5 0.0];   % green
c_gmp  = [0.1  0.1 0.85];  % blue
mk_a   = 0.40;              % scatter alpha

scr = get(0, 'ScreenSize');
fig = figure('Position', scr, 'Name', 'DPD Comparison');
set(fig, 'WindowState', 'maximized');

% ---- (1) DPD+PA AM/AM ----
subplot(2,3,1);
scatter(r_sorted, abs(Y_dpd_iter(si)), 2, c_iter, 'filled', 'MarkerFaceAlpha', mk_a); hold on;
scatter(r_sorted, abs(Y_dpd_bsp(si)),  2, c_bsp,  'filled', 'MarkerFaceAlpha', mk_a);
scatter(r_sorted, abs(Y_dpd_gmp(si)),  2, c_gmp,  'filled', 'MarkerFaceAlpha', mk_a);
plot(r_ax, G_target*r_ax, 'k--', 'LineWidth', 1.5);
plot(r_ax, hardclip_ax,   'm:',  'LineWidth', 2.0);
legend('ILC-DPD','B-spline','GMP','Ideal','Hard-clip','Location','northwest');
xlabel('|x|'); ylabel('|Y|');
title('DPD+PA AM/AM'); grid on;

% ---- (2) DPD+PA AM/PM ----
ph_iter = angle(Y_dpd_iter(si) .* conj(x(si))) * 180/pi;
ph_bsp  = angle(Y_dpd_bsp(si)  .* conj(x(si))) * 180/pi;
ph_gmp  = angle(Y_dpd_gmp(si)  .* conj(x(si))) * 180/pi;
ph_iter(r_sorted < 0.02) = NaN;
ph_bsp(r_sorted  < 0.02) = NaN;
ph_gmp(r_sorted  < 0.02) = NaN;

subplot(2,3,2);
scatter(r_sorted, ph_iter, 2, c_iter, 'filled', 'MarkerFaceAlpha', mk_a); hold on;
scatter(r_sorted, ph_bsp,  2, c_bsp,  'filled', 'MarkerFaceAlpha', mk_a);
scatter(r_sorted, ph_gmp,  2, c_gmp,  'filled', 'MarkerFaceAlpha', mk_a);
yline(0, 'k--', 'LineWidth', 1.5);
legend('ILC-DPD','B-spline','GMP','Ideal (0 deg)');
xlabel('|x|'); ylabel('Phase (deg)');
title('DPD+PA AM/PM (residual)'); grid on;

% ---- (3) Convergence ----
subplot(2,3,3);
semilogy(0:n_iter, mse_iter, 'o-', 'Color', c_iter, 'MarkerFaceColor', c_iter, ...
         'MarkerSize', 4, 'LineWidth', 1.5); hold on;
semilogy(0:n_iter, mse_bsp,  's-', 'Color', c_bsp,  'MarkerFaceColor', c_bsp, ...
         'MarkerSize', 4, 'LineWidth', 1.5);
semilogy(0:n_iter, mse_gmp,  'd-', 'Color', c_gmp,  'MarkerFaceColor', c_gmp, ...
         'MarkerSize', 4, 'LineWidth', 1.5);
yline(mse_no_dpd, 'k--', 'LineWidth', 1.0);
legend('ILC-DPD','B-spline','GMP','No DPD','Location','northeast');
xlabel('Iteration'); ylabel('MSE');
title('Convergence'); grid on;

% ---- (4) DPD gain |G_DPD| vs |x| ----
g_iter = abs(P_iter(si))   ./ max(r_sorted, 1e-12);
g_bsp  = abs(xd_bsp(si))   ./ max(r_sorted, 1e-12);
g_gmp  = abs(xd_gmp(si))   ./ max(r_sorted, 1e-12);
g_iter(r_sorted < 0.02) = NaN;
g_bsp(r_sorted  < 0.02) = NaN;
g_gmp(r_sorted  < 0.02) = NaN;

subplot(2,3,4);
scatter(r_sorted, g_iter, 2, c_iter, 'filled', 'MarkerFaceAlpha', mk_a); hold on;
scatter(r_sorted, g_bsp,  2, c_bsp,  'filled', 'MarkerFaceAlpha', mk_a);
scatter(r_sorted, g_gmp,  2, c_gmp,  'filled', 'MarkerFaceAlpha', mk_a);
yline(1, 'k--', 'LineWidth', 1.0);
legend('ILC-DPD','B-spline','GMP','Unity gain');
xlabel('|x|'); ylabel('|G_{DPD}|');
title('DPD Gain vs Amplitude'); grid on;

% ---- (5) DPD phase vs |x| ----
dp_iter = angle(P_iter(si)  .* conj(x(si))) * 180/pi;
dp_bsp  = angle(xd_bsp(si)  .* conj(x(si))) * 180/pi;
dp_gmp  = angle(xd_gmp(si)  .* conj(x(si))) * 180/pi;
dp_iter(r_sorted < 0.02) = NaN;
dp_bsp(r_sorted  < 0.02) = NaN;
dp_gmp(r_sorted  < 0.02) = NaN;

subplot(2,3,5);
scatter(r_sorted, dp_iter, 2, c_iter, 'filled', 'MarkerFaceAlpha', mk_a); hold on;
scatter(r_sorted, dp_bsp,  2, c_bsp,  'filled', 'MarkerFaceAlpha', mk_a);
scatter(r_sorted, dp_gmp,  2, c_gmp,  'filled', 'MarkerFaceAlpha', mk_a);
plot(r_ax, -ampm_rad(r_ax)*180/pi, 'k--', 'LineWidth', 1.5);
legend('ILC-DPD','B-spline','GMP','Ideal: -PA AM/PM');
xlabel('|x|'); ylabel('Phase (deg)');
title('DPD Phase vs Amplitude'); grid on;

% ---- (6) Summary bar chart ----
subplot(2,3,6);
mse_vals = [mse_iter(end), mse_bsp(end), mse_gmp(end)];
db_vals  = 10*log10(mse_no_dpd ./ mse_vals);
b = bar(db_vals, 0.6);
b.FaceColor = 'flat';
b.CData(1,:) = c_iter;
b.CData(2,:) = c_bsp;
b.CData(3,:) = c_gmp;
set(gca, 'XTickLabel', {'ILC-DPD','B-spline','GMP'});
ylabel('Improvement (dB)');
title('Final MSE Improvement');
grid on;
for i = 1:3
    text(i, db_vals(i)+0.3, sprintf('%.1f dB\n(%.1e)', db_vals(i), mse_vals(i)), ...
         'HorizontalAlignment','center', 'FontSize', 9, 'FontWeight', 'bold');
end

sgtitle(sprintf(['DPD Comparison  |  Memoryless PA  |  ' ...
    '%d iters, SNR=%d dB  |  No-DPD MSE=%.2e'], n_iter, SNR_dB, mse_no_dpd));

saveas(gcf, 'DPD_Comparison.png');
fprintf('\nSaved -> DPD_Comparison.png\n');

%% ========================================================================
%  Local Functions
%% ========================================================================

function B = buildB(x, t_kv, deg, n_bas)
    x = x(:);
    B = zeros(numel(x), n_bas);
    for i = 1:n_bas
        B(:,i) = bspF(x, t_kv, i, deg);
    end
end

function b = bspF(x, t, i, d)
    if d == 0
        if t(i) == t(i+1)
            b = zeros(size(x));
        else
            b = double(x >= t(i) & x < t(i+1));
            if i == numel(t)-1
                b(x == t(end)) = 1;
            end
        end
        return
    end
    b  = zeros(size(x));
    d1 = t(i+d)   - t(i);
    d2 = t(i+d+1) - t(i+1);
    if d1 ~= 0
        b = b + (x - t(i))     / d1 .* bspF(x, t, i,   d-1);
    end
    if d2 ~= 0
        b = b + (t(i+d+1) - x) / d2 .* bspF(x, t, i+1, d-1);
    end
end
