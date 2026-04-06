% IterativeDPD.m — ILC-DPD (Iterative Learning Control Digital Predistortion)
%
% ==========================================================================
% Algorithm Overview
% ==========================================================================
%
% Model-free, sample-by-sample waveform predistortion.
% No basis functions, no WLS, no model coefficients.
%
% ILC-DPD (Chani-Cahuana et al., IEEE MTT 2016; Schoukens et al., 2017):
%   We want PA(p) = y_d, where y_d is the desired linear output.
%   Since PA^{-1} is unknown, we iterate (multiplicative form):
%
%     p_0 = x                                  (start with original signal)
%     For iteration k = 1,2,...:
%       y_k   = PA(p_k) + noise                (measure PA output)
%       y_avg = average over N_avg captures     (I/Q averaging to denoise)
%       y_d(n) = min(G*|x(n)|, Psat) * e^{j*angle(x(n))}  (target output)
%       p_{k+1}(n) = p_k(n) * y_d(n) / y_avg(n)           (complex ratio)
%       clip |p_{k+1}| <= r_drive_max           (prevent PA overdrive)
%
% ==========================================================================
% PA Model Selection
% ==========================================================================
%
% 'memoryless' — Rapp AM/AM + AM/PM (output depends only on current input)
%   PA(x(n)) = f(|x(n)|) * x(n)
%
% 'memory' — Wiener model: FIR filter → memoryless Rapp
%   x → [h0, h1, h2, h3] → x_filtered → Rapp AM/AM+AM/PM → y
%   Output y(n) depends on x(n), x(n-1), x(n-2), x(n-3)
%   The FIR filter models bandwidth limitation / electrical memory
%   before the nonlinear transistor (standard PA modeling approach).
%   Memory causes AM/AM to "spread" — same |x| gives different |y|
%   depending on signal history.
%
% ==========================================================================
% SNR Handling
% ==========================================================================
%
% I/Q averaging: average N_avg noisy captures per iteration.
% Reduces noise by sqrt(N_avg). Essential for sample-by-sample DPD
% (no model to smooth). Learning rate mu < 1 adds further stability.
%
% ==========================================================================

clear; clc; close all;
rng(42);

%% ========== USER CONFIG ==========
pa_mode = 'memory';     % <<< 'memoryless' or 'memory'
% ==================================

%% 1. Load waveform
load('RefX.mat');
x = RefX;
r = abs(x);
N = length(x);
fprintf('Signal: N=%d  max|x|=%.4f\n', N, max(r));

%% 2. PA models

% ---- Memoryless PA: Rapp AM/AM + AM/PM ----
G_pa   = 1.0;
A_sat  = 0.85;    p_rapp = 2.5;
phi1   = 22;      phi2   = 18;
r1     = 0.70;    r2     = 0.82;

amam     = @(rv) G_pa*rv ./ (1 + (rv/A_sat).^(2*p_rapp)).^(1/(2*p_rapp));
ampm_rad = @(rv) (phi1*rv.^2./(r1^2+rv.^2) + phi2*rv.^6./(r2^6+rv.^6)) * (pi/180);
pa_memoryless = @(xv) amam(abs(xv)) ./ max(abs(xv),1e-12) .* xv ...
                      .* exp(1j*ampm_rad(abs(xv)));

% ---- Memory PA: Wiener model (FIR → memoryless PA) ----
%
% FIR taps: h = [1.0, h1, h2, h3]
%   h1, h2, h3 are small complex coefficients (~5-10% of main tap)
%   They introduce inter-sample coupling (memory).
%   Complex phases model group delay / dispersion.
%
% Physical meaning:
%   h(0) = 1.0            : main signal path (dominant)
%   h(1) = 0.15*e^{j0.5}  : 1-sample delayed echo (~15%, electrical memory)
%   h(2) = -0.10*e^{-j0.3}: 2-sample echo (~10%, impedance mismatch reflection)
%   h(3) = 0.05*e^{j0.7}  : 3-sample echo (~5%, thermal memory tail)

h_mem = [1.0; 0.15*exp(1j*0.5); -0.10*exp(-1j*0.3); 0.05*exp(1j*0.7)];
pa_memory = @(xv) pa_memoryless(filter(h_mem, 1, xv));

% ---- Select PA model ----
switch pa_mode
    case 'memoryless'
        pa_model = pa_memoryless;
        n_iter   = 20;
        fprintf('PA mode: MEMORYLESS (Rapp AM/AM + AM/PM)\n');
    case 'memory'
        pa_model = pa_memory;
        n_iter   = 30;    % memory needs more iterations
        fprintf('PA mode: MEMORY (Wiener: FIR[%d taps] -> Rapp)\n', length(h_mem));
        fprintf('  FIR taps: ');
        for i = 1:length(h_mem)
            fprintf('h(%d)=%.3f%+.3fj  ', i-1, real(h_mem(i)), imag(h_mem(i)));
        end
        fprintf('\n');
    otherwise
        error('Unknown pa_mode: %s (use ''memoryless'' or ''memory'')', pa_mode);
end

%% 3. AWGN noise model
SNR_dB    = 30;
sig_pow   = mean(abs(x).^2);
noise_std = sqrt(sig_pow / 10^(SNR_dB/10) / 2);
add_noise = @(y) y + noise_std * (randn(size(y)) + 1j*randn(size(y)));
fprintf('Noise: SNR=%d dB  noise_std=%.4e\n', SNR_dB, noise_std);

%% 4. Psat probe
%
% Drive the PA at r_cap and observe max output.
% For memory PA, the filter slightly changes effective drive,
% but main tap = 1.0 so the difference is small.

G_target = 1.0;
r_cap    = 1.8;
r_fine   = linspace(0, r_cap, 100000)';
x_probe  = (r_cap ./ max(r, 1e-12)) .* x;

Psat       = max(abs(pa_model(x_probe))) * 0.99;
[~, idx]   = min(abs(amam(r_fine) - Psat));
r_drive_max = r_fine(idx);
fprintf('Psat=%.4f  r_drive_max=%.4f\n', Psat, r_drive_max);

%% 5. Iterative Direct DPD

mu    = 0.8;
N_avg = 8;

fprintf('\n=== Iterative Direct DPD ===\n');
fprintf('Config: %d iters, mu=%.2f, N_avg=%d, SNR=%d dB\n', n_iter, mu, N_avg, SNR_dB);

P = x;
mse_log = zeros(1, n_iter+1);

% Baseline: PA without DPD
Y_pa_only  = pa_model(x);
mse_no_dpd = mean(abs(Y_pa_only - G_target*x).^2);
fprintf('No DPD        MSE = %.3e\n', mse_no_dpd);

mse_log(1) = mse_no_dpd;
fprintf('Iter  0       MSE = %.3e\n', mse_log(1));

Psat_cur = Psat;
for k = 1:n_iter
    % (a) I/Q averaged capture
    M_avg = zeros(N, 1);
    for j = 1:N_avg
        M_avg = M_avg + add_noise(pa_model(P));
    end
    M_avg = M_avg / N_avg;

    % Update Psat from observation
    Psat_cur = max(Psat_cur, max(abs(M_avg)) * 0.99);

    % (b) Target: linear up to Psat, hard-clip beyond
    r_tgt = min(G_target * r, Psat_cur);
    T     = r_tgt ./ max(r, 1e-12) .* x;

    % (c) Complex ratio correction: corr = T / M_avg
    corr = T ./ max(abs(M_avg), 1e-12) .* exp(-1j*angle(M_avg));
    P    = P .* ((1 - mu) + mu * corr);

    % (d) Clamp to r_drive_max
    P = P .* min(1, r_drive_max ./ max(abs(P), 1e-12));

    % (e) Evaluate MSE (noise-free)
    Y_check      = pa_model(P);
    mse_log(k+1) = mean(abs(Y_check - G_target*x).^2);
    fprintf('Iter %2d       MSE = %.3e\n', k, mse_log(k+1));
end

P_final = P;
Y_dpd   = pa_model(P_final);

%% 6. Summary
fprintf('\n=== Summary ===\n');
fprintf('PA mode:     %s\n', pa_mode);
fprintf('No-DPD MSE = %.3e\n', mse_no_dpd);
fprintf('DPD    MSE = %.3e\n', mse_log(end));
fprintf('Improvement: %.1f dB\n', 10*log10(mse_no_dpd / mse_log(end)));

%% 7. Visualization

r_ax = linspace(0, 1, 500)';

% Sort data by input amplitude for cleaner plots
[r_sorted, sort_idx] = sort(r);
Y_dpd_sorted     = Y_dpd(sort_idx);
Y_pa_only_sorted = Y_pa_only(sort_idx);
P_sorted         = P_final(sort_idx);
x_sorted         = x(sort_idx);

hardclip_ax = min(G_target * r_ax, Psat);

scr = get(0, 'ScreenSize');
fig = figure('Position', scr, 'Name', sprintf('ILC-DPD — %s PA', pa_mode));
set(fig, 'WindowState', 'maximized');

% ---- (1) PA AM/AM ----
subplot(2,3,1);
if strcmp(pa_mode, 'memoryless')
    plot(r_ax, amam(r_ax), 'b-', 'LineWidth', 2.0); hold on;
else
    scatter(r_sorted, abs(Y_pa_only_sorted), 2, [0.1 0.3 0.8], ...
            'filled', 'MarkerFaceAlpha', 0.35); hold on;
end
plot(r_ax, G_target*r_ax, 'k--', 'LineWidth', 1.5);
xlabel('|x|'); ylabel('|y|');
title('PA AM/AM'); grid on;
if strcmp(pa_mode, 'memory')
    legend('PA output (memory spread)','Ideal linear','Location','northwest');
else
    legend('PA output','Ideal linear','Location','northwest');
end

% ---- (2) PA AM/PM ----
subplot(2,3,2);
if strcmp(pa_mode, 'memoryless')
    plot(r_ax, ampm_rad(r_ax)*180/pi, 'b-', 'LineWidth', 2.0); hold on;
else
    pa_phase = angle(Y_pa_only_sorted .* conj(x_sorted)) * 180/pi;
    pa_phase(r_sorted < 0.02) = NaN;
    scatter(r_sorted, pa_phase, 2, [0.1 0.3 0.8], ...
            'filled', 'MarkerFaceAlpha', 0.35); hold on;
end
yline(0, 'k--', 'LineWidth', 1.0);
xlabel('|x|'); ylabel('Phase (deg)');
title('PA AM/PM'); grid on;

% ---- (3) DPD+PA AM/AM ----
subplot(2,3,3);
scatter(r_sorted, abs(Y_pa_only_sorted), 2, [0.1 0.3 0.8], ...
        'filled', 'MarkerFaceAlpha', 0.20); hold on;
scatter(r_sorted, abs(Y_dpd_sorted), 2, [0.85 0.1 0.1], ...
        'filled', 'MarkerFaceAlpha', 0.40);
plot(r_ax, G_target*r_ax,   'k--', 'LineWidth', 1.5);
plot(r_ax, hardclip_ax,     'm:',  'LineWidth', 2.0);
legend('PA only','DPD+PA','Ideal linear','Hard-clip','Location','northwest');
xlabel('|x|'); ylabel('|Y|');
title('DPD+PA  AM/AM'); grid on;

% ---- (4) DPD+PA AM/PM ----
phase_dpd     = angle(Y_dpd_sorted .* conj(x_sorted)) * 180/pi;
phase_pa_only = angle(Y_pa_only_sorted .* conj(x_sorted)) * 180/pi;
phase_dpd(r_sorted < 0.02)     = NaN;
phase_pa_only(r_sorted < 0.02) = NaN;

subplot(2,3,4);
scatter(r_sorted, phase_pa_only, 2, [0.1 0.3 0.8], ...
        'filled', 'MarkerFaceAlpha', 0.20); hold on;
scatter(r_sorted, phase_dpd, 2, [0.85 0.1 0.1], ...
        'filled', 'MarkerFaceAlpha', 0.40);
yline(0, 'k--', 'LineWidth', 1.5);
legend('PA only','DPD+PA','Ideal (0 deg)');
xlabel('|x|'); ylabel('Phase (deg)');
title('DPD+PA  AM/PM  (residual)'); grid on;

% ---- (5) Convergence ----
subplot(2,3,5);
semilogy(0:n_iter, mse_log, 'ko-', 'MarkerFaceColor','k', 'LineWidth', 1.5); hold on;
yline(mse_no_dpd, 'b--', 'LineWidth', 1.0);
xlabel('Iteration'); ylabel('MSE');
title('Convergence'); grid on;
legend('DPD MSE','No-DPD MSE','Location','northeast');

% ---- (6) DPD gain & phase ----
dpd_gain  = abs(P_sorted) ./ max(r_sorted, 1e-12);
dpd_phase = angle(P_sorted .* conj(x_sorted)) * 180/pi;
dpd_gain(r_sorted < 0.02)  = NaN;
dpd_phase(r_sorted < 0.02) = NaN;

subplot(2,3,6);
yyaxis left;
scatter(r_sorted, dpd_gain, 2, [0.0 0.5 0.0], 'filled', 'MarkerFaceAlpha', 0.35);
ylabel('|G_{DPD}|  (gain)');
ylim([0.8 max(dpd_gain(~isnan(dpd_gain)))*1.1]);
yyaxis right;
scatter(r_sorted, dpd_phase, 2, [0.7 0.3 0.0], 'filled', 'MarkerFaceAlpha', 0.35);
ylabel('Phase (deg)');
xlabel('|x|');
title('DPD: gain & phase vs amplitude'); grid on;

sgtitle(sprintf(['ILC-DPD  |  PA: %s  |  ' ...
    '%d iters, \\mu=%.2f, N_{avg}=%d, SNR=%d dB  |  ' ...
    'MSE: %.2e \\rightarrow %.2e  (%.1f dB)'], ...
    pa_mode, n_iter, mu, N_avg, SNR_dB, mse_no_dpd, mse_log(end), ...
    10*log10(mse_no_dpd/mse_log(end))));

saveas(gcf, sprintf('IterativeDPD_%s.png', pa_mode));
fprintf('Saved -> IterativeDPD_%s.png\n', pa_mode);
