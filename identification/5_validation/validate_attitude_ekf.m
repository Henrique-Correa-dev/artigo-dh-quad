% validate_attitude_ekf.m — Valida o EKF de atitude (roll/pitch) vs ATT do log
%
% Compara 3 curvas de phi/theta:
%   1. ATT do log         (EKF do ArduPilot — referência "verdade")
%   2. EKF próprio        (attitude_ekf.m — gyro+accel fundidos)
%   3. Integração open-loop (cinemática de Euler só com gyro — o baseline ruim)
%
% Mostra que o EKF próprio cura a deriva da integração pura.
% Roda nos sensores MEDIDOS (gyro+accel do IMU), não no modelo.

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ===================== CONFIG =====================
LOG_FILE = 'logs_concat.mat';     % ou um log cru individual
t_window = [75, 110];             % [s] janela de validação
CALM_WIN = [131, 150];            % [s] janela calma p/ estimar bias de gyro
                                  %     (use [] pra não estimar — bias=0)

%% ===================== CARREGAR =====================
L = load_log_data(fullfile(paths.data, LOG_FILE));
g = 9.81;

% Bias de gyro (de janela calma) — melhora a predição
bias_gyro = [0;0;0];
if ~isempty(CALM_WIN)
    Bc = estimate_bias(L, struct('windows_manual', CALM_WIN, 'verbose', false));
    if all(isfinite(Bc.gyro)), bias_gyro = Bc.gyro; end
    fprintf('Bias gyro (janela calma [%g,%g]): [%+.4f %+.4f %+.4f] rad/s\n', ...
        CALM_WIN(1), CALM_WIN(2), bias_gyro);
end

% Recortar janela na grade do IMU (gyro+accel mesma base de tempo)
idx = (L.time_IMU >= t_window(1)) & (L.time_IMU <= t_window(2));
t     = L.time_IMU(idx);
gyro  = [L.gyrX_raw(idx), L.gyrY_raw(idx), L.gyrZ_raw(idx)];
accel = [L.accX_raw(idx), L.accY_raw(idx), L.accZ_raw(idx)];
N = numel(t);

% ATT do log (referência) interpolado na grade do IMU
roll_ref  = interp1(L.time_ATT, L.roll_deg,  t, 'linear', 'extrap');
pitch_ref = interp1(L.time_ATT, L.pitch_deg, t, 'linear', 'extrap');

fprintf('Janela [%g,%g]s: %d amostras IMU (rate≈%.0f Hz)\n', ...
    t_window(1), t_window(2), N, N/(t(end)-t(1)));

%% ===================== EKF PRÓPRIO =====================
est = attitude_ekf(t, gyro, accel, struct('g', g, 'bias_gyro', bias_gyro));

%% ===================== BASELINE: integração open-loop =====================
% Cinemática de Euler só com gyro (sem correção do accel) — o que diverge.
gyro_c = gyro - bias_gyro(:)';
phi_ol = zeros(N,1); theta_ol = zeros(N,1);
phi_ol(1)   = est.phi(1);     % mesma IC do EKF (justo)
theta_ol(1) = est.theta(1);
for k = 1:N-1
    dt = t(k+1) - t(k);  if dt<=0||~isfinite(dt), dt=1e-3; end
    p = gyro_c(k,1); q = gyro_c(k,2); r = gyro_c(k,3);
    ph = phi_ol(k); th = theta_ol(k);
    phi_ol(k+1)   = ph + dt*(p + (q*sin(ph)+r*cos(ph))*tan(th));
    theta_ol(k+1) = th + dt*(q*cos(ph) - r*sin(ph));
end
phi_ol = rad2deg(phi_ol);  theta_ol = rad2deg(theta_ol);

%% ===================== MÉTRICAS =====================
rmse = @(a,b) sqrt(mean((a-b).^2));
fprintf('\n=== RMSE vs ATT do log (graus) ===\n');
fprintf('  %-12s  %8s  %8s\n', '', 'phi', 'theta');
fprintf('  %-12s  %8.3f  %8.3f\n', 'EKF proprio', ...
    rmse(est.phi_deg, roll_ref), rmse(est.theta_deg, pitch_ref));
fprintf('  %-12s  %8.3f  %8.3f\n', 'Open-loop', ...
    rmse(phi_ol, roll_ref), rmse(theta_ol, pitch_ref));

%% ===================== PLOT =====================
fig = figure('Name','validate_attitude_ekf','Position',[80 80 1300 700]);

subplot(2,1,1); hold on; grid on;
plot(t, roll_ref,    'b-',  'LineWidth', 1.4, 'DisplayName', 'ATT log (ref)');
plot(t, est.phi_deg, 'r--', 'LineWidth', 1.4, 'DisplayName', 'EKF proprio');
plot(t, phi_ol,      'm:',  'LineWidth', 1.2, 'DisplayName', 'open-loop (gyro)');
ylabel('\phi (°)'); legend('Location','best');
title(sprintf('Roll — RMSE EKF=%.2f° vs open-loop=%.2f°', ...
    rmse(est.phi_deg, roll_ref), rmse(phi_ol, roll_ref)));
ylim([-60 60]);

subplot(2,1,2); hold on; grid on;
plot(t, pitch_ref,     'b-',  'LineWidth', 1.4, 'DisplayName', 'ATT log (ref)');
plot(t, est.theta_deg, 'r--', 'LineWidth', 1.4, 'DisplayName', 'EKF proprio');
plot(t, theta_ol,      'm:',  'LineWidth', 1.2, 'DisplayName', 'open-loop (gyro)');
ylabel('\theta (°)'); xlabel('t [s]'); legend('Location','best');
title(sprintf('Pitch — RMSE EKF=%.2f° vs open-loop=%.2f°', ...
    rmse(est.theta_deg, pitch_ref), rmse(theta_ol, pitch_ref)));
ylim([-60 60]);

sgtitle(sprintf('attitude\\_ekf vs ATT do log — janela [%g, %g]s', t_window(1), t_window(2)));
saveas(fig, fullfile(paths.images, 'validate_attitude_ekf.png'));
fprintf('\nFigura: %s\n', fullfile(paths.images, 'validate_attitude_ekf.png'));
