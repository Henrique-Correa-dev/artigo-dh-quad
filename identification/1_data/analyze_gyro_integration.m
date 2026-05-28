%ANALYZE_GYRO_INTEGRATION  Integra p, q, r medidos via cinemática de Euler e
%                          compara com a atitude do EKF (Pixhawk).
%
% Pergunta: se eu integrar o gyro puro, bate com o que o EKF do Pixhawk
%           reporta como atitude?
%
% Resposta esperada: NÃO bate perfeitamente. Diferença = (drift do gyro) +
%                    (correção da fusão sensorial accel+mag+GPS do EKF).
%
% Cinemática usada (mesma do vtol_dynamics.m):
%   phi_dot   = p + (q·sin(phi) + r·cos(phi)) · tan(theta)
%   theta_dot = q·cos(phi) - r·sin(phi)
%   psi_dot   = (q·sin(phi) + r·cos(phi)) / cos(theta)
%
% Integração: RK4 com sub-step (igual identify_plant.m).
%
% USO:
%   >> analyze_gyro_integration

clear; clc; close all;

addpath(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), '1_data'));

% ╔══════════════════════════════════════════════════════════════════╗
% ║  CONFIGURAÇÃO                                                    ║
% ╚══════════════════════════════════════════════════════════════════╝
LOG_FILE  = 'logs_concat.mat';
T_WINDOW  = [34, 124];        % segundos (time-base do concat)
DT_INTEG  = 0.005;            % passo de integração (s) — fino p/ RK4

% Bias do gyro (rad/s) — se você já rodou estimate_bias, cole aqui:
GYRO_BIAS = [0; 0; 0];        % [bias_x; bias_y; bias_z]
% Exemplo: GYRO_BIAS = [-0.0038; +0.0035; -0.0015];

%% ====== Carrega log ======
L = load_log_data(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '1_data', LOG_FILE));

% Recorta janela
idx_imu = (L.time_IMU >= T_WINDOW(1)) & (L.time_IMU <= T_WINDOW(2));
idx_att = (L.time_ATT >= T_WINDOW(1)) & (L.time_ATT <= T_WINDOW(2));

t_imu = L.time_IMU(idx_imu);
p_meas = L.gyrX_raw(idx_imu) - GYRO_BIAS(1);
q_meas = L.gyrY_raw(idx_imu) - GYRO_BIAS(2);
r_meas = L.gyrZ_raw(idx_imu) - GYRO_BIAS(3);

t_att = L.time_ATT(idx_att);
phi_ekf   = deg2rad(L.roll_deg(idx_att));
theta_ekf = deg2rad(L.pitch_deg(idx_att));
psi_ekf   = deg2rad(L.yaw_deg(idx_att));

fprintf('Janela: [%.1f, %.1f] s (%.1f s)\n', T_WINDOW(1), T_WINDOW(2), diff(T_WINDOW));
fprintf('Amostras IMU: %d  |  ATT: %d\n', numel(t_imu), numel(t_att));
fprintf('Bias do gyro aplicado: [%+.5f, %+.5f, %+.5f] rad/s\n', GYRO_BIAS);

%% ====== Grade de integração e interpolação dos sinais ======
t_int = (T_WINDOW(1):DT_INTEG:T_WINDOW(2))';
N = numel(t_int);

% Interpola gyro pra grade de integração
p_int = interp1(t_imu, p_meas, t_int, 'linear', 'extrap');
q_int = interp1(t_imu, q_meas, t_int, 'linear', 'extrap');
r_int = interp1(t_imu, r_meas, t_int, 'linear', 'extrap');

% Condição inicial = atitude do EKF no instante T_WINDOW(1)
phi_0   = interp1(t_att, phi_ekf,   t_int(1), 'linear', 'extrap');
theta_0 = interp1(t_att, theta_ekf, t_int(1), 'linear', 'extrap');
psi_0   = interp1(t_att, psi_ekf,   t_int(1), 'linear', 'extrap');

fprintf('\nCondição inicial (do EKF em t=%.2fs):\n', t_int(1));
fprintf('  phi_0   = %+.3f deg\n', rad2deg(phi_0));
fprintf('  theta_0 = %+.3f deg\n', rad2deg(theta_0));
fprintf('  psi_0   = %+.3f deg\n', rad2deg(psi_0));

%% ====== Integração RK4 da cinemática de Euler ======
phi_sim   = zeros(N, 1);  phi_sim(1)   = phi_0;
theta_sim = zeros(N, 1);  theta_sim(1) = theta_0;
psi_sim   = zeros(N, 1);  psi_sim(1)   = psi_0;

% Função p/ derivada (mesma fórmula do vtol_dynamics.m)
%   pqr_at(t) interpolado linearmente entre amostras
pqr_func = @(t) [interp1(t_int, p_int, t, 'linear', 'extrap'); ...
                 interp1(t_int, q_int, t, 'linear', 'extrap'); ...
                 interp1(t_int, r_int, t, 'linear', 'extrap')];

f_euler = @(t, y) euler_kinematics(t, y, pqr_func);

for k = 1:N-1
    y = [phi_sim(k); theta_sim(k); psi_sim(k)];
    t = t_int(k);
    h = DT_INTEG;
    k1 = f_euler(t,         y);
    k2 = f_euler(t + h/2,   y + h/2*k1);
    k3 = f_euler(t + h/2,   y + h/2*k2);
    k4 = f_euler(t + h,     y + h*k3);
    y_next = y + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    phi_sim(k+1)   = y_next(1);
    theta_sim(k+1) = y_next(2);
    psi_sim(k+1)   = y_next(3);
end

%% ====== Métricas de drift ======
phi_ekf_int   = interp1(t_att, phi_ekf,   t_int, 'linear', 'extrap');
theta_ekf_int = interp1(t_att, theta_ekf, t_int, 'linear', 'extrap');
psi_ekf_int   = interp1(t_att, psi_ekf,   t_int, 'linear', 'extrap');

err_phi   = rad2deg(phi_sim   - phi_ekf_int);
err_theta = rad2deg(theta_sim - theta_ekf_int);
err_psi   = rad2deg(psi_sim   - psi_ekf_int);

fprintf('\n=== Drift acumulado em %.1f s ===\n', diff(T_WINDOW));
fprintf('  phi    (roll):  inicial=%+.2f°  final=%+.2f°  delta=%+.2f°  RMS=%.2f°\n', ...
    err_phi(1), err_phi(end), err_phi(end)-err_phi(1), rms(err_phi));
fprintf('  theta (pitch):  inicial=%+.2f°  final=%+.2f°  delta=%+.2f°  RMS=%.2f°\n', ...
    err_theta(1), err_theta(end), err_theta(end)-err_theta(1), rms(err_theta));
fprintf('  psi    (yaw):   inicial=%+.2f°  final=%+.2f°  delta=%+.2f°  RMS=%.2f°\n', ...
    err_psi(1), err_psi(end), err_psi(end)-err_psi(1), rms(err_psi));

%% ====== Plots ======
fig = figure('Color','w','Position',[80 60 1200 800]);

subplot(3,2,1); hold on; grid on;
plot(t_int, rad2deg(phi_ekf_int), 'b-', 'LineWidth', 1.2, 'DisplayName','EKF (Pixhawk)');
plot(t_int, rad2deg(phi_sim),     'r--','LineWidth', 1.2, 'DisplayName','Integrando gyro');
ylabel('\phi (roll) [°]'); legend('Location','best');
title('Roll: integração de p vs EKF');

subplot(3,2,3); hold on; grid on;
plot(t_int, rad2deg(theta_ekf_int), 'b-', 'LineWidth', 1.2, 'DisplayName','EKF');
plot(t_int, rad2deg(theta_sim),     'r--','LineWidth', 1.2, 'DisplayName','Integrando gyro');
ylabel('\theta (pitch) [°]');
title('Pitch: integração de q vs EKF');

subplot(3,2,5); hold on; grid on;
plot(t_int, rad2deg(psi_ekf_int), 'b-', 'LineWidth', 1.2, 'DisplayName','EKF');
plot(t_int, rad2deg(psi_sim),     'r--','LineWidth', 1.2, 'DisplayName','Integrando gyro');
ylabel('\psi (yaw) [°]'); xlabel('t [s]');
title('Yaw: integração de r vs EKF');

subplot(3,2,2); hold on; grid on;
plot(t_int, err_phi, 'm-', 'LineWidth', 1);
yline(0,'k--'); ylabel('\Delta\phi [°]');
title(sprintf('Drift de roll (RMS=%.2f°, final=%+.2f°)', rms(err_phi), err_phi(end)));

subplot(3,2,4); hold on; grid on;
plot(t_int, err_theta, 'm-', 'LineWidth', 1);
yline(0,'k--'); ylabel('\Delta\theta [°]');
title(sprintf('Drift de pitch (RMS=%.2f°, final=%+.2f°)', rms(err_theta), err_theta(end)));

subplot(3,2,6); hold on; grid on;
plot(t_int, err_psi, 'm-', 'LineWidth', 1);
yline(0,'k--'); ylabel('\Delta\psi [°]'); xlabel('t [s]');
title(sprintf('Drift de yaw (RMS=%.2f°, final=%+.2f°)', rms(err_psi), err_psi(end)));

sgtitle(sprintf('Integração do gyro vs EKF — janela [%g, %g]s (drift = bias × t + ruído)', ...
    T_WINDOW(1), T_WINDOW(2)));

img_path = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'outputs', 'images', 'gyro_integration_vs_ekf.png');
saveas(fig, img_path);
fprintf('\nFigura salva em: %s\n', img_path);

%% ====== Interpretação automática ======
fprintf('\n=== INTERPRETAÇÃO ===\n');
dur = diff(T_WINDOW);
implied_bias = [(rad2deg(phi_sim(end)) - rad2deg(phi_ekf_int(end))) - ...
                (rad2deg(phi_sim(1))   - rad2deg(phi_ekf_int(1))); ...
                (rad2deg(theta_sim(end))-rad2deg(theta_ekf_int(end))) - ...
                (rad2deg(theta_sim(1)) - rad2deg(theta_ekf_int(1))); ...
                (rad2deg(psi_sim(end)) - rad2deg(psi_ekf_int(end))) - ...
                (rad2deg(psi_sim(1))   - rad2deg(psi_ekf_int(1)))] / dur;
fprintf('Bias do gyro IMPLICITO (drift/duração):\n');
fprintf('  bias_x ≈ %+.5f rad/s  (%+.3f°/s)\n', deg2rad(implied_bias(1)), implied_bias(1));
fprintf('  bias_y ≈ %+.5f rad/s  (%+.3f°/s)\n', deg2rad(implied_bias(2)), implied_bias(2));
fprintf('  bias_z ≈ %+.5f rad/s  (%+.3f°/s)\n', deg2rad(implied_bias(3)), implied_bias(3));

fprintf('\nNOTA: o drift NÃO é só bias — também inclui:\n');
fprintf('  - Correção do EKF via accel (vetor gravidade) sobre roll/pitch\n');
fprintf('  - Correção do EKF via mag/GPS sobre yaw\n');
fprintf('  - Erros de Euler (gimbal lock em theta=±90°) se atitude for grande\n');
fprintf('\nPra ver bias PURO, use estimate_bias.m em janela com drone PARADO.\n');


%% ====== Subfunção ======
function dy = euler_kinematics(t, y, pqr_func)
    phi   = y(1);
    theta = y(2);
    pqr = pqr_func(t);
    p = pqr(1); q = pqr(2); r = pqr(3);

    sin_phi = sin(phi); cos_phi = cos(phi);
    cos_theta = cos(theta);
    if abs(cos_theta) < 1e-7
        cos_theta = 1e-7 * sign(cos_theta + 1e-12);   % evita div-0
    end
    tan_theta = sin(theta) / cos_theta;

    phi_dot   = p + (q*sin_phi + r*cos_phi) * tan_theta;
    theta_dot = q*cos_phi - r*sin_phi;
    psi_dot   = (q*sin_phi + r*cos_phi) / cos_theta;

    dy = [phi_dot; theta_dot; psi_dot];
end
