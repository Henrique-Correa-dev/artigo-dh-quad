%ANALYZE_MZ  Análise passo-a-passo do Mz (dinâmica de yaw).
%
% Hipótese investigada: o modelo está convergindo k_Q1 e k_Q3 pra valores
% normais, mas zerando k_Q2 e k_Q4 (LB). Isso pode indicar:
%   - Motores físicos diferentes (M1+M3 vs M2+M4)
%   - Fórmula CW/CCW errada
%   - Algum acoplamento não modelado
%
% Esta análise usa a janela 0-45s do logs_concat (dominada por manobras
% de yaw — fica mais fácil ver a relação entre PWM e r).
%
% PASSO 1 (atual): apenas visualizar
%   - r medido (gyro Z) ao longo do tempo
%   - PWM dos 4 motores ao longo do tempo
%
% Próximos passos virão depois.

clear; clc; close all;

addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();

% ╔══════════════════════════════════════════════════════════════════╗
% ║  CONFIGURAÇÃO                                                    ║
% ╚══════════════════════════════════════════════════════════════════╝
LOG_FILE = 'logs_concat.mat';
T_WINDOW = [0, 45];

%% ====== Carrega dados ======
L = load_log_data(fullfile(setup_paths().data, LOG_FILE));

idx_imu  = (L.time_IMU  >= T_WINDOW(1)) & (L.time_IMU  <= T_WINDOW(2));
idx_rcou = (L.time_RCOU >= T_WINDOW(1)) & (L.time_RCOU <= T_WINDOW(2));

t_imu  = L.time_IMU(idx_imu);
r_meas = L.gyrZ_raw(idx_imu);

t_rcou = L.time_RCOU(idx_rcou);
pwm1 = L.pwm1_raw(idx_rcou);
pwm2 = L.pwm2_raw(idx_rcou);
pwm3 = L.pwm3_raw(idx_rcou);
pwm4 = L.pwm4_raw(idx_rcou);

fprintf('Janela: [%.1f, %.1f] s (%.1f s)\n', T_WINDOW(1), T_WINDOW(2), diff(T_WINDOW));
fprintf('Amostras IMU: %d  |  RCOU: %d\n', numel(t_imu), numel(t_rcou));
fprintf('r range:    [%+.3f, %+.3f] rad/s\n', min(r_meas), max(r_meas));
fprintf('PWM ranges:\n');
fprintf('  M1: [%4d, %4d]\n', round(min(pwm1)), round(max(pwm1)));
fprintf('  M2: [%4d, %4d]\n', round(min(pwm2)), round(max(pwm2)));
fprintf('  M3: [%4d, %4d]\n', round(min(pwm3)), round(max(pwm3)));
fprintf('  M4: [%4d, %4d]\n', round(min(pwm4)), round(max(pwm4)));

%% ====== Plot ======
fig = figure('Color','w','Position',[80 80 1300 700]);

subplot(2,1,1); hold on; grid on;
plot(t_imu, r_meas, 'b-', 'LineWidth', 1.2);
yline(0, 'k--', 'HandleVisibility','off');
xlabel('t [s]'); ylabel('r (yaw rate) [rad/s]');
title(sprintf('r medido — janela [%g, %g]s', T_WINDOW(1), T_WINDOW(2)));

subplot(2,1,2); hold on; grid on;
plot(t_rcou, pwm1, 'LineWidth', 1.2, 'DisplayName', 'M1 (FR)');
plot(t_rcou, pwm2, 'LineWidth', 1.2, 'DisplayName', 'M2 (RL)');
plot(t_rcou, pwm3, 'LineWidth', 1.2, 'DisplayName', 'M3 (FL)');
plot(t_rcou, pwm4, 'LineWidth', 1.2, 'DisplayName', 'M4 (RR)');
xlabel('t [s]'); ylabel('PWM [\mus]');
title('PWM dos motores (RCOU)');
legend('Location','best');

sgtitle('analyze\_mz — Passo 1: r e PWM');

img_path = fullfile(setup_paths().images, 'analyze_mz_step1.png');
saveas(fig, img_path);
fprintf('\nFigura salva: %s\n', img_path);
