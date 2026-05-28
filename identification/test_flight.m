%% test_flight.m
%  Junta múltiplos logs de voo num único time-base contínuo e plota:
%  RCOU (motores), pqr (gyro), Acc (IMU), ATT (atitude), RCIN (rádio).
%
%  Cada log tem tempo absoluto independente (segundos desde boot da Pixhawk).
%  Aqui, cada log é shiftado pra começar imediatamente após o anterior + GAP.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '1_data'));

% ╔══════════════════════════════════════════════════════════════════╗
% ║  CONFIGURAÇÃO                                                    ║
% ╚══════════════════════════════════════════════════════════════════╝
LOG_DIR = fullfile(fileparts(mfilename('fullpath')), '1_data');

LOGS = { ...
    '4 25-05-2026 09-31-48.log-132954.mat', ...
    '6 31-12-1979 21-00-00.bin-116760.mat', ...
    '7 31-12-1979 21-00-00.bin-113223.mat'};

GAP = 5;   % segundos de "lacuna" entre logs (separação visual)

%% ====================== Junção dos logs ============================
J = struct();   % junto
J.time_IMU  = []; J.gyrX = []; J.gyrY = []; J.gyrZ = [];
J.accX = []; J.accY = []; J.accZ = [];
J.time_ATT  = []; J.roll = []; J.pitch = []; J.yaw = [];
J.time_RCOU = []; J.pwm1 = []; J.pwm2 = []; J.pwm3 = []; J.pwm4 = [];
J.time_RCIN = []; J.rc_roll = []; J.rc_pitch = []; J.rc_thr = []; J.rc_yaw = [];

boundaries = zeros(numel(LOGS)-1, 1);   % tempo (no eixo unificado) onde cada log termina
log_starts = zeros(numel(LOGS), 1);     % onde cada log começa no eixo unificado

t_cursor = 0;
for i = 1:numel(LOGS)
    L = load_log_data(fullfile(LOG_DIR, LOGS{i}));

    % Normaliza: subtrai t_start do próprio log e adiciona o cursor
    t0_log = min([L.time_IMU(1), L.time_ATT(1), L.time_RCOU(1)]);
    t_shift = t_cursor - t0_log;

    log_starts(i) = t_cursor;

    J.time_IMU  = [J.time_IMU;  L.time_IMU  + t_shift];
    J.gyrX = [J.gyrX; L.gyrX_raw]; J.gyrY = [J.gyrY; L.gyrY_raw]; J.gyrZ = [J.gyrZ; L.gyrZ_raw];
    J.accX = [J.accX; L.accX_raw]; J.accY = [J.accY; L.accY_raw]; J.accZ = [J.accZ; L.accZ_raw];

    J.time_ATT  = [J.time_ATT;  L.time_ATT  + t_shift];
    J.roll  = [J.roll;  L.roll_deg];
    J.pitch = [J.pitch; L.pitch_deg];
    J.yaw   = [J.yaw;   L.yaw_deg];

    J.time_RCOU = [J.time_RCOU; L.time_RCOU + t_shift];
    J.pwm1 = [J.pwm1; L.pwm1_raw]; J.pwm2 = [J.pwm2; L.pwm2_raw];
    J.pwm3 = [J.pwm3; L.pwm3_raw]; J.pwm4 = [J.pwm4; L.pwm4_raw];

    if ~isempty(L.time_RCIN)
        J.time_RCIN = [J.time_RCIN; L.time_RCIN + t_shift];
        J.rc_roll  = [J.rc_roll;  L.rcin_roll];
        J.rc_pitch = [J.rc_pitch; L.rcin_pitch];
        J.rc_thr   = [J.rc_thr;   L.rcin_throttle];
        J.rc_yaw   = [J.rc_yaw;   L.rcin_yaw];
    end

    % Avança cursor pro fim deste log + gap
    t_end_log = max([L.time_IMU(end), L.time_ATT(end), L.time_RCOU(end)]);
    t_cursor = t_end_log + t_shift + GAP;

    if i < numel(LOGS)
        boundaries(i) = t_end_log + t_shift;
    end

    fprintf('Log %d → janela unificada [%.1f, %.1f] s  (dur %.1f s)\n', ...
        i, log_starts(i), t_end_log + t_shift, t_end_log + t_shift - log_starts(i));
end

t_total = t_cursor - GAP;
fprintf('\nTotal: %.1f s de dados concatenados de %d logs.\n', t_total, numel(LOGS));

%% ====================== Salva como .mat reutilizável =================
% Formato compatível com a saída de load_log_data() — pra que qualquer
% script (identify_plant, validate_params, etc.) possa fazer:
%   tmp = load(LOG_FILE);  L = tmp.L;
% no lugar de:
%   L = load_log_data(LOG_FILE);
L = struct();
L.format = 'concat';
L.source = strjoin(LOGS, ' + ');
L.time_IMU = J.time_IMU;
L.gyrX_raw = J.gyrX; L.gyrY_raw = J.gyrY; L.gyrZ_raw = J.gyrZ;
L.accX_raw = J.accX; L.accY_raw = J.accY; L.accZ_raw = J.accZ;
L.time_ATT = J.time_ATT;
L.roll_deg = J.roll; L.pitch_deg = J.pitch; L.yaw_deg = J.yaw;
L.time_RCOU = J.time_RCOU;
L.pwm1_raw = J.pwm1; L.pwm2_raw = J.pwm2;
L.pwm3_raw = J.pwm3; L.pwm4_raw = J.pwm4;
L.time_RCIN = J.time_RCIN;
L.rcin_roll = J.rc_roll; L.rcin_pitch = J.rc_pitch;
L.rcin_throttle = J.rc_thr; L.rcin_yaw = J.rc_yaw;
L.time_GPS = [];           % concat não preserva GPS (cada log tem origem distinta)
L.boundaries = boundaries; % tempos onde cada log termina (no time-base unificado)
L.log_starts = log_starts; % tempos onde cada log começa
L.log_names  = LOGS;
L.gap        = GAP;

out_file = fullfile(LOG_DIR, 'logs_concat.mat');
save(out_file, 'L', '-v7.3');
fprintf('Salvo: %s\n', out_file);

%% ====================== Helper p/ marcar boundaries =================
mark_logs = @() arrayfun(@(t) xline(t, 'w--', 'Alpha', 0.5, ...
    'HandleVisibility','off'), boundaries);
mark_labels = @() arrayfun(@(i) text(log_starts(i)+5, 0, ...
    sprintf('LOG %d', i+3), 'Color', 'y', 'FontWeight', 'bold', ...
    'VerticalAlignment','top','HandleVisibility','off'), 1:numel(LOGS));

%% =========================================================================
%  Fig 1 — RCOU (PWM dos motores)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 400]);
plot(J.time_RCOU, J.pwm1, 'LineWidth', 0.8); hold on
plot(J.time_RCOU, J.pwm2, 'LineWidth', 0.8);
plot(J.time_RCOU, J.pwm3, 'LineWidth', 0.8);
plot(J.time_RCOU, J.pwm4, 'LineWidth', 0.8);
mark_logs();
xlabel('Tempo unificado [s]'); ylabel('PWM [\mus]');
title('RCOU — PWM enviado aos motores (todos os logs concatenados)');
legend('M1 (FR)','M2 (RL)','M3 (FL)','M4 (RR)','Location','best');
grid on;

%% =========================================================================
%  Fig 2 — RCIN (comando do piloto via rádio)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 600]);
subplot(4,1,1);
plot(J.time_RCIN, J.rc_roll, 'b'); grid on
mark_logs(); ylabel('Roll [\mus]'); ylim([1000 2000]); yline(1500,'k--');
title('RCIN — comando do piloto');
subplot(4,1,2);
plot(J.time_RCIN, J.rc_pitch, 'r'); grid on
mark_logs(); ylabel('Pitch [\mus]'); ylim([1000 2000]); yline(1500,'k--');
subplot(4,1,3);
plot(J.time_RCIN, J.rc_thr, 'Color',[0.1 0.6 0.1]); grid on
mark_logs(); ylabel('Throttle [\mus]'); ylim([1000 2000]); yline(1500,'k--');
subplot(4,1,4);
plot(J.time_RCIN, J.rc_yaw, 'm'); grid on
mark_logs(); ylabel('Yaw [\mus]'); ylim([1000 2000]); yline(1500,'k--');
xlabel('Tempo unificado [s]');

%% =========================================================================
%  Fig 3 — pqr (gyro)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 600]);
subplot(3,1,1)
plot(J.time_IMU, J.gyrX, 'b'); grid on; mark_logs();
title('pqr — velocidades angulares [rad/s]'); legend('p');
subplot(3,1,2)
plot(J.time_IMU, J.gyrY, 'r'); grid on; mark_logs(); legend('q');
subplot(3,1,3)
plot(J.time_IMU, J.gyrZ, 'Color',[0.1 0.6 0.1]); grid on; mark_logs(); legend('r');
xlabel('Tempo unificado [s]');

%% =========================================================================
%  Fig 4 — Acc (acelerações)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 600]);
subplot(3,1,1)
plot(J.time_IMU, J.accX, 'b'); grid on; mark_logs();
title('Aceleração [m/s²]'); legend('AccX');
subplot(3,1,2)
plot(J.time_IMU, J.accY, 'r'); grid on; mark_logs(); legend('AccY');
subplot(3,1,3)
plot(J.time_IMU, J.accZ, 'Color',[0.1 0.6 0.1]); grid on; mark_logs(); legend('AccZ');
xlabel('Tempo unificado [s]');

%% =========================================================================
%  Fig 5 — ATT (atitude em graus)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 600]);
subplot(3,1,1)
plot(J.time_ATT, J.roll, 'b'); grid on; mark_logs();
title('Atitude [deg]'); legend('Roll \phi');
subplot(3,1,2)
plot(J.time_ATT, J.pitch, 'r'); grid on; mark_logs(); legend('Pitch \theta');
subplot(3,1,3)
plot(J.time_ATT, J.yaw, 'Color',[0.1 0.6 0.1]); grid on; mark_logs(); legend('Yaw \psi');
xlabel('Tempo unificado [s]');

%% =========================================================================
%  Resumo
%  =========================================================================
fprintf('\n================== RESUMO LOGS CONCATENADOS ====================\n');
fprintf('Total: %.1f s (%.1f min) | %d logs | gap entre logs = %d s\n', ...
    t_total, t_total/60, numel(LOGS), GAP);
fprintf('IMU:  %d amostras  | dt_med=%.4f s\n', numel(J.time_IMU), median(diff(J.time_IMU(J.time_IMU<boundaries(1)))));
fprintf('ATT:  %d amostras\n', numel(J.time_ATT));
fprintf('RCOU: %d amostras\n', numel(J.time_RCOU));
fprintf('RCIN: %d amostras\n', numel(J.time_RCIN));
fprintf('===============================================================\n');
