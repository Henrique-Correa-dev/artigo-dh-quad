%% motor_mapping_analysis.m
% =========================================================================
% Analise de Mapeamento de Motores - DH Hybrid-Drone VTOL
% =========================================================================
% Objetivo: Identificar qual canal PWM do log ArduPilot (C1-C4) corresponde
%           a qual posicao fisica do motor no frame H.
%
% Geometria do frame (vista de cima, nariz para cima):
%
%             NARIZ (frente)
%                |
%       M? o----+----o M?        <- frontais (311.18 mm do CG)
%          |    |    |
%     232mm|    CG   |232mm
%          |         |
%       M? o---------o M?        <- traseiros (342.87 mm do CG)
%
% Metodo:
%   1. Carregar dados brutos do log ArduPilot
%   2. Interpolar todos os sinais para grade temporal comum (10 Hz)
%   3. Calcular derivadas angulares (p_dot, q_dot, r_dot)
%   4. Correlacao cruzada entre variacao de PWM e derivadas angulares
%   5. Interpretar: qual motor esta em qual posicao
%
% Convencoes (NED, corpo):
%   - p_dot positivo = aceleracao de roll para direita (asa direita desce)
%     -> Causado por aumento de empuxo no lado ESQUERDO
%   - q_dot positivo = aceleracao de pitch para cima (nariz sobe)
%     -> Causado por aumento de empuxo TRASEIRO
%   - r_dot positivo = aceleracao de yaw para direita (nariz vira para direita)
%     -> Causado por desbalanco de torque reativo
%
% Autor: Henrique / Claude
% Data: 2026-03-15
% =========================================================================

clear; close all; clc;

%% ========================================================================
%  1. CARREGAR DADOS BRUTOS
% =========================================================================
fprintf('=== Carregando dados do log ArduPilot ===\n');
load(fullfile('..', 'identification', 'log_data.mat'));

% Converter timestamps de microsegundos para segundos
ATT.TimeS  = double(ATT.TimeUS) / 1e6;
IMU.TimeS  = double(IMU.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;
GPS.TimeS  = double(GPS.TimeUS) / 1e6;

fprintf('  IMU:  %d amostras, %.1f - %.1f s\n', height(IMU), min(IMU.TimeS), max(IMU.TimeS));
fprintf('  ATT:  %d amostras, %.1f - %.1f s\n', height(ATT), min(ATT.TimeS), max(ATT.TimeS));
fprintf('  RCOU: %d amostras, %.1f - %.1f s\n', height(RCOU), min(RCOU.TimeS), max(RCOU.TimeS));
fprintf('  GPS:  %d amostras, %.1f - %.1f s\n', height(GPS), min(GPS.TimeS), max(GPS.TimeS));

%% ========================================================================
%  2. EXTRAIR SINAIS BRUTOS
% =========================================================================

% IMU - usar instancia 0
idx_imu = IMU.I == 0;
time_IMU  = IMU.TimeS(idx_imu);
gyrX_raw  = IMU.GyrX(idx_imu);   % p [rad/s]
gyrY_raw  = IMU.GyrY(idx_imu);   % q [rad/s]
gyrZ_raw  = IMU.GyrZ(idx_imu);   % r [rad/s]
accX_raw  = IMU.AccX(idx_imu);   % [m/s^2]
accY_raw  = IMU.AccY(idx_imu);
accZ_raw  = IMU.AccZ(idx_imu);

% ATT - atitude do EKF
time_ATT   = ATT.TimeS;
roll_raw   = ATT.Roll;    % [deg]
pitch_raw  = ATT.Pitch;   % [deg]
yaw_raw    = ATT.Yaw;     % [deg]

% RCOU - sinais PWM dos motores
time_RCOU  = RCOU.TimeS;
pwm1_raw   = double(RCOU.C1);  % Canal 1 [us]
pwm2_raw   = double(RCOU.C2);  % Canal 2 [us]
pwm3_raw   = double(RCOU.C3);  % Canal 3 [us]
pwm4_raw   = double(RCOU.C4);  % Canal 4 [us]

% GPS
time_GPS = GPS.TimeS;
lat_raw  = GPS.Lat;
lon_raw  = GPS.Lng;
alt_raw  = GPS.Alt;

%% ========================================================================
%  3. INTERPOLACAO PARA GRADE TEMPORAL COMUM (10 Hz)
% =========================================================================
fprintf('\n=== Interpolacao para grade comum (10 Hz) ===\n');

% Definir tempo comum: sobreposicao de todos os sensores
t_start = max([min(time_IMU), min(time_ATT), min(time_GPS), min(time_RCOU)]);
t_end   = min([max(time_IMU), max(time_ATT), max(time_GPS), max(time_RCOU)]);
dt = 0.1;  % 10 Hz
t = (t_start : dt : t_end)';

fprintf('  Intervalo: %.1f - %.1f s (%.0f amostras)\n', t_start, t_end, length(t));

% IMU
p = interp1(time_IMU, gyrX_raw, t, 'linear');  % roll rate [rad/s]
q = interp1(time_IMU, gyrY_raw, t, 'linear');  % pitch rate [rad/s]
r = interp1(time_IMU, gyrZ_raw, t, 'linear');  % yaw rate [rad/s]
accX = interp1(time_IMU, accX_raw, t, 'linear');
accY = interp1(time_IMU, accY_raw, t, 'linear');
accZ = interp1(time_IMU, accZ_raw, t, 'linear');

% ATT
roll  = interp1(time_ATT, roll_raw,  t, 'linear');  % [deg]
pitch = interp1(time_ATT, pitch_raw, t, 'linear');  % [deg]
yaw   = interp1(time_ATT, yaw_raw,   t, 'linear');  % [deg]

% RCOU
pwm1 = interp1(time_RCOU, pwm1_raw, t, 'linear');
pwm2 = interp1(time_RCOU, pwm2_raw, t, 'linear');
pwm3 = interp1(time_RCOU, pwm3_raw, t, 'linear');
pwm4 = interp1(time_RCOU, pwm4_raw, t, 'linear');

% GPS -> NED
lat = interp1(time_GPS, lat_raw, t, 'linear');
lon = interp1(time_GPS, lon_raw, t, 'linear');
alt = interp1(time_GPS, alt_raw, t, 'linear');

%% ========================================================================
%  4. DEFINIR TRECHO DE ANALISE
% =========================================================================
% Janela de analise manual — ajuste aqui para isolar trechos com mais
% excitacao ou melhor qualidade de dados.
t_analise_start = 150;   % [s] inicio da janela
t_analise_end   = 200;   % [s] fim da janela

fprintf('\n=== Trecho de analise: %.0f - %.0f s ===\n', t_analise_start, t_analise_end);

t_vtol_start = t_analise_start;
t_vtol_end   = t_analise_end;

%% ========================================================================
%  5. VISAO GERAL DOS DADOS (TODO O LOG)
% =========================================================================
fprintf('\n=== Gerando figuras de visao geral ===\n');

figure('Name', 'Visao Geral - PWM', 'Position', [50 50 1200 700]);

subplot(5,1,1);
plot(t, pwm1, 'LineWidth', 1); ylabel('C1 [us]'); grid on;
title('Sinais PWM do Log (canais C1-C4 do RCOU)');
xline(t_vtol_start, 'g--', 'VTOL inicio'); xline(t_vtol_end, 'r--', 'VTOL fim');

subplot(5,1,2);
plot(t, pwm2, 'LineWidth', 1); ylabel('C2 [us]'); grid on;
xline(t_vtol_start, 'g--'); xline(t_vtol_end, 'r--');

subplot(5,1,3);
plot(t, pwm3, 'LineWidth', 1); ylabel('C3 [us]'); grid on;
xline(t_vtol_start, 'g--'); xline(t_vtol_end, 'r--');

subplot(5,1,4);
plot(t, pwm4, 'LineWidth', 1); ylabel('C4 [us]'); grid on;
xline(t_vtol_start, 'g--'); xline(t_vtol_end, 'r--');

subplot(5,1,5);
plot(t, roll, 'b', t, pitch, 'r', t, yaw, 'Color', [0.4 0.7 0.2], 'LineWidth', 1);
ylabel('[deg]'); xlabel('Tempo [s]'); grid on;
legend('Roll', 'Pitch', 'Yaw', 'Location', 'best');
xline(t_vtol_start, 'g--'); xline(t_vtol_end, 'r--');

saveas(gcf, 'fig01_visao_geral_pwm.png');

%% ========================================================================
%  6. ANALISE DE CORRELACAO - MAPEAMENTO DE MOTORES
% =========================================================================
fprintf('\n=== Analise de Correlacao para Mapeamento de Motores ===\n');

% Usar trecho completo de voo VTOL
idx_vtol = (t >= t_vtol_start) & (t <= t_vtol_end);

t_v   = t(idx_vtol);
p_v   = p(idx_vtol);
q_v   = q(idx_vtol);
r_v   = r(idx_vtol);
pwm_v = [pwm1(idx_vtol), pwm2(idx_vtol), pwm3(idx_vtol), pwm4(idx_vtol)];

% Derivadas angulares suavizadas
smooth_win = 5;
p_dot = gradient(movmean(p_v, smooth_win), dt);
q_dot = gradient(movmean(q_v, smooth_win), dt);
r_dot = gradient(movmean(r_v, smooth_win), dt);

% Remover media dos PWM (foco na variacao)
pwm_det = pwm_v - mean(pwm_v);

% Matriz de correlacao: 4 canais x 3 eixos
corr_matrix = zeros(4, 3);
labels_ch  = {'C1', 'C2', 'C3', 'C4'};
labels_axis = {'p_dot (Roll)', 'q_dot (Pitch)', 'r_dot (Yaw)'};

for i = 1:4
    corr_matrix(i, 1) = corr(pwm_det(:,i), p_dot);
    corr_matrix(i, 2) = corr(pwm_det(:,i), q_dot);
    corr_matrix(i, 3) = corr(pwm_det(:,i), r_dot);
end

% Exibir resultados
fprintf('\n  Correlacao normalizada: PWM vs derivadas angulares\n');
fprintf('  Trecho VTOL: %.1f - %.1f s\n\n', t_vtol_start, t_vtol_end);
fprintf('  %6s  |  %14s  %14s  %14s\n', '', labels_axis{:});
fprintf('  %s\n', repmat('-', 1, 62));
for i = 1:4
    fprintf('  %6s  |  %+14.4f  %+14.4f  %+14.4f\n', ...
        labels_ch{i}, corr_matrix(i,1), corr_matrix(i,2), corr_matrix(i,3));
end

%% ========================================================================
%  7. INTERPRETACAO AUTOMATICA
% =========================================================================
fprintf('\n\n=== INTERPRETACAO ===\n');
fprintf('  Convencoes (NED, corpo rígido):\n');
fprintf('  - corr(PWM, p_dot) > 0 -> motor ESQUERDO (empuxo esq sobe -> roll dir +)\n');
fprintf('  - corr(PWM, p_dot) < 0 -> motor DIREITO\n');
fprintf('  - corr(PWM, q_dot) > 0 -> motor TRASEIRO (empuxo tras sobe -> pitch up +)\n');
fprintf('  - corr(PWM, q_dot) < 0 -> motor FRONTAL\n');
fprintf('  - corr(r_dot) -> indica sentido de rotacao CW/CCW\n\n');

pos_labels = cell(1,4);
for i = 1:4
    lat_str = 'ESQ';
    lon_str = 'FRONTAL';

    if corr_matrix(i,1) < 0
        lat_str = 'DIR';
    end
    if corr_matrix(i,2) < 0
        lon_str = 'TRASEIRO';
    end

    if corr_matrix(i,3) > 0
        yaw_str = 'CW';
    else
        yaw_str = 'CCW';
    end

    pos_labels{i} = sprintf('%s-%s', lon_str, lat_str);
    fprintf('  %s -> %-16s | Rotacao: %s (visto de cima)\n', ...
        labels_ch{i}, pos_labels{i}, yaw_str);
end

%% ========================================================================
%  8. FIGURAS DE DIAGNOSTICO
% =========================================================================

% --- Figura: Heatmap de correlacao ---
figure('Name', 'Correlacao Heatmap', 'Position', [200 200 550 400]);
imagesc(corr_matrix');
colorbar;
caxis([-1 1]);

% Colormap vermelho-branco-azul
n_cmap = 256;
r_c = [linspace(0,1,n_cmap/2), ones(1,n_cmap/2)];
g_c = [linspace(0,1,n_cmap/2), linspace(1,0,n_cmap/2)];
b_c = [ones(1,n_cmap/2), linspace(1,0,n_cmap/2)];
colormap([r_c', g_c', b_c']);

set(gca, 'XTick', 1:4, 'XTickLabel', labels_ch);
set(gca, 'YTick', 1:3, 'YTickLabel', {'Roll (p\_dot)', 'Pitch (q\_dot)', 'Yaw (r\_dot)'});
title('Correlacao: \DeltaPWM vs Derivadas Angulares');

for i = 1:4
    for j = 1:3
        text(i, j, sprintf('%.3f', corr_matrix(i,j)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
            'FontSize', 12, 'Color', 'k');
    end
end
saveas(gcf, 'fig02_correlacao_heatmap.png');

% --- Figura: Sinais PWM vs derivadas angulares ---
figure('Name', 'PWM vs Derivadas', 'Position', [50 50 1200 900]);
colors = lines(4);

subplot(4,1,1);
plot(t_v, pwm_v, 'LineWidth', 1.2);
legend('C1','C2','C3','C4', 'Location', 'best');
ylabel('PWM [us]'); grid on;
title('Sinais PWM durante voo VTOL');

subplot(4,1,2);
yyaxis left;
plot(t_v, p_dot, 'k-', 'LineWidth', 1.5); ylabel('p\_dot [rad/s^2]');
yyaxis right;
for i = 1:4
    plot(t_v, pwm_det(:,i), '--', 'Color', colors(i,:), 'LineWidth', 0.8); hold on;
end
ylabel('\DeltaPWM [us]'); grid on;
legend('p\_dot','C1','C2','C3','C4', 'Location', 'best');
title('Roll: p\_dot vs \DeltaPWM');

subplot(4,1,3);
yyaxis left;
plot(t_v, q_dot, 'k-', 'LineWidth', 1.5); ylabel('q\_dot [rad/s^2]');
yyaxis right;
for i = 1:4
    plot(t_v, pwm_det(:,i), '--', 'Color', colors(i,:), 'LineWidth', 0.8); hold on;
end
ylabel('\DeltaPWM [us]'); grid on;
legend('q\_dot','C1','C2','C3','C4', 'Location', 'best');
title('Pitch: q\_dot vs \DeltaPWM');

subplot(4,1,4);
yyaxis left;
plot(t_v, r_dot, 'k-', 'LineWidth', 1.5); ylabel('r\_dot [rad/s^2]');
yyaxis right;
for i = 1:4
    plot(t_v, pwm_det(:,i), '--', 'Color', colors(i,:), 'LineWidth', 0.8); hold on;
end
ylabel('\DeltaPWM [us]'); xlabel('Tempo [s]'); grid on;
legend('r\_dot','C1','C2','C3','C4', 'Location', 'best');
title('Yaw: r\_dot vs \DeltaPWM');

sgtitle('Diagnostico de Mapeamento de Motores');
saveas(gcf, 'fig03_pwm_vs_derivadas.png');

% --- Figura: PWM e atitude no trecho VTOL ---
figure('Name', 'VTOL Overview', 'Position', [50 50 1200 800]);

subplot(3,1,1);
plot(t_v, pwm_v, 'LineWidth', 1);
legend('C1','C2','C3','C4', 'Location', 'best');
ylabel('PWM [us]'); grid on;
title('PWM dos Motores VTOL');

subplot(3,1,2);
plot(t_v, p_v, 'b', t_v, q_v, 'r', t_v, r_v, 'Color', [0.4 0.7 0.2], 'LineWidth', 1);
legend('p (roll)', 'q (pitch)', 'r (yaw)', 'Location', 'best');
ylabel('[rad/s]'); grid on;
title('Velocidades Angulares');

subplot(3,1,3);
roll_v  = roll(idx_vtol);
pitch_v = pitch(idx_vtol);
yaw_v   = yaw(idx_vtol);
plot(t_v, roll_v, 'b', t_v, pitch_v, 'r', t_v, yaw_v, 'Color', [0.4 0.7 0.2], 'LineWidth', 1);
legend('Roll', 'Pitch', 'Yaw', 'Location', 'best');
ylabel('[deg]'); xlabel('Tempo [s]'); grid on;
title('Atitude (EKF)');

sgtitle('Visao Geral do Trecho VTOL');
saveas(gcf, 'fig04_vtol_overview.png');

% --- Figura: Diagrama do frame com resultados ---
figure('Name', 'Frame Diagram', 'Position', [300 200 600 700]);
axis([-0.35 0.35 -0.45 0.42]);
hold on; axis equal; grid on;
title('Mapeamento: Canais PWM -> Posicao Fisica');
xlabel('Lateral [m] (+ = Direita)');
ylabel('Longitudinal [m] (+ = Frente/Nariz)');

% Geometria do frame
ly = 0.232;      % braço lateral [m]
lx_front = 0.311; % braço frontal [m]
lx_rear  = 0.343; % braço traseiro [m]

% Posicoes fisicas dos motores
motor_pos = [
    -ly,  lx_front;   % Frontal-Esquerdo
     ly,  lx_front;   % Frontal-Direito
    -ly, -lx_rear;    % Traseiro-Esquerdo
     ly, -lx_rear;    % Traseiro-Direito
];

motor_labels_phys = {'FRONTAL-ESQ', 'FRONTAL-DIR', 'TRASEIRO-ESQ', 'TRASEIRO-DIR'};

% Desenhar frame
plot([motor_pos(1,1), motor_pos(2,1)], [motor_pos(1,2), motor_pos(2,2)], 'k-', 'LineWidth', 2);
plot([motor_pos(3,1), motor_pos(4,1)], [motor_pos(3,2), motor_pos(4,2)], 'k-', 'LineWidth', 2);
plot([motor_pos(1,1), motor_pos(3,1)], [motor_pos(1,2), motor_pos(3,2)], 'k-', 'LineWidth', 2);
plot([motor_pos(2,1), motor_pos(4,1)], [motor_pos(2,2), motor_pos(4,2)], 'k-', 'LineWidth', 2);

% CG
plot(0, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3);
text(0.02, -0.02, 'CG', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');

% Nariz
annotation_arrow_y = lx_front + 0.05;
plot(0, annotation_arrow_y, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
text(0.02, annotation_arrow_y, 'NARIZ', 'FontSize', 10, 'FontWeight', 'bold');

% Plotar motores com canal correspondente
for m = 1:4
    % Encontrar qual canal corresponde a esta posicao
    ch_idx = 0;
    for c = 1:4
        if contains(pos_labels{c}, 'FRONTAL') && contains(motor_labels_phys{m}, 'FRONTAL') && ...
           contains(pos_labels{c}, 'ESQ') && contains(motor_labels_phys{m}, 'ESQ')
            ch_idx = c; break;
        elseif contains(pos_labels{c}, 'FRONTAL') && contains(motor_labels_phys{m}, 'FRONTAL') && ...
               contains(pos_labels{c}, 'DIR') && contains(motor_labels_phys{m}, 'DIR')
            ch_idx = c; break;
        elseif contains(pos_labels{c}, 'TRASEIRO') && contains(motor_labels_phys{m}, 'TRASEIRO') && ...
               contains(pos_labels{c}, 'ESQ') && contains(motor_labels_phys{m}, 'ESQ')
            ch_idx = c; break;
        elseif contains(pos_labels{c}, 'TRASEIRO') && contains(motor_labels_phys{m}, 'TRASEIRO') && ...
               contains(pos_labels{c}, 'DIR') && contains(motor_labels_phys{m}, 'DIR')
            ch_idx = c; break;
        end
    end

    x = motor_pos(m, 1);
    y = motor_pos(m, 2);
    plot(x, y, 'ko', 'MarkerSize', 20, 'MarkerFaceColor', colors(ch_idx,:), 'LineWidth', 2);

    if ch_idx > 0
        text(x, y, sprintf('C%d', ch_idx), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'w');
        text(x, y - 0.04, motor_labels_phys{m}, ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
    end
end

% Dimensoes
plot([0 ly], [lx_front+0.02 lx_front+0.02], 'b-', 'LineWidth', 1);
text(ly/2, lx_front+0.04, '232mm', 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'b');
plot([-ly-0.02 -ly-0.02], [0 lx_front], 'Color', [0.6 0 0], 'LineWidth', 1);
text(-ly-0.06, lx_front/2, '311mm', 'HorizontalAlignment', 'center', 'FontSize', 9, ...
    'Color', [0.6 0 0], 'Rotation', 90);
plot([-ly-0.02 -ly-0.02], [-lx_rear 0], 'Color', [0 0.5 0], 'LineWidth', 1);
text(-ly-0.06, -lx_rear/2, '343mm', 'HorizontalAlignment', 'center', 'FontSize', 9, ...
    'Color', [0 0.5 0], 'Rotation', 90);

saveas(gcf, 'fig05_frame_diagram.png');

%% ========================================================================
%  9. RESUMO FINAL
% =========================================================================
fprintf('\n\n');
fprintf('==============================================================\n');
fprintf('  RESUMO: MAPEAMENTO DE CANAIS -> POSICAO FISICA\n');
fprintf('==============================================================\n\n');
fprintf('             NARIZ (frente)\n');
fprintf('                |\n');

% Montar diagrama ASCII com os canais nas posicoes corretas
frame_pos = {'FRONTAL-ESQ', 'FRONTAL-DIR', 'TRASEIRO-ESQ', 'TRASEIRO-DIR'};
ch_at_pos = {'??', '??', '??', '??'};
for c = 1:4
    for fp = 1:4
        if contains(pos_labels{c}, strsplit(frame_pos{fp}, '-'))
            ch_at_pos{fp} = sprintf('C%d', c);
        end
    end
end

fprintf('       %s o----+----o %s      <- 311.18 mm\n', ch_at_pos{1}, ch_at_pos{2});
fprintf('          |    |    |\n');
fprintf('    232mm |    CG   | 232mm\n');
fprintf('          |         |\n');
fprintf('       %s o---------o %s      <- 342.87 mm\n', ch_at_pos{3}, ch_at_pos{4});
fprintf('\n');

for c = 1:4
    fprintf('  %s = %s\n', labels_ch{c}, pos_labels{c});
end

fprintf('\n  Figuras salvas em: %s\n', pwd);
fprintf('==============================================================\n');

fprintf('\nScript finalizado com sucesso.\n');
