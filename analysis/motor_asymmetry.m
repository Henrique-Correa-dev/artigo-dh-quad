%% motor_asymmetry.m
% =========================================================================
% Analise de Assimetria dos Motores — k_T e k_Q individuais
% =========================================================================
%
% Ideia: O modelo de bancada fornece T_ref(pwm) e Q_ref(pwm) para UM motor.
%        Na pratica, cada motor pode diferir (variacao de fabricacao,
%        helice, instalacao, bateria, etc.). Modelamos isso com fatores de
%        escala individuais:
%
%          T_i = k_Ti * T_ref(pwm_i)
%          Q_i = k_Qi * Q_ref(pwm_i)
%
%        Se todos fossem iguais: k_Ti = k_Qi = 1.0
%
% Metodo: Regressao linear (Equation Error Method)
%
%   Das equacoes de Newton (rotacional, simplificado perto de hover):
%
%     p_dot ≈ (1/Jx) * tau_x  = (1/Jx) * sum(-yi * k_Ti * Tref_i)
%     q_dot ≈ (1/Jy) * tau_y  = (1/Jy) * sum(+xi * k_Ti * Tref_i)
%     r_dot ≈ (1/Jz) * tau_z  = (1/Jz) * sum(di  * k_Qi * Qref_i)
%
%   Definindo parametros combinados (nao separamos k de J):
%
%     Roll:  p_dot = a1*(-Tref1) + a2*(+Tref2) + a3*(+Tref3) + a4*(-Tref4) + Bp
%            onde ai = ly * k_Ti / Jx
%
%     Pitch: q_dot = b1*(+Tref1) + b2*(-Tref2) + b3*(+Tref3) + b4*(-Tref4) + Bq
%            onde b1,b3 = lxf * k_Ti / Jy,  b2,b4 = lxr * k_Ti / Jy
%
%     Yaw:   r_dot = c1*(+Qref1) + c2*(+Qref2) + c3*(-Qref3) + c4*(-Qref4) + Br
%            onde ci = k_Qi / Jz
%
%   Cada equacao eh uma regressao linear -> solucao por minimos quadrados.
%   Os coeficientes relativos revelam a assimetria entre motores.
%
% Autor: Henrique / Claude
% Data: 2026-03-15
% =========================================================================

clear; close all; clc;

%% ========================================================================
%  1. DADOS DE BANCADA
% =========================================================================
pwm_bench = [1000; 1200; 1400; 1600; 1800; 2000];
thrust_g  = [   0;  143;  328;  532;  784;  843];
torque_Nm = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];

g_acc = 9.80665;
thrust_N = thrust_g * g_acc / 1000;

% Polinomios grau 3
coeffs_T = polyfit(pwm_bench, thrust_N, 3);
coeffs_Q = polyfit(pwm_bench, torque_Nm, 3);

idx_T0 = find(thrust_N > 1e-9, 1, 'first');
idx_Q0 = find(torque_Nm > 1e-9, 1, 'first');
pwm_min_T = pwm_bench(idx_T0);
pwm_min_Q = pwm_bench(idx_Q0);

T_ref = @(pwm) (pwm >= pwm_min_T) .* max(0, polyval(coeffs_T, pwm));
Q_ref = @(pwm) (pwm >= pwm_min_Q) .* max(0, polyval(coeffs_Q, pwm));

%% ========================================================================
%  2. GEOMETRIA
% =========================================================================
lx_f = 0.31118;  lx_r = 0.34287;  ly = 0.232;

% Posicao e sentido de cada canal (C1..C4)
x_m = [+lx_f; -lx_r; +lx_f; -lx_r];  % longitudinal
y_m = [+ly;   -ly;   -ly;   +ly  ];   % lateral
d_m = [+1;    +1;    -1;    -1   ];   % CW=+1, CCW=-1

%% ========================================================================
%  3. CARREGAR E INTERPOLAR DADOS DE VOO
% =========================================================================
fprintf('=== Carregando dados ===\n');
load(fullfile('..', 'identification', 'log_data.mat'));

IMU.TimeS  = double(IMU.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;

idx_imu = IMU.I == 0;
time_IMU  = IMU.TimeS(idx_imu);
time_RCOU = RCOU.TimeS;

dt = 0.1;
t_s = max(min(time_IMU), min(time_RCOU));
t_e = min(max(time_IMU), max(time_RCOU));
t = (t_s:dt:t_e)';

p = interp1(time_IMU, IMU.GyrX(idx_imu), t, 'linear');
q = interp1(time_IMU, IMU.GyrY(idx_imu), t, 'linear');
r = interp1(time_IMU, IMU.GyrZ(idx_imu), t, 'linear');

pwm = zeros(length(t), 4);
pwm(:,1) = interp1(time_RCOU, double(RCOU.C1), t, 'linear');
pwm(:,2) = interp1(time_RCOU, double(RCOU.C2), t, 'linear');
pwm(:,3) = interp1(time_RCOU, double(RCOU.C3), t, 'linear');
pwm(:,4) = interp1(time_RCOU, double(RCOU.C4), t, 'linear');

%% ========================================================================
%  4. SELECIONAR TRECHO DE ANALISE
% =========================================================================
t_ini = 150; t_fim = 200;
idx = (t >= t_ini) & (t <= t_fim);

t_v   = t(idx);
p_v   = p(idx);
q_v   = q(idx);
r_v   = r(idx);
pwm_v = pwm(idx,:);
N = sum(idx);

fprintf('  Trecho: %d-%ds (%d amostras)\n', t_ini, t_fim, N);

% Empuxo e torque de referencia (k=1) para cada motor
Tref = zeros(N, 4);
Qref = zeros(N, 4);
for i = 1:4
    Tref(:,i) = T_ref(pwm_v(:,i));
    Qref(:,i) = Q_ref(pwm_v(:,i));
end

% Derivadas angulares suavizadas
sw = 5;
p_dot = gradient(movmean(p_v, sw), dt);
q_dot = gradient(movmean(q_v, sw), dt);
r_dot = gradient(movmean(r_v, sw), dt);

%% ========================================================================
%  5. REGRESSAO LINEAR — ROLL (p_dot)
% =========================================================================
%  p_dot = a1*(-Tref1) + a2*(+Tref2) + a3*(+Tref3) + a4*(-Tref4) + Bp
%  onde ai = ly * k_Ti / Jx
%
%  Matriz de regressao: cada coluna = sinal * Tref_i

sign_roll = [-1, +1, +1, -1];  % sinais de -yi/ly

A_roll = Tref .* sign_roll;      % N x 4
A_roll = [A_roll, ones(N,1)];    % adicionar bias

theta_roll = A_roll \ p_dot;     % minimos quadrados

a_roll = theta_roll(1:4);  % ai = ly * k_Ti / Jx
Bp     = theta_roll(5);

p_dot_pred = A_roll * theta_roll;
R2_roll = 1 - sum((p_dot - p_dot_pred).^2) / sum((p_dot - mean(p_dot)).^2);

% k_T relativo (normalizado pela media)
kT_from_roll = a_roll / mean(a_roll);

fprintf('\n=== ROLL (p_dot) ===\n');
fprintf('  R² = %.4f\n', R2_roll);
fprintf('  Bias Bp = %.4f rad/s²\n', Bp);
fprintf('  Coeficientes ai (= ly*k_Ti/Jx):\n');
for i = 1:4
    fprintf('    a%d = %.6f  ->  k_T%d relativo = %.3f\n', i, a_roll(i), i, kT_from_roll(i));
end

%% ========================================================================
%  6. REGRESSAO LINEAR — PITCH (q_dot)
% =========================================================================
%  q_dot = b1*(+Tref1) + b2*(-Tref2) + b3*(+Tref3) + b4*(-Tref4) + Bq
%  onde b1,b3 = lxf*k_Ti/Jy  e  b2,b4 = lxr*k_Ti/Jy

sign_pitch = [+1, -1, +1, -1];  % sinais de +xi / |xi|

A_pitch = Tref .* sign_pitch;
A_pitch = [A_pitch, ones(N,1)];

theta_pitch = A_pitch \ q_dot;

b_pitch = theta_pitch(1:4);
Bq      = theta_pitch(5);

q_dot_pred = A_pitch * theta_pitch;
R2_pitch = 1 - sum((q_dot - q_dot_pred).^2) / sum((q_dot - mean(q_dot)).^2);

% Para pitch, os bracos sao diferentes (lxf vs lxr),
% entao o coeficiente bi ja inclui o braco:
%   b1 = lxf*k_T1/Jy,  b2 = lxr*k_T2/Jy
% Para obter k_T relativo, normalizamos pelo braco:
kT_from_pitch_raw = [b_pitch(1)/lx_f; b_pitch(2)/lx_r; b_pitch(3)/lx_f; b_pitch(4)/lx_r];
kT_from_pitch = kT_from_pitch_raw / mean(kT_from_pitch_raw);

fprintf('\n=== PITCH (q_dot) ===\n');
fprintf('  R² = %.4f\n', R2_pitch);
fprintf('  Bias Bq = %.4f rad/s²\n', Bq);
fprintf('  Coeficientes bi:\n');
for i = 1:4
    if i == 1 || i == 3, arm_str = 'lxf'; else, arm_str = 'lxr'; end
    fprintf('    b%d = %.6f  (%s*k_T%d/Jy)  ->  k_T%d relativo = %.3f\n', ...
        i, b_pitch(i), arm_str, i, i, kT_from_pitch(i));
end

%% ========================================================================
%  7. REGRESSAO LINEAR — YAW (r_dot)
% =========================================================================
%  r_dot = c1*(+Qref1) + c2*(+Qref2) + c3*(-Qref3) + c4*(-Qref4) + Br
%  onde ci = k_Qi / Jz

sign_yaw = [+1, +1, -1, -1];  % di

A_yaw = Qref .* sign_yaw;
A_yaw = [A_yaw, ones(N,1)];

theta_yaw = A_yaw \ r_dot;

c_yaw = theta_yaw(1:4);
Br    = theta_yaw(5);

r_dot_pred = A_yaw * theta_yaw;
R2_yaw = 1 - sum((r_dot - r_dot_pred).^2) / sum((r_dot - mean(r_dot)).^2);

kQ_from_yaw = c_yaw / mean(c_yaw);

fprintf('\n=== YAW (r_dot) ===\n');
fprintf('  R² = %.4f\n', R2_yaw);
fprintf('  Bias Br = %.4f rad/s²\n', Br);
fprintf('  Coeficientes ci (= k_Qi/Jz):\n');
for i = 1:4
    fprintf('    c%d = %.6f  ->  k_Q%d relativo = %.3f\n', i, c_yaw(i), i, kQ_from_yaw(i));
end

%% ========================================================================
%  8. RESUMO COMPARATIVO
% =========================================================================
fprintf('\n');
fprintf('==============================================================\n');
fprintf('  RESUMO: FATORES k_T e k_Q RELATIVOS (media = 1.0)\n');
fprintf('==============================================================\n');
fprintf('\n');
fprintf('  Motor    Posicao           k_T(roll)  k_T(pitch)  k_Q(yaw)\n');
fprintf('  %-6s   %-18s %+9.3f   %+9.3f   %+9.3f\n', ...
    'C1', 'Frontal-Dir (CW)',  kT_from_roll(1), kT_from_pitch(1), kQ_from_yaw(1));
fprintf('  %-6s   %-18s %+9.3f   %+9.3f   %+9.3f\n', ...
    'C2', 'Traseiro-Esq (CW)', kT_from_roll(2), kT_from_pitch(2), kQ_from_yaw(2));
fprintf('  %-6s   %-18s %+9.3f   %+9.3f   %+9.3f\n', ...
    'C3', 'Frontal-Esq (CCW)', kT_from_roll(3), kT_from_pitch(3), kQ_from_yaw(3));
fprintf('  %-6s   %-18s %+9.3f   %+9.3f   %+9.3f\n', ...
    'C4', 'Traseiro-Dir (CCW)', kT_from_roll(4), kT_from_pitch(4), kQ_from_yaw(4));

% Media entre roll e pitch para k_T
kT_combined = (kT_from_roll + kT_from_pitch) / 2;
fprintf('\n  k_T combinado (media roll+pitch):\n');
for i = 1:4
    fprintf('    k_T%d = %.3f\n', i, kT_combined(i));
end
fprintf('    k_Q  = [%.3f, %.3f, %.3f, %.3f]\n', kQ_from_yaw);

%% ========================================================================
%  9. FIGURAS
% =========================================================================

% --- Bar chart dos k relativos ---
figure('Name', 'k_T e k_Q Relativos', 'Position', [100 100 900 450], 'Color', 'w');

subplot(1,2,1);
bar_data_kT = [kT_from_roll(:), kT_from_pitch(:), kT_combined(:)];
b = bar(bar_data_kT);
b(1).FaceColor = [0.2 0.5 0.8];
b(2).FaceColor = [0.8 0.3 0.3];
b(3).FaceColor = [0.3 0.7 0.3];
set(gca, 'XTickLabel', {'C1 (F-Dir)', 'C2 (T-Esq)', 'C3 (F-Esq)', 'C4 (T-Dir)'});
ylabel('k_T relativo'); grid on;
yline(1, 'k--', 'Ideal', 'LineWidth', 1.2);
legend('via Roll', 'via Pitch', 'Combinado', 'Location', 'best');
title('Fator de Escala de Empuxo k_T');

subplot(1,2,2);
bar(kQ_from_yaw, 'FaceColor', [0.6 0.3 0.7]);
set(gca, 'XTickLabel', {'C1 (F-Dir)', 'C2 (T-Esq)', 'C3 (F-Esq)', 'C4 (T-Dir)'});
ylabel('k_Q relativo'); grid on;
yline(1, 'k--', 'Ideal', 'LineWidth', 1.2);
title('Fator de Escala de Torque k_Q');

saveas(gcf, 'fig09_kT_kQ_relativos.png');

% --- Predicao vs medido ---
figure('Name', 'Predicao vs Medido', 'Position', [50 50 1200 700], 'Color', 'w');

subplot(3,1,1);
plot(t_v, p_dot, 'b-', 'LineWidth', 0.8); hold on;
plot(t_v, p_dot_pred, 'r--', 'LineWidth', 1.2);
ylabel('p_{dot} [rad/s^2]'); grid on;
legend('Medido', sprintf('Predito (R^2=%.3f)', R2_roll), 'Location', 'best');
title('Roll: Regressao com k_{Ti} individuais');

subplot(3,1,2);
plot(t_v, q_dot, 'b-', 'LineWidth', 0.8); hold on;
plot(t_v, q_dot_pred, 'r--', 'LineWidth', 1.2);
ylabel('q_{dot} [rad/s^2]'); grid on;
legend('Medido', sprintf('Predito (R^2=%.3f)', R2_pitch), 'Location', 'best');
title('Pitch: Regressao com k_{Ti} individuais');

subplot(3,1,3);
plot(t_v, r_dot, 'b-', 'LineWidth', 0.8); hold on;
plot(t_v, r_dot_pred, 'r--', 'LineWidth', 1.2);
ylabel('r_{dot} [rad/s^2]'); xlabel('Tempo [s]'); grid on;
legend('Medido', sprintf('Predito (R^2=%.3f)', R2_yaw), 'Location', 'best');
title('Yaw: Regressao com k_{Qi} individuais');

sgtitle('Equacao de Erro: Derivadas preditas vs medidas (k individuais)');
saveas(gcf, 'fig10_predicao_kT_kQ.png');

% --- Comparacao: k=1 (ideal) vs k individuais ---
% Recalcular predicao com k=1 (sem regressao)
% Precisamos de J, que nao conhecemos. Usamos a regressao com k=1 forcado.
A_roll_k1 = [sum(Tref .* sign_roll, 2), ones(N,1)];
theta_roll_k1 = A_roll_k1 \ p_dot;
p_dot_k1 = A_roll_k1 * theta_roll_k1;
R2_roll_k1 = 1 - sum((p_dot - p_dot_k1).^2) / sum((p_dot - mean(p_dot)).^2);

A_pitch_k1 = [sum(Tref .* sign_pitch .* [lx_f, lx_r, lx_f, lx_r], 2), ones(N,1)];
theta_pitch_k1 = A_pitch_k1 \ q_dot;
q_dot_k1 = A_pitch_k1 * theta_pitch_k1;
R2_pitch_k1 = 1 - sum((q_dot - q_dot_k1).^2) / sum((q_dot - mean(q_dot)).^2);

A_yaw_k1 = [sum(Qref .* sign_yaw, 2), ones(N,1)];
theta_yaw_k1 = A_yaw_k1 \ r_dot;
r_dot_k1 = A_yaw_k1 * theta_yaw_k1;
R2_yaw_k1 = 1 - sum((r_dot - r_dot_k1).^2) / sum((r_dot - mean(r_dot)).^2);

fprintf('\n');
fprintf('==============================================================\n');
fprintf('  COMPARACAO: k=1 (identicos) vs k individuais\n');
fprintf('==============================================================\n');
fprintf('\n');
fprintf('  Eixo     R²(k=1)   R²(k_ind)   Melhora\n');
fprintf('  Roll    %7.4f    %7.4f     %+.4f\n', R2_roll_k1, R2_roll, R2_roll - R2_roll_k1);
fprintf('  Pitch   %7.4f    %7.4f     %+.4f\n', R2_pitch_k1, R2_pitch, R2_pitch - R2_pitch_k1);
fprintf('  Yaw     %7.4f    %7.4f     %+.4f\n', R2_yaw_k1, R2_yaw, R2_yaw - R2_yaw_k1);
fprintf('\n');

fprintf('  Se a melhora eh significativa, motores individuais diferem.\n');
fprintf('  Se a melhora eh pequena, k=1 para todos eh aceitavel.\n');
fprintf('\n==============================================================\n');
fprintf('\nScript finalizado.\n');
