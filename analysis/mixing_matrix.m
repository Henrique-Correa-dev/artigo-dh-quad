%% mixing_matrix.m
% =========================================================================
% Derivacao da Mixing Matrix — DH Hybrid-Drone VTOL (H-frame)
% =========================================================================
%
% Objetivo: Partindo dos dados de bancada e da geometria do frame,
%           derivar a mixing matrix que mapeia PWM dos 4 motores para
%           as forcas e momentos atuando no drone.
%
% Cadeia completa:
%
%   PWM(us) --[bancada]--> T(N), Q(Nm) --[geometria]--> [T, tx, ty, tz]
%
%   ┌  T  ┐       ┌ T1(pwm1) ┐         ┌ pwm1 ┐
%   │ tx  │ = M · │ T2(pwm2) │   <--   │ pwm2 │  via bancada
%   │ ty  │       │ T3(pwm3) │         │ pwm3 │
%   └ tz  ┘       └ T4(pwm4) ┘         └ pwm4 ┘
%
% Convencoes:
%   - Frame NED (x=frente, y=direita, z=baixo)
%   - Empuxo aponta em -k_b (para cima)
%   - Roll positivo  (tx > 0): asa direita desce
%   - Pitch positivo (ty > 0): nariz sobe
%   - Yaw positivo   (tz > 0): nariz vira para direita (CW visto de cima)
%
% Mapeamento (resultado da analise de correlacao):
%   C1 = Frontal-Direito  (CW)
%   C2 = Traseiro-Esquerdo (CW)
%   C3 = Frontal-Esquerdo  (CCW)
%   C4 = Traseiro-Direito  (CCW)
%
% Autor: Henrique / Claude
% Data: 2026-03-15
% =========================================================================

clear; close all; clc;

%% ========================================================================
%  1. DADOS DE BANCADA (dinamometro)
% =========================================================================
% Teste realizado com o mesmo motor/helice usado no voo.
% PWM de 1000 a 2000 us (0% a 100% throttle).

pwm_bench     = [1000; 1200; 1400; 1600; 1800; 2000];  % [us]
thrust_g      = [   0;  143;  328;  532;  784;  843];   % [gramas-forca]
torque_Nm     = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176]; % [Nm]

% Converter empuxo: gramas-forca -> Newtons
g = 9.80665;  % [m/s^2]
thrust_N = thrust_g * g / 1000;

fprintf('=== DADOS DE BANCADA ===\n');
fprintf('  PWM(us)  Empuxo(g)  Empuxo(N)  Torque(Nm)\n');
fprintf('  %6.0f   %7.0f    %7.3f    %8.4f\n', ...
    [pwm_bench, thrust_g, thrust_N, torque_Nm]');

%% ========================================================================
%  2. MODELOS POLINOMIAIS: PWM -> T(N) e PWM -> Q(Nm)
% =========================================================================
% Ajuste polinomial de grau 3 (cubico) aos dados de bancada.
% Estes polinomios representam a curva de referencia do motor.

poly_degree = 3;

% Ajustar polinomios
coeffs_T = polyfit(pwm_bench, thrust_N, poly_degree);
coeffs_Q = polyfit(pwm_bench, torque_Nm, poly_degree);

% PWM minimo para empuxo/torque > 0 (dead zone)
idx_T_active = find(thrust_N > 1e-9, 1, 'first');
idx_Q_active = find(torque_Nm > 1e-9, 1, 'first');
pwm_min_T = pwm_bench(idx_T_active);
pwm_min_Q = pwm_bench(idx_Q_active);

% Funcoes do modelo (com dead zone e clamp >= 0)
T_ref = @(pwm) (pwm >= pwm_min_T) .* max(0, polyval(coeffs_T, pwm));
Q_ref = @(pwm) (pwm >= pwm_min_Q) .* max(0, polyval(coeffs_Q, pwm));

fprintf('\n=== MODELOS POLINOMIAIS (grau %d) ===\n', poly_degree);
fprintf('  T_ref(pwm) = %.6e*pwm^3 + %.6e*pwm^2 + %.6e*pwm + %.6e  [N]\n', coeffs_T);
fprintf('  Q_ref(pwm) = %.6e*pwm^3 + %.6e*pwm^2 + %.6e*pwm + %.6e  [Nm]\n', coeffs_Q);
fprintf('  Dead zone: T=0 para PWM < %.0f, Q=0 para PWM < %.0f\n', pwm_min_T, pwm_min_Q);

% --- Figura: ajuste dos polinomios ---
pwm_fine = linspace(1000, 2000, 200);

figure('Name', 'Modelo de Bancada', 'Position', [100 100 900 380], 'Color', 'w');

subplot(1,2,1);
plot(pwm_fine, T_ref(pwm_fine), 'b-', 'LineWidth', 1.5); hold on;
plot(pwm_bench, thrust_N, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('PWM (\mus)'); ylabel('Empuxo (N)');
title('Empuxo vs PWM'); grid on;
legend('Polinomio 3^o grau', 'Dados bancada', 'Location', 'northwest');
xlim([950 2050]);

subplot(1,2,2);
plot(pwm_fine, Q_ref(pwm_fine), 'b-', 'LineWidth', 1.5); hold on;
plot(pwm_bench, torque_Nm, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('PWM (\mus)'); ylabel('Torque (N\cdotm)');
title('Torque Reativo vs PWM'); grid on;
legend('Polinomio 3^o grau', 'Dados bancada', 'Location', 'northwest');
xlim([950 2050]);

saveas(gcf, 'fig06_modelo_bancada.png');

%% ========================================================================
%  3. GEOMETRIA DO FRAME H
% =========================================================================
%
%  Mapeamento confirmado (motor_mapping_analysis.m):
%
%              NARIZ (frente, +x)
%                  |
%     C3(CCW) o----+----o C1(CW)       <- lx_f = 311.18 mm
%              |   CG    |
%     C2(CW)  o---------o C4(CCW)      <- lx_r = 342.87 mm
%           ly              ly
%         232mm            232mm
%
%  Posicoes dos motores (x,y) relativas ao CG:
%    C1: (+lx_f, +ly)   Frontal-Direito   CW
%    C2: (-lx_r, -ly)   Traseiro-Esquerdo CW
%    C3: (+lx_f, -ly)   Frontal-Esquerdo  CCW
%    C4: (-lx_r, +ly)   Traseiro-Direito  CCW

lx_f = 0.31118;   % braco frontal [m]
lx_r = 0.34287;   % braco traseiro [m]
ly   = 0.232;     % braco lateral [m]

% Posicao (x, y) de cada canal relativo ao CG
%         x        y      sentido (CW=+1, CCW=-1)
motors = [
    +lx_f,  +ly,   +1;    % C1: Frontal-Direito,  CW
    -lx_r,  -ly,   +1;    % C2: Traseiro-Esquerdo, CW
    +lx_f,  -ly,   -1;    % C3: Frontal-Esquerdo,  CCW
    -lx_r,  +ly,   -1;    % C4: Traseiro-Direito,  CCW
];

x_m = motors(:,1);   % posicao longitudinal [m]
y_m = motors(:,2);   % posicao lateral [m]
d_m = motors(:,3);   % direcao yaw: +1=CW, -1=CCW

fprintf('\n=== GEOMETRIA DO FRAME ===\n');
fprintf('  Braco frontal: %.3f m\n', lx_f);
fprintf('  Braco traseiro: %.3f m\n', lx_r);
fprintf('  Braco lateral: %.3f m\n', ly);
fprintf('\n  Motor   Posicao (x,y) [m]      Rotacao\n');
for i = 1:4
    if d_m(i) > 0, dir_str = 'CW'; else, dir_str = 'CCW'; end
    fprintf('  C%d     (%+.3f, %+.3f)        %s\n', i, x_m(i), y_m(i), dir_str);
end

%% ========================================================================
%  4. DERIVACAO DA MIXING MATRIX
% =========================================================================
%
%  Para um motor na posicao (xi, yi) com empuxo Ti ao longo de -k_b:
%
%    Forca do motor:  F_i = (0, 0, -Ti)
%    Posicao:         r_i = (xi, yi, 0)
%
%    Momento = r_i x F_i = | i     j     k  |
%                          | xi    yi    0  |
%                          | 0     0    -Ti |
%
%    tau_i = (yi*(-Ti) - 0, 0 - xi*(-Ti), 0)
%          = (-yi*Ti, +xi*Ti, 0)
%
%  Portanto:
%    tau_x (roll)  = -yi * Ti     <- braco lateral
%    tau_y (pitch) = +xi * Ti     <- braco longitudinal
%    tau_z (yaw)   = di * Qi      <- torque reativo (CW=+, CCW=-)
%
%  Somando os 4 motores:
%
%    ┌  T  ┐   ┌   1     1     1     1  ┐ ┌ T1 ┐   ┌ 0  0  0  0 ┐ ┌ Q1 ┐
%    │ tx  │ = │ -y1   -y2   -y3   -y4  │ │ T2 │ + │ 0  0  0  0 │ │ Q2 │
%    │ ty  │   │ +x1   +x2   +x3   +x4  │ │ T3 │   │ 0  0  0  0 │ │ Q3 │
%    └ tz  ┘   └   0     0     0     0  ┘ └ T4 ┘   └ d1 d2 d3 d4 ┘ └ Q4 ┘
%
%  Ou de forma compacta:  u = M_T * T_vec + M_Q * Q_vec

% Matriz de empuxo
M_T = [  ones(1,4);       % T  = T1 + T2 + T3 + T4
        -y_m';             % tx = sum(-yi * Ti)
        +x_m';             % ty = sum(+xi * Ti)
         zeros(1,4)];      % tz (contribuicao de empuxo = 0)

% Matriz de torque reativo
M_Q = [ zeros(1,4);        % T  (sem contribuicao)
        zeros(1,4);         % tx (sem contribuicao)
        zeros(1,4);         % ty (sem contribuicao)
        d_m'];              % tz = sum(di * Qi)

% Mixing matrix combinada (em termos de Ti e Qi separados)
fprintf('\n=== MIXING MATRIX (empuxo) ===\n');
fprintf('  ┌  T  ┐   ┌ %+7.3f  %+7.3f  %+7.3f  %+7.3f ┐ ┌ T1 ┐\n', M_T(1,:));
fprintf('  │ tx  │ = │ %+7.3f  %+7.3f  %+7.3f  %+7.3f │ │ T2 │\n', M_T(2,:));
fprintf('  │ ty  │   │ %+7.3f  %+7.3f  %+7.3f  %+7.3f │ │ T3 │\n', M_T(3,:));
fprintf('  └ tz  ┘   └ %+7.3f  %+7.3f  %+7.3f  %+7.3f ┘ └ T4 ┘\n', M_T(4,:));

fprintf('\n=== MIXING MATRIX (torque reativo) ===\n');
fprintf('  tz = %+.0f*Q1 %+.0f*Q2 %+.0f*Q3 %+.0f*Q4\n', d_m');

fprintf('\n=== EQUACOES EXPANDIDAS ===\n');
fprintf('  T  =  T1 + T2 + T3 + T4\n');
fprintf('  tx = %+.3f*T1 %+.3f*T2 %+.3f*T3 %+.3f*T4\n', -y_m');
fprintf('  ty = %+.3f*T1 %+.3f*T2 %+.3f*T3 %+.3f*T4\n', +x_m');
fprintf('  tz = %+.0f*Q1 %+.0f*Q2 %+.0f*Q3 %+.0f*Q4\n', d_m');

%% ========================================================================
%  5. VERIFICACAO COM DADOS DE VOO
% =========================================================================
% Carregar o log e calcular forcas/momentos para um trecho de voo,
% depois comparar com as derivadas angulares medidas.

fprintf('\n=== VERIFICACAO COM DADOS DE VOO ===\n');
load(fullfile('..', 'identification', 'log_data.mat'));

% Interpolar para grade comum
ATT.TimeS  = double(ATT.TimeUS) / 1e6;
IMU.TimeS  = double(IMU.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;

idx_imu = IMU.I == 0;
time_IMU  = IMU.TimeS(idx_imu);
time_RCOU = RCOU.TimeS;

dt = 0.1;
t_start = max(min(time_IMU), min(time_RCOU));
t_end   = min(max(time_IMU), max(time_RCOU));
t = (t_start:dt:t_end)';

p = interp1(time_IMU, IMU.GyrX(idx_imu), t, 'linear');
q = interp1(time_IMU, IMU.GyrY(idx_imu), t, 'linear');
r = interp1(time_IMU, IMU.GyrZ(idx_imu), t, 'linear');

pwm = zeros(length(t), 4);
pwm(:,1) = interp1(time_RCOU, double(RCOU.C1), t, 'linear');
pwm(:,2) = interp1(time_RCOU, double(RCOU.C2), t, 'linear');
pwm(:,3) = interp1(time_RCOU, double(RCOU.C3), t, 'linear');
pwm(:,4) = interp1(time_RCOU, double(RCOU.C4), t, 'linear');

% Trecho de analise: 150-200s
idx_v = (t >= 150) & (t <= 200);
t_v   = t(idx_v);
p_v   = p(idx_v);
q_v   = q(idx_v);
r_v   = r(idx_v);
pwm_v = pwm(idx_v, :);

% Calcular empuxo e torque de cada motor via modelo de bancada
T_motors = zeros(size(pwm_v));
Q_motors = zeros(size(pwm_v));
for i = 1:4
    T_motors(:,i) = T_ref(pwm_v(:,i));
    Q_motors(:,i) = Q_ref(pwm_v(:,i));
end

% Calcular forcas e momentos via mixing matrix
T_total = sum(T_motors, 2);
tau_x   = T_motors * (-y_m);    % roll
tau_y   = T_motors * (+x_m);    % pitch
tau_z   = Q_motors * d_m;       % yaw

% Derivadas angulares medidas (para comparacao)
sw = 5;
p_dot = gradient(movmean(p_v, sw), dt);
q_dot = gradient(movmean(q_v, sw), dt);
r_dot = gradient(movmean(r_v, sw), dt);

% Correlacao: momentos calculados vs derivadas medidas
corr_roll  = corr(tau_x, p_dot);
corr_pitch = corr(tau_y, q_dot);
corr_yaw   = corr(tau_z, r_dot);

fprintf('  Correlacao (momento calculado vs derivada medida):\n');
fprintf('    Roll:  corr(tau_x, p_dot) = %+.4f\n', corr_roll);
fprintf('    Pitch: corr(tau_y, q_dot) = %+.4f\n', corr_pitch);
fprintf('    Yaw:   corr(tau_z, r_dot) = %+.4f\n', corr_yaw);
fprintf('  (Valores positivos confirmam que a mixing matrix esta correta)\n');

% --- Figura: momentos calculados vs derivadas medidas ---
figure('Name', 'Verificacao Mixing Matrix', 'Position', [50 50 1200 800], 'Color', 'w');

subplot(4,1,1);
plot(t_v, T_total, 'k-', 'LineWidth', 1.2);
ylabel('T_{total} [N]'); grid on;
title('Empuxo Total (soma dos 4 motores)');

subplot(4,1,2);
yyaxis left;
plot(t_v, tau_x, 'b-', 'LineWidth', 1.2); ylabel('\tau_x [Nm]');
yyaxis right;
plot(t_v, p_dot, 'r--', 'LineWidth', 1); ylabel('p_{dot} [rad/s^2]');
title(sprintf('Roll: \\tau_x vs \\dot{p}  (corr = %+.3f)', corr_roll));
legend('\tau_x (calculado)', '\dot{p} (medido)', 'Location', 'best');
grid on;

subplot(4,1,3);
yyaxis left;
plot(t_v, tau_y, 'b-', 'LineWidth', 1.2); ylabel('\tau_y [Nm]');
yyaxis right;
plot(t_v, q_dot, 'r--', 'LineWidth', 1); ylabel('q_{dot} [rad/s^2]');
title(sprintf('Pitch: \\tau_y vs \\dot{q}  (corr = %+.3f)', corr_pitch));
legend('\tau_y (calculado)', '\dot{q} (medido)', 'Location', 'best');
grid on;

subplot(4,1,4);
yyaxis left;
plot(t_v, tau_z, 'b-', 'LineWidth', 1.2); ylabel('\tau_z [Nm]');
yyaxis right;
plot(t_v, r_dot, 'r--', 'LineWidth', 1); ylabel('r_{dot} [rad/s^2]');
title(sprintf('Yaw: \\tau_z vs \\dot{r}  (corr = %+.3f)', corr_yaw));
legend('\tau_z (calculado)', '\dot{r} (medido)', 'Location', 'best');
xlabel('Tempo [s]'); grid on;

sgtitle('Verificacao: Momentos (Mixing Matrix + Bancada) vs Derivadas Medidas');
saveas(gcf, 'fig07_verificacao_mixing.png');

% --- Figura: empuxo e torque individuais ---
figure('Name', 'Forcas por Motor', 'Position', [50 50 1200 500], 'Color', 'w');

subplot(1,2,1);
plot(t_v, T_motors, 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('Empuxo [N]');
title('Empuxo por Motor'); grid on;
legend('C1 (F-Dir)', 'C2 (T-Esq)', 'C3 (F-Esq)', 'C4 (T-Dir)', 'Location', 'best');

subplot(1,2,2);
plot(t_v, Q_motors, 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('Torque [Nm]');
title('Torque Reativo por Motor'); grid on;
legend('C1 (F-Dir)', 'C2 (T-Esq)', 'C3 (F-Esq)', 'C4 (T-Dir)', 'Location', 'best');

saveas(gcf, 'fig08_forcas_por_motor.png');

%% ========================================================================
%  6. RESUMO
% =========================================================================
fprintf('\n');
fprintf('==============================================================\n');
fprintf('  RESUMO DA MIXING MATRIX\n');
fprintf('==============================================================\n');
fprintf('\n');
fprintf('  ┌  T  ┐   ┌ %+7.3f  %+7.3f  %+7.3f  %+7.3f ┐ ┌ T1 ┐   ┌ 0  0  0  0 ┐ ┌ Q1 ┐\n', M_T(1,:));
fprintf('  │ tx  │ = │ %+7.3f  %+7.3f  %+7.3f  %+7.3f │ │ T2 │ + │ 0  0  0  0 │ │ Q2 │\n', M_T(2,:));
fprintf('  │ ty  │   │ %+7.3f  %+7.3f  %+7.3f  %+7.3f │ │ T3 │   │ 0  0  0  0 │ │ Q3 │\n', M_T(3,:));
fprintf('  └ tz  ┘   └ %+7.3f  %+7.3f  %+7.3f  %+7.3f ┘ └ T4 ┘   └%+2.0f %+2.0f %+2.0f %+2.0f ┘ └ Q4 ┘\n', M_T(4,:), d_m');
fprintf('\n');
fprintf('  Onde Ti = T_ref(pwm_i) e Qi = Q_ref(pwm_i)  [polinomios de bancada]\n');
fprintf('\n');
fprintf('  Expandido:\n');
fprintf('    T  = T1 + T2 + T3 + T4\n');
fprintf('    tx = -%.3f*(T1+T4) + %.3f*(T2+T3)      [roll]\n', ly, ly);
fprintf('    ty = +%.3f*(T1+T3) - %.3f*(T2+T4)      [pitch]\n', lx_f, lx_r);
fprintf('    tz = +(Q1+Q2) - (Q3+Q4)                  [yaw]\n');
fprintf('\n');
fprintf('==============================================================\n');
fprintf('\nScript finalizado.\n');
