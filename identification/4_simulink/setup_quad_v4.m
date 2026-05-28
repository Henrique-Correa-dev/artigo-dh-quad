%SETUP_QUAD_V4  Prepara workspace pra simular quad_model_v4.slx
%
% v4.slx ATUALIZADO (sem biases p/q/r, sem drag Xu/Yv/Zw, sem Bz, com Lx/Ly
% assimétricos). Espelha vtol_dynamics.m + accelerometer_model.m atuais.
%
% Variáveis populadas no base workspace:
%   pwm1_in..pwm4_in              (timeseries)
%   p_out, q_out, r_out           (timeseries de referência — IMU gyro)
%   phi_out, theta_out, psi_out   (timeseries de referência — EKF, GRAUS)
%   accX_out, accY_out, accZ_out  (timeseries de referência — IMU acc)
%   phi0, theta0, psi0            (escalares, ICs)
%   mass, g_acc, tau_motor        (escalares)
%   P_estimated (20×1)            (G-formulation, via P_J_to_simulink)
%   Lx_r, Lx_l, Ly_f, Ly_r        (braços de momento, via P_J_to_simulink)
%   r_imu                         (offset CG→IMU, via P_J_to_simulink)
%   bias_acc                      (bias DC do acelerômetro, via P_J_to_simulink)
%
% Janela e fonte de parâmetros configuráveis abaixo.

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % raiz pra setup_paths
paths = setup_paths();

% ===================================================================
%  ESCOLHA DA JANELA E FONTE DE PARÂMETROS
% ===================================================================
LOG_FILE = '4 25-05-2026 09-31-48.log-132954.mat';   % log a usar
t_window = [473, 483];                                % [t_start, t_end] em segundos (10s teste rápido)

% Fonte do P_J:
%   'p0'      → chute inicial (parameters().P0_J)
%   'p_final' → resultado da identificação (P_identified.mat)
%   'manual'  → P_MANUAL definido abaixo
P_SOURCE = 'p_final';

% P_MANUAL (15 elementos — usado se P_SOURCE = 'manual')
P_MANUAL = [ ...
    0.072; 0.200; 0.150; 0.000933;    ... % Jx Jy Jz Jxz
    0.71;  0.71;  0.66;  0.67;         ... % k_T1..4
    0.79;  0.79;  0.87;  0.86;         ... % k_Q1..4
    2.27;  1.44;  0.47];                  % Dp Dq Dr

% Dinâmica do motor (1ª ordem entre Polyval e K_T no v4):
%   T_motor(s) = 1/(tau_motor·s + 1)
%   Pra desabilitar (motor "instantâneo"): tau_motor = 1e-9
tau_motor = 1e-9;   % matches parameters().tau_motor = 0 — sem lag


%% ===================================================================
%  1) Carregar parâmetros e P_J
%  ===================================================================
proj_params = parameters();
mass  = proj_params.m;
g_acc = proj_params.g;

switch lower(P_SOURCE)
    case 'p0'
        P_J = proj_params.P0_J;
        fprintf('setup_quad_v4: usando P0 (chute inicial).\n');
    case 'p_final'
        P_file = fullfile(paths.outputs, 'P_identified.mat');
        if exist(P_file, 'file')
            Pdat = load(P_file);
            P_J = Pdat.P_final;
            fprintf('setup_quad_v4: usando P_final de P_identified.mat.\n');
        else
            warning('P_identified.mat não encontrado — fallback pra P0.');
            P_J = proj_params.P0_J;
        end
    case 'manual'
        P_J = P_MANUAL(:);
        fprintf('setup_quad_v4: usando P_MANUAL.\n');
    otherwise
        error('P_SOURCE inválido: %s', P_SOURCE);
end
if numel(P_J) ~= 15
    error('P_J deve ter 15 elementos (recebeu %d).', numel(P_J));
end


%% ===================================================================
%  2) Converter P_J (15) → P_estimated (20) e popular workspace
%      (Lx_r, Lx_l, Ly_f, Ly_r, r_imu, bias_acc também são exportados)
%  ===================================================================
P_estimated = P_J_to_simulink(P_J);


%% ===================================================================
%  3) Carregar log e reamostrar na grade comum (dt=0.1s)
%  ===================================================================
log_path = fullfile(paths.data, LOG_FILE);
L = load_log_data(log_path);

t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
t_common = (t_lo:0.1:t_hi)';

if t_window(1) < t_lo || t_window(2) > t_hi
    error('Janela [%g, %g] fora do log [%g, %g].', ...
        t_window(1), t_window(2), t_lo, t_hi);
end

% Recortar janela
idx_w = (t_common >= t_window(1)) & (t_common <= t_window(2));
time_w   = t_common(idx_w);
time_rel = time_w - time_w(1);   % v4 começa em t=0
t_sim    = time_rel(end);

% Reamostrar sinais na grade comum, recortar janela
pwm1 = interp1(L.time_RCOU, L.pwm1_raw, t_common, 'linear');
pwm2 = interp1(L.time_RCOU, L.pwm2_raw, t_common, 'linear');
pwm3 = interp1(L.time_RCOU, L.pwm3_raw, t_common, 'linear');
pwm4 = interp1(L.time_RCOU, L.pwm4_raw, t_common, 'linear');
gyrX = interp1(L.time_IMU, L.gyrX_raw, t_common, 'linear');
gyrY = interp1(L.time_IMU, L.gyrY_raw, t_common, 'linear');
gyrZ = interp1(L.time_IMU, L.gyrZ_raw, t_common, 'linear');
accX = interp1(L.time_IMU, L.accX_raw, t_common, 'linear');
accY = interp1(L.time_IMU, L.accY_raw, t_common, 'linear');
accZ = interp1(L.time_IMU, L.accZ_raw, t_common, 'linear');
roll_deg  = interp1(L.time_ATT, L.roll_deg,  t_common, 'linear');
pitch_deg = interp1(L.time_ATT, L.pitch_deg, t_common, 'linear');
yaw_deg   = interp1(L.time_ATT, L.yaw_deg,   t_common, 'linear');


%% ===================================================================
%  4) Condições iniciais dos State-Space do lag de motor
%      O estado é a saída do Polyval CRUA (antes do K_T) no instante 0
%  ===================================================================
[func_T_ref, func_Q_ref] = motor_models();   % spline igual ao .m e ao .slx
pwm0_each = [pwm1(find(idx_w,1,'first')), pwm2(find(idx_w,1,'first')), ...
             pwm3(find(idx_w,1,'first')), pwm4(find(idx_w,1,'first'))];
T_eff_init = zeros(4,1);
Q_eff_init = zeros(4,1);
for i = 1:4
    T_eff_init(i) = func_T_ref(pwm0_each(i));
    Q_eff_init(i) = func_Q_ref(pwm0_each(i));
end


%% ===================================================================
%  5) Timeseries pros From Workspace blocks
%  ===================================================================
pwm1_in = timeseries(pwm1(idx_w), time_rel, 'Name','pwm1_in');
pwm2_in = timeseries(pwm2(idx_w), time_rel, 'Name','pwm2_in');
pwm3_in = timeseries(pwm3(idx_w), time_rel, 'Name','pwm3_in');
pwm4_in = timeseries(pwm4(idx_w), time_rel, 'Name','pwm4_in');

% Referências (medições — usadas em Scopes pra overlay)
p_out = timeseries(gyrX(idx_w), time_rel, 'Name','p_out');     % rad/s
q_out = timeseries(gyrY(idx_w), time_rel, 'Name','q_out');
r_out = timeseries(gyrZ(idx_w), time_rel, 'Name','r_out');

phi_out   = timeseries(roll_deg(idx_w),  time_rel, 'Name','phi_out');     % graus
theta_out = timeseries(pitch_deg(idx_w), time_rel, 'Name','theta_out');
psi_out   = timeseries(yaw_deg(idx_w),   time_rel, 'Name','psi_out');

accX_out = timeseries(accX(idx_w), time_rel, 'Name','accX_out');          % m/s²
accY_out = timeseries(accY(idx_w), time_rel, 'Name','accY_out');
accZ_out = timeseries(accZ(idx_w), time_rel, 'Name','accZ_out');

% Compatibilidade com nomes antigos do v4 (u_dot_out etc.)
u_dot_out = accX_out;
v_dot_out = accY_out;
w_dot_out = accZ_out;


%% ===================================================================
%  6) Condições iniciais de atitude (escalares — Integrators usam internal IC)
%  ===================================================================
phi0   = deg2rad(roll_deg(find(idx_w, 1, 'first')));
theta0 = deg2rad(pitch_deg(find(idx_w, 1, 'first')));
psi0   = deg2rad(yaw_deg(find(idx_w, 1, 'first')));


%% ===================================================================
%  7) Struct ref auxiliar (pra comparação posterior em script)
%  ===================================================================
ref = struct();
ref.time     = time_rel;
ref.time_abs = time_w;
ref.pqr      = [gyrX(idx_w), gyrY(idx_w), gyrZ(idx_w)];
ref.acc      = [accX(idx_w), accY(idx_w), accZ(idx_w)];
ref.att_deg  = [roll_deg(idx_w), pitch_deg(idx_w), yaw_deg(idx_w)];
ref.pwm      = [pwm1(idx_w), pwm2(idx_w), pwm3(idx_w), pwm4(idx_w)];
ref.P_J      = P_J;


%% ===================================================================
%  Resumo
%  ===================================================================
fprintf('\nsetup_quad_v4 pronto:\n');
fprintf('  janela [%g, %g] s  (sim t=0 ate %.2f s, %d amostras dt=0.1)\n', ...
        t_window(1), t_window(2), t_sim, sum(idx_w));
fprintf('  IC: phi0=%.2f° theta0=%.2f° psi0=%.2f°\n', ...
        rad2deg(phi0), rad2deg(theta0), rad2deg(psi0));
fprintf('  mass=%.4f kg  g=%.4f  tau_motor=%.4f s\n', mass, g_acc, tau_motor);
fprintf('  P_estimated(1:9) = [G1..G8, invJy] = %s\n', mat2str(P_estimated(1:9)', 4));
fprintf('  k_T = [%s]\n', strjoin(arrayfun(@(x) sprintf('%.3f', x), P_estimated(10:13), 'UniformOutput', false), ', '));
fprintf('  k_Q = [%s]\n', strjoin(arrayfun(@(x) sprintf('%.3f', x), P_estimated(14:17), 'UniformOutput', false), ', '));
fprintf('  Dp=%.3f Dq=%.3f Dr=%.3f\n', P_estimated(18), P_estimated(19), P_estimated(20));
fprintf('\nPra simular:\n');
fprintf('  simOut = sim(''quad_model_v4'', ''StopTime'', num2str(t_sim));\n');
