%SETUP_QUAD_V4  Prepara workspace para simular quad_model_v4.slx em uma janela
%               de tempo configuravel do log_data.mat.
%
% Variaveis populadas no base workspace:
%   pwm1_in..pwm4_in            (timeseries Nx1 cada)
%   p_out, q_out, r_out         (timeseries de referencia, rad/s)
%   phi_out, theta_out, psi_out (timeseries de referencia, GRAUS)
%   u_dot_out, v_dot_out, w_dot_out (timeseries de referencia, m/s^2)
%   phi0, theta0, psi0          (timeseries com IC; v4 usa external IC nos
%                                Integrators de pqr2euler1)
%   P_estimated (23x1)          (G-formulation; convertido de P_final via
%                                P_J_to_simulink se houver P_identified.mat)
%   p_bias, q_bias, r_bias      (escalares)
%
% Janela e modo configuraveis abaixo.

% ===================================================================
%  ESCOLHA A JANELA
% ===================================================================
t_window = [147, 157];   %  [t_inicio, t_fim] em segundos
%   Exemplos:
%     [147, 157]  ← janela de validacao
%     [157, 167]  ← primeiro segmento de treino
%     [167, 177]  ← segundo segmento de treino
%     [177, 187]  ← terceiro segmento de treino

% ===================================================================
%  ESCOLHA OS PARAMETROS
% ===================================================================
%   false -> usa P_final (resultado identificado em P_identified.mat)
%   true  -> usa P0 (chute inicial)
use_P0 = false;

% DINAMICA DO MOTOR (lag de 1a ordem entre Polyval e K_T/K_Q no v4):
%   T_motor(s) = 1 / (tau_motor * s + 1)   [ganho estatico = 1]
%   Para desabilitar (motor "instantaneo"): tau_motor = 1e-9
%   Valores tipicos BLDC + ESC: 0.03 ~ 0.10 s
tau_motor = 0.05;   % segundos


%% ===================================================================
%  1) Carregar log e definir P_J (J-formulation)
%  ===================================================================
addpath(fileparts(fileparts(mfilename('fullpath'))));   % raiz pra setup_paths
paths = setup_paths();
load(fullfile(paths.data, 'log_data.mat'));

% Tudo vem de parameters.m (centralizado em 2_model/)
proj_params = parameters();
P0_J = proj_params.P0_J;

% v4.slx lê 'mass' e 'g_acc' do base workspace (blocos Massa e Gravidade)
mass  = proj_params.m;
g_acc = proj_params.g;

P_file = fullfile(paths.outputs, 'P_identified.mat');
if use_P0
    P_J = P0_J;
    fprintf('setup_quad_v4: usando P0 (chute inicial)\n');
elseif exist(P_file, 'file')
    Pdat = load(P_file);
    P_J  = Pdat.P_final;
    fprintf('setup_quad_v4: usando P_final de P_identified.mat\n');
else
    P_J = P0_J;
    fprintf('setup_quad_v4: P_identified.mat nao achado, usando P0 (fallback)\n');
end


%% ===================================================================
%  2) Converter J -> G (formato esperado pelo quad_model_v4)
%      P_estimated(1:8)   = G1..G8
%      P_estimated(9)     = invJy
%      P_estimated(10:13) = k_T1..k_T4  (slots reservados; v4 nao usa)
%      P_estimated(14:17) = k_Q1..k_Q4  (slots reservados; v4 nao usa)
%      P_estimated(18:20) = Dp, Dq, Dr
%      P_estimated(21:23) = Bp, Bq, Br  (slots reservados; v4 le p_bias/q_bias/r_bias)
%  ===================================================================
P_estimated = P_J_to_simulink(P_J);  % wrapper ja existente

% v4 le os biases de variaveis avulsas (nao de P_estimated)
p_bias = P_J(16);
q_bias = P_J(17);
r_bias = P_J(18);


%% ===================================================================
%  3) Reamostrar log na grade comum e cortar para t_window
%  ===================================================================
ATT.TimeS  = double(ATT.TimeUS)  / 1e6;
IMU.TimeS  = double(IMU.TimeUS)  / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;

idx = IMU.I == 0;
gyrX = IMU.GyrX(idx); gyrY = IMU.GyrY(idx); gyrZ = IMU.GyrZ(idx);
accX = IMU.AccX(idx); accY = IMU.AccY(idx); accZ = IMU.AccZ(idx);
time_IMU = IMU.TimeS(idx);

t_start  = max([min(time_IMU), min(ATT.TimeS), min(RCOU.TimeS)]);
t_end    = min([max(time_IMU), max(ATT.TimeS), max(RCOU.TimeS)]);
t_common = (t_start:0.1:t_end)';

gyrX = interp1(time_IMU, gyrX, t_common, 'linear');
gyrY = interp1(time_IMU, gyrY, t_common, 'linear');
gyrZ = interp1(time_IMU, gyrZ, t_common, 'linear');
accX = interp1(time_IMU, accX, t_common, 'linear');
accY = interp1(time_IMU, accY, t_common, 'linear');
accZ = interp1(time_IMU, accZ, t_common, 'linear');
roll_deg  = interp1(ATT.TimeS,  ATT.Roll,  t_common, 'linear');
pitch_deg = interp1(ATT.TimeS,  ATT.Pitch, t_common, 'linear');
yaw_deg   = interp1(ATT.TimeS,  ATT.Yaw,   t_common, 'linear');
pwm1 = interp1(RCOU.TimeS, double(RCOU.C1), t_common, 'linear');
pwm2 = interp1(RCOU.TimeS, double(RCOU.C2), t_common, 'linear');
pwm3 = interp1(RCOU.TimeS, double(RCOU.C3), t_common, 'linear');
pwm4 = interp1(RCOU.TimeS, double(RCOU.C4), t_common, 'linear');

if t_window(1) < t_common(1) || t_window(2) > t_common(end)
    error('Janela [%g, %g] fora do log [%g, %g].', ...
        t_window(1), t_window(2), t_common(1), t_common(end));
end

idx_w   = (t_common >= t_window(1)) & (t_common <= t_window(2));
time_w  = t_common(idx_w);
time_rel = time_w - time_w(1);   % v4 comeca em t=0
t_sim   = time_rel(end);


%% ===================================================================
%  3b) Tabela de bancada PWM->T,Q (consumida via spline cubica)
%      No v4: Lookup_n-D blocks com InterpMethod='Akima spline'
%      No .m: interp1 com 'makima' (numericamente identico ao Akima do Simulink)
%  ===================================================================
pwm_bench = [1000; 1200; 1400; 1600; 1800; 2000];
thrust_g  = [0; 143; 328; 532; 784; 843];                  % gramas
torque_Nm = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];     % Nm
T_bench   = thrust_g * 9.80665 / 1000;                      % converter para N
Q_bench   = torque_Nm;

%% ===================================================================
%  3c) Condicoes iniciais dos State-Space do lag de motor (T_eff/Q_eff)
%      lidas no v4 pelos blocos MotorLag_T*/MotorLag_Q*
%      Usamos o valor de regime estacionario no PWM inicial da janela
%  ===================================================================
% Modelos de motor estaticos (mesma spline do .m e do .slx)
T_ref_at = @(pwm) max(0, interp1(pwm_bench, T_bench, ...
                                  min(max(pwm,1000),2000), 'makima','extrap'));
Q_ref_at = @(pwm) max(0, interp1(pwm_bench, Q_bench, ...
                                  min(max(pwm,1000),2000), 'makima','extrap'));

% IC dos State-Space: o estado eh a saida do Polyval CRUA (antes do K_T,
% porque no v4 o K_T vem DEPOIS da State-Space)
T_eff_init = zeros(4,1);
Q_eff_init = zeros(4,1);
pwm0_each = [pwm1(find(idx_w,1,'first')), pwm2(find(idx_w,1,'first')), ...
             pwm3(find(idx_w,1,'first')), pwm4(find(idx_w,1,'first'))];
for i = 1:4
    T_eff_init(i) = T_ref_at(pwm0_each(i));   % CRU (sem k_T)
    Q_eff_init(i) = Q_ref_at(pwm0_each(i));
end

%% ===================================================================
%  4) PWMs separados (1 timeseries por motor)
%  ===================================================================
pwm1_in = timeseries(pwm1(idx_w), time_rel, 'Name','pwm1_in');
pwm2_in = timeseries(pwm2(idx_w), time_rel, 'Name','pwm2_in');
pwm3_in = timeseries(pwm3(idx_w), time_rel, 'Name','pwm3_in');
pwm4_in = timeseries(pwm4(idx_w), time_rel, 'Name','pwm4_in');


%% ===================================================================
%  5) Referencias (saidas medidas) - usadas pelos Scopes do v4
%  ===================================================================
p_out = timeseries(gyrX(idx_w), time_rel, 'Name','p_out');   % rad/s
q_out = timeseries(gyrY(idx_w), time_rel, 'Name','q_out');
r_out = timeseries(gyrZ(idx_w), time_rel, 'Name','r_out');

phi_out   = timeseries(roll_deg(idx_w),  time_rel, 'Name','phi_out');    % graus
theta_out = timeseries(pitch_deg(idx_w), time_rel, 'Name','theta_out');
psi_out   = timeseries(yaw_deg(idx_w),   time_rel, 'Name','psi_out');

u_dot_out = timeseries(accX(idx_w), time_rel, 'Name','u_dot_out');   % m/s^2
v_dot_out = timeseries(accY(idx_w), time_rel, 'Name','v_dot_out');
w_dot_out = timeseries(accZ(idx_w), time_rel, 'Name','w_dot_out');


%% ===================================================================
%  6) Condicoes iniciais para pqr2euler1 (apos Bug #4: ESCALARES,
%     integradores agora usam internal IC)
%  ===================================================================
phi0   = deg2rad(roll_deg(find(idx_w, 1, 'first')));
theta0 = deg2rad(pitch_deg(find(idx_w, 1, 'first')));
psi0   = deg2rad(yaw_deg(find(idx_w, 1, 'first')));


%% ===================================================================
%  7) Struct ref auxiliar (para comparacao posterior em script)
%  ===================================================================
ref = struct();
ref.time     = time_rel;
ref.time_abs = time_w;
ref.pqr      = [gyrX(idx_w), gyrY(idx_w), gyrZ(idx_w)];
ref.acc      = [accX(idx_w), accY(idx_w), accZ(idx_w)];
ref.att_deg  = [roll_deg(idx_w), pitch_deg(idx_w), yaw_deg(idx_w)];
ref.pwm      = [pwm1(idx_w), pwm2(idx_w), pwm3(idx_w), pwm4(idx_w)];

% Tempo de simulacao do .slx
% (v4 espera StopTime no parametro do modelo; use sim('quad_model_v4', t_sim))


%% ===================================================================
%  Resumo
%  ===================================================================
fprintf('\nsetup_quad_v4 pronto:\n');
fprintf('  janela t_window = [%g, %g] s  (sim t=0 ate t=%.2f)\n', t_window(1), t_window(2), t_sim);
fprintf('  amostras na janela: %d (dt=0.1 s)\n', sum(idx_w));
fprintf('  IC: phi0=%.2f deg theta0=%.2f deg psi0=%.2f deg\n', ...
    rad2deg(phi0), rad2deg(theta0), rad2deg(psi0));
fprintf('  P_estimated(1:9) = [G1..G8, invJy] = %s\n', mat2str(P_estimated(1:9)', 4));
fprintf('  p_bias=%.4f q_bias=%.4f r_bias=%.4f\n', p_bias, q_bias, r_bias);
fprintf('\nPara simular o v4:\n');
fprintf('  simOut = sim(''quad_model_v4'', ''StopTime'', num2str(t_sim));\n');
