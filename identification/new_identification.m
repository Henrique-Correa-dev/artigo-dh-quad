% new_identification.m - VTOL multi-maneuver identification (22 params)
%
% Multi-maneuver: usa múltiplas janelas de treino simultaneamente.
% O custo é a soma dos resíduos de todos os segmentos, conforme
% recomendado por Klein & Morelli e Jategaonkar.
%
% Fase A: EEM rotacional (algébrico, dados concatenados) → P_eem
% Fase B: OEM progressivo (1s → 2s → 5s → full, multi-segmento) → P_final
%
% P(1:4)   = inércias [Jx, Jy, Jz, Jxz] → G1-G8, InvJy computados internamente
% P(5:8)   = k_T1..k_T4
% P(9:12)  = k_Q1..k_Q4
% P(13:15) = Dp, Dq, Dr
% P(16:18) = Bp, Bq, Br
% P(19:20) = CG offsets [dx_cg, dy_cg]
% P(21:24) = translacionais [Xu_m, Yv_m, Zw_m, Bz]

%% 1. Configuração de janelas
% EDITE AQUI: adicione quantas janelas de treino quiser.
% Cada linha é [t_inicio, t_fim] em segundos.
% A validação deve ser uma janela que NÃO aparece no treino.
t_trains = {[157, 167]; ...
            [167, 177]; ...
            [177, 187]}

t_val = [147, 157];

%% 2. Data loading and preprocessing
load("log_data.mat")

ATT.TimeS  = double(ATT.TimeUS) / 1e6;
IMU.TimeS  = double(IMU.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;
GPS.TimeS  = double(GPS.TimeUS) / 1e6;

idx = IMU.I == 0;
gyrX_raw = IMU.GyrX(idx); gyrY_raw = IMU.GyrY(idx); gyrZ_raw = IMU.GyrZ(idx);
accX_raw = IMU.AccX(idx); accY_raw = IMU.AccY(idx); accZ_raw = IMU.AccZ(idx);
time_IMU = IMU.TimeS(idx);

time_ATT  = ATT.TimeS;
time_GPS  = GPS.TimeS;
time_RCOU = RCOU.TimeS;

pwm1_raw = double(RCOU.C1); pwm2_raw = double(RCOU.C2);
pwm3_raw = double(RCOU.C3); pwm4_raw = double(RCOU.C4);

t_start  = max([min(time_IMU), min(time_ATT), min(time_GPS), min(time_RCOU)]);
t_end    = min([max(time_IMU), max(time_ATT), max(time_GPS), max(time_RCOU)]);
t_common = t_start:0.1:t_end;
dt = 0.1;

gyrX_interp = interp1(time_IMU, gyrX_raw, t_common, 'linear');
gyrY_interp = interp1(time_IMU, gyrY_raw, t_common, 'linear');
gyrZ_interp = interp1(time_IMU, gyrZ_raw, t_common, 'linear');
accX_interp = interp1(time_IMU, accX_raw, t_common, 'linear');
accY_interp = interp1(time_IMU, accY_raw, t_common, 'linear');
accZ_interp = interp1(time_IMU, accZ_raw, t_common, 'linear');
roll_interp  = interp1(time_ATT, ATT.Roll, t_common, 'linear');
pitch_interp = interp1(time_ATT, ATT.Pitch, t_common, 'linear');
yaw_interp   = interp1(time_ATT, ATT.Yaw, t_common, 'linear');
pwm1_interp = interp1(time_RCOU, pwm1_raw, t_common, 'linear');
pwm2_interp = interp1(time_RCOU, pwm2_raw, t_common, 'linear');
pwm3_interp = interp1(time_RCOU, pwm3_raw, t_common, 'linear');
pwm4_interp = interp1(time_RCOU, pwm4_raw, t_common, 'linear');

%% 3. Motor reference models
pwm_values_exp = [1000; 1200; 1400; 1600; 1800; 2000];
thrust_grams_exp = [0; 143; 328; 532; 784; 843];
torque_Nm_exp = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];
poly_degree = 3;
func_T_ref = create_thrust_model(pwm_values_exp, thrust_grams_exp, poly_degree);
func_Q_ref = create_torque_model(pwm_values_exp, torque_Nm_exp, poly_degree);

%% 4. Extract training segments
n_seg = length(t_trains);
segs = cell(n_seg, 1);

fprintf('\n  Segmentos de treino: %d\n', n_seg);
for s = 1:n_seg
    idx_s = (t_common >= t_trains{s}(1)) & (t_common <= t_trains{s}(2));
    seg.time = t_common(idx_s)';
    seg.N    = sum(idx_s);
    seg.pwm  = [pwm1_interp(idx_s)', pwm2_interp(idx_s)', ...
                pwm3_interp(idx_s)', pwm4_interp(idx_s)'];
    seg.pqr  = [gyrX_interp(idx_s)', gyrY_interp(idx_s)', gyrZ_interp(idx_s)'];
    seg.acc  = [accX_interp(idx_s)', accY_interp(idx_s)', accZ_interp(idx_s)'];
    seg.att  = [roll_interp(idx_s)', pitch_interp(idx_s)', yaw_interp(idx_s)'];
    seg.att_rad = deg2rad(seg.att);
    seg.T_ref = zeros(seg.N, 4);
    seg.Q_ref = zeros(seg.N, 4);
    for j = 1:4
        seg.T_ref(:,j) = func_T_ref(seg.pwm(:,j));
        seg.Q_ref(:,j) = func_Q_ref(seg.pwm(:,j));
    end
    segs{s} = seg;
    fprintf('    [%d] %d-%ds  (%d pontos)\n', s, t_trains{s}(1), t_trains{s}(2), seg.N);
end

% Validation data
idx_val = (t_common >= t_val(1)) & (t_common <= t_val(2));
time_vl = t_common(idx_val)';
pwm_vl  = [pwm1_interp(idx_val)', pwm2_interp(idx_val)', ...
            pwm3_interp(idx_val)', pwm4_interp(idx_val)'];
pqr_vl  = [gyrX_interp(idx_val)', gyrY_interp(idx_val)', gyrZ_interp(idx_val)'];
acc_vl  = [accX_interp(idx_val)', accY_interp(idx_val)', accZ_interp(idx_val)'];
att_vl  = [roll_interp(idx_val)', pitch_interp(idx_val)', yaw_interp(idx_val)'];

%% 5. Parameters: initial guess and physical bounds (24 params)
% P(1:4)   = inércias [Jx, Jy, Jz, Jxz]
% P(5:8)   = k_T1..k_T4
% P(9:12)  = k_Q1..k_Q4
% P(13:15) = Dp, Dq, Dr
% P(16:18) = Bp, Bq, Br
% P(19:20) = dx_cg, dy_cg  (offset do CG vs CAD, afeta braços de momento)
% P(21:24) = Xu_m, Yv_m, Zw_m, Bz

Jx0  = 63.244/1000;
Jy0  = 250.554/1000;
Jz0  = 116.192/1000;
Jxz0 = 1.571/1000;

P0 = [Jx0; Jy0; Jz0; Jxz0;           ... % Inércias
      0.55; 0.45; 1.0; 0.75;          ... % k_T
      0.55; 0.45; 1.0; 0.75;          ... % k_Q
      10; 5; 0.5;                  ... % Dp, Dq, Dr
      0.7; 1.4; 0.3;                  ... % Bp, Bq, Br
      0.0; 0.0;                        ... % dx_cg, dy_cg (chute: CG no CAD)
      -4.0; -4.0; -0.1; -0.05];            % Xu_m, Yv_m, Zw_m, Bz

lb = [0.032;  0.125;  0.058;  0.0001;   ... % Jx, Jy, Jz, Jxz (±50% do CAD)
      0.05; 0.05; 0.05; 0.05;         ... % k_T
      0.10; 0.10; 0.10; 0.10;         ... % k_Q (mín 10% do ref — fisicamente razoável)
         0;    0;    0;                ... % Dp, Dq, Dr
       -10;  -10;  -10;               ... % Bp, Bq, Br
     -0.08; -0.05;                     ... % dx_cg, dy_cg (±8cm, ±5cm)
       -30;  -30;  -2;  -5];            % Xu_m, Yv_m, Zw_m, Bz

ub = [ 0.095;  0.376;  0.174;  0.006;  ... % Jx, Jy, Jz, Jxz (±50% do CAD)
         5;    5;    5;    5;          ... % k_T
       3.0;  3.0;  3.0;  3.0;         ... % k_Q
        20;    10;    10;                ... % Dp, Dq, Dr
        10;   10;   10;                ... % Bp, Bq, Br
      0.08;  0.05;                     ... % dx_cg, dy_cg
         0;    0;    0;   5];              % Xu_m, Yv_m, Zw_m, Bz

n_params = 24;
n_rot = 20;  % EEM otimiza P(1:20) — rotacional + CG offsets
param_names = {'Jx','Jy','Jz','Jxz', ...
    'k_T1','k_T2','k_T3','k_T4','k_Q1','k_Q2','k_Q3','k_Q4', ...
    'Dp','Dq','Dr','Bp','Bq','Br', ...
    'dx_cg','dy_cg', ...
    'Xu_m','Yv_m','Zw_m','Bz'};

R2_func = @(y_e, y_s) 1 - sum((y_e - y_s).^2) / max(sum((y_e - mean(y_e)).^2), 1e-12);
constants_sim.m = 1.6011;
constants_sim.g = 9.81;
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);  % Mesmo que Simulink

% Criar pasta de imagens se não existir
img_dir = fullfile('C:/Users/Henrique/ARTIGO/identification/images');
if ~exist(img_dir, 'dir'), mkdir(img_dir); end

%% ========================================================================
%  6. VALIDAÇÃO COM CHUTE INICIAL (P0) - primeiro segmento + validação
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  VALIDAÇÃO DO MODELO COM P0 (chute inicial)\n');
fprintf('==========================================================\n');

seg1 = segs{1};

% === DIAGNÓSTICO: comparar valores com Simulink ===
diag_P = P0;
diag_Jx = diag_P(1); diag_Jy = diag_P(2); diag_Jz = diag_P(3); diag_Jxz = diag_P(4);
diag_gamma = diag_Jx*diag_Jz - diag_Jxz^2;
fprintf('\n  --- Constantes do modelo (compare com Simulink) ---\n');
fprintf('  Jx=%.6f  Jy=%.6f  Jz=%.6f  Jxz=%.6f\n', diag_Jx, diag_Jy, diag_Jz, diag_Jxz);
fprintf('  gamma0 = %.6f\n', diag_gamma);
fprintf('  G1=%.6f  G2=%.6f  G3=%.6f  G4=%.6f\n', ...
    diag_Jxz*(diag_Jx-diag_Jy+diag_Jz)/diag_gamma, ...
    (diag_Jz*(diag_Jz-diag_Jy)+diag_Jxz^2)/diag_gamma, ...
    diag_Jz/diag_gamma, diag_Jxz/diag_gamma);
fprintf('  G5=%.6f  G6=%.6f  G7=%.6f  G8=%.6f\n', ...
    (diag_Jz-diag_Jx)/diag_Jy, diag_Jxz/diag_Jy, ...
    (diag_Jx*(diag_Jx-diag_Jy)+diag_Jxz^2)/diag_gamma, diag_Jx/diag_gamma);
fprintf('  invJy = %.6f\n', 1/diag_Jy);
fprintf('  k_T = [%.4f, %.4f, %.4f, %.4f]\n', diag_P(5:8));
fprintf('  k_Q = [%.4f, %.4f, %.4f, %.4f]\n', diag_P(9:12));
fprintf('  Dp=%.4f  Dq=%.4f  Dr=%.4f\n', diag_P(13:15));
fprintf('  Bp=%.4f  Bq=%.4f  Br=%.4f\n', diag_P(16:18));
fprintf('  dx_cg=%.4f  dy_cg=%.4f\n', diag_P(19:20));
fprintf('  Braços: Lx_r=%.6f  Lx_l=%.6f  Ly_f=%.6f  Ly_r=%.6f\n', ...
    0.232-diag_P(20), 0.232+diag_P(20), 0.311185-diag_P(19), 0.342865+diag_P(19));

% Derivadas no instante t=0 (para comparação pontual com Simulink)
t0 = seg1.time(1);
y0_rot = seg1.pqr(1,:)';
dy0 = vtol_dynamics(t0, y0_rot, P0, seg1.time, seg1.pwm, func_T_ref, func_Q_ref);
pwm0 = seg1.pwm(1,:);
T0 = diag_P(5:8)' .* [func_T_ref(pwm0(1)), func_T_ref(pwm0(2)), func_T_ref(pwm0(3)), func_T_ref(pwm0(4))];
Q0 = diag_P(9:12)' .* [func_Q_ref(pwm0(1)), func_Q_ref(pwm0(2)), func_Q_ref(pwm0(3)), func_Q_ref(pwm0(4))];
Lx_r = 0.232-diag_P(20); Lx_l = 0.232+diag_P(20);
Ly_f = 0.311185-diag_P(19); Ly_r = 0.342865+diag_P(19);
Mx0 = -(Lx_r*T0(1) - Lx_l*T0(2) - Lx_l*T0(3) + Lx_r*T0(4));
My0 = Ly_f*T0(1) - Ly_r*T0(2) + Ly_f*T0(3) - Ly_r*T0(4);
Mz0 = Q0(1) + Q0(2) - Q0(3) - Q0(4);
fprintf('\n  --- Estado inicial (t=%.1fs) ---\n', t0);
fprintf('  PWM = [%.0f, %.0f, %.0f, %.0f] us\n', pwm0);
fprintf('  T_ref(PWM) = [%.4f, %.4f, %.4f, %.4f] N\n', ...
    func_T_ref(pwm0(1)), func_T_ref(pwm0(2)), func_T_ref(pwm0(3)), func_T_ref(pwm0(4)));
fprintf('  T_mr (k_T*T_ref) = [%.4f, %.4f, %.4f, %.4f] N\n', T0);
fprintf('  Q_mr (k_Q*Q_ref) = [%.6f, %.6f, %.6f, %.6f] N·m\n', Q0);
fprintf('  Mx=%.6f  My=%.6f  Mz=%.6f  N·m\n', Mx0, My0, Mz0);
fprintf('  p0=%.4f  q0=%.4f  r0=%.4f  rad/s\n', y0_rot);
fprintf('  p_dot=%.6f  q_dot=%.6f  r_dot=%.6f  rad/s²\n', dy0);
fprintf('  -----------------------------------------------\n\n');

[res_P0] = simulate_full(P0, seg1.time, seg1.pwm, seg1.pqr, seg1.acc, seg1.att, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, ode_opts);

print_R2('P0', res_P0, seg1.pqr, pqr_vl, seg1.acc, acc_vl, t_trains{1}, t_val, R2_func);
plot_all_results('CHUTE INICIAL (P0)', res_P0, seg1.time, seg1.pqr, seg1.acc, ...
    time_vl, pqr_vl, acc_vl, R2_func, t_trains{1}, t_val, att_vl);

%% ========================================================================
%  7. FASE A: EEM rotacional (dados concatenados de todos os segmentos)
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  FASE A: EEM rotacional (%d segmentos, %d params)\n', n_seg, n_rot);
fprintf('==========================================================\n');

smooth_win = 5;
p_all = []; q_all = []; r_all = [];
pd_all = []; qd_all = []; rd_all = [];
Tr_all = []; Qr_all = [];

for s = 1:n_seg
    sg = segs{s};
    p_s = sg.pqr(:,1); q_s = sg.pqr(:,2); r_s = sg.pqr(:,3);
    p_all = [p_all; p_s]; %#ok<AGROW>
    q_all = [q_all; q_s]; %#ok<AGROW>
    r_all = [r_all; r_s]; %#ok<AGROW>
    pd_all = [pd_all; gradient(movmean(p_s, smooth_win), dt)]; %#ok<AGROW>
    qd_all = [qd_all; gradient(movmean(q_s, smooth_win), dt)]; %#ok<AGROW>
    rd_all = [rd_all; gradient(movmean(r_s, smooth_win), dt)]; %#ok<AGROW>
    Tr_all = [Tr_all; sg.T_ref]; %#ok<AGROW>
    Qr_all = [Qr_all; sg.Q_ref]; %#ok<AGROW>
end

var_pd = var(pd_all); if var_pd < 1e-12, var_pd = 1; end
var_qd = var(qd_all); if var_qd < 1e-12, var_qd = 1; end
var_rd = var(rd_all); if var_rd < 1e-12, var_rd = 1; end
weights_eem = [1/var_pd; 1/var_qd; 1/var_rd];

cost_eem = @(Prot) eem_cost_function(ones(n_rot,1), Prot, weights_eem, ...
    p_all, q_all, r_all, pd_all, qd_all, rd_all, Tr_all, Qr_all);

opts_eem = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter', ...
    'MaxIterations', 1500, ...
    'MaxFunctionEvaluations', 60000, ...
    'StepTolerance', 1e-14, ...
    'FunctionTolerance', 1e-14);

[P_eem_rot, rn_eem, ~, ef_eem] = lsqnonlin(cost_eem, P0(1:n_rot), lb(1:n_rot), ub(1:n_rot), opts_eem);
P_eem = [P_eem_rot; P0(n_rot+1:n_params)];

fprintf('\n  EEM Resnorm: %.4f  |  Exit flag: %d\n', rn_eem, ef_eem);
for i = 1:n_rot
    fprintf('    %-6s: %12.6f  (chute: %12.6f)\n', param_names{i}, P_eem(i), P0(i));
end

%% ========================================================================
%  8. FASE B: OEM Progressivo Multi-Segmento (1s → 2s → 5s → full)
%  ========================================================================

% Pesos OEM: computados sobre dados concatenados de todos os segmentos
pqr_cat = cell2mat(cellfun(@(s) s.pqr, segs, 'UniformOutput', false));
acc_cat = cell2mat(cellfun(@(s) s.acc, segs, 'UniformOutput', false));

var_p = var(pqr_cat(:,1)); if var_p < 1e-12, var_p = 1; end
var_q = var(pqr_cat(:,2)); if var_q < 1e-12, var_q = 1; end
var_r = var(pqr_cat(:,3)); if var_r < 1e-12, var_r = 1; end
weights_pqr = [1/var_p; 1/var_q; 1/var_r];

var_ax = var(acc_cat(:,1)); if var_ax < 1e-12, var_ax = 1; end
var_ay = var(acc_cat(:,2)); if var_ay < 1e-12, var_ay = 1; end
var_az = var(acc_cat(:,3)); if var_az < 1e-12, var_az = 1; end
weights_acc = [1/var_ax; 1/var_ay; 1/var_az];

win_durations = [1.0, 2.0, 3.0];
P_current = P_eem;
best_R2_mean = -Inf;
P_best = P_eem;
best_stage_name = 'EEM';

for stage = 1:length(win_durations)
    win_sec = win_durations(stage);
    if isinf(win_sec)
        stage_name = 'FULL';
    else
        stage_name = sprintf('%.0fs', win_sec);
    end

    fprintf('\n==========================================================\n');
    fprintf('  FASE B.%d: OEM [%s] (%d segmentos)\n', stage, stage_name, n_seg);
    fprintf('==========================================================\n');

    cost_oem = @(P) oem_multi_seg_cost(P, segs, win_sec, ...
        constants_sim.m, constants_sim.g, dt, weights_pqr, weights_acc);

    max_iter = 500;
    if isinf(win_sec), max_iter = 1000; end

    opts_oem = optimoptions('lsqnonlin', ...
        'Algorithm', 'trust-region-reflective', ...
        'Display', 'iter', ...
        'MaxIterations', max_iter, ...
        'MaxFunctionEvaluations', 80000, ...
        'StepTolerance', 1e-14, ...
        'FunctionTolerance', 1e-14);

    [P_current, rn, ~, ef] = lsqnonlin(cost_oem, P_current, lb, ub, opts_oem);
    fprintf('  [%s] Resnorm: %.4f | Exit: %d\n', stage_name, rn, ef);

    % Avaliar validação para escolher o melhor estágio
    sg_eval = segs{1};
    res_eval = simulate_full(P_current, sg_eval.time, sg_eval.pwm, sg_eval.pqr, sg_eval.acc, sg_eval.att, ...
        time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
        func_T_ref, func_Q_ref, constants_sim, ode_opts);

    if res_eval.pqr_vl_ok
        R2_p = R2_func(pqr_vl(:,1), res_eval.p_s_vl);
        R2_q = R2_func(pqr_vl(:,2), res_eval.q_s_vl);
        R2_r = R2_func(pqr_vl(:,3), res_eval.r_s_vl);
        R2_mean_stage = mean([R2_p, R2_q, R2_r]);
        fprintf('  [%s] Val R² p=%.4f | q=%.4f | r=%.4f | média=%.4f\n', ...
            stage_name, R2_p, R2_q, R2_r, R2_mean_stage);

        if R2_mean_stage > best_R2_mean
            best_R2_mean = R2_mean_stage;
            P_best = P_current;
            best_stage_name = stage_name;
        end
    else
        fprintf('  [%s] Val -> DIVERGIU\n', stage_name);
    end
end

fprintf('\n  >>> Melhor estágio: [%s] com R² médio = %.4f\n', best_stage_name, best_R2_mean);

P_final = P_best;

fprintf('\n  --- Parâmetros Finais ---\n');
for i = 1:n_params
    fprintf('    %-6s: %12.6f  (EEM: %12.6f)  (P0: %12.6f)\n', ...
        param_names{i}, P_final(i), P_eem(i), P0(i));
end

% Mostrar G's resultantes das inércias identificadas
Jx_f = P_final(1); Jy_f = P_final(2); Jz_f = P_final(3); Jxz_f = P_final(4);
gam0_f = Jx_f*Jz_f - Jxz_f^2;
fprintf('\n  --- Constantes G (derivadas das inércias) ---\n');
fprintf('    G1=%8.4f  G2=%8.4f  G3=%8.4f  G4=%8.4f\n', ...
    Jxz_f*(Jx_f-Jy_f+Jz_f)/gam0_f, (Jz_f*(Jz_f-Jy_f)+Jxz_f^2)/gam0_f, ...
    Jz_f/gam0_f, Jxz_f/gam0_f);
fprintf('    G5=%8.4f  G6=%8.4f  G7=%8.4f  G8=%8.4f  InvJy=%8.4f\n', ...
    (Jz_f-Jx_f)/Jy_f, Jxz_f/Jy_f, (Jx_f*(Jx_f-Jy_f)+Jxz_f^2)/gam0_f, ...
    Jx_f/gam0_f, 1/Jy_f);

%% ========================================================================
%  9. VALIDAÇÃO COM P_final
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  VALIDAÇÃO COM P_final (Rot + Trans separados)\n');
fprintf('==========================================================\n');

sg = segs{1};
[res_final] = simulate_full(P_final, sg.time, sg.pwm, sg.pqr, sg.acc, sg.att, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, ode_opts);

print_R2('P_final', res_final, sg.pqr, pqr_vl, sg.acc, acc_vl, ...
    t_trains{1}, t_val, R2_func);

plot_all_results('P_final', ...
    res_final, sg.time, sg.pqr, sg.acc, time_vl, pqr_vl, acc_vl, R2_func, t_trains{1}, t_val, att_vl);

% Diagnóstico: Torques Mx, My, Mz na janela de validação
plot_torques(P_final, time_vl, pwm_vl, func_T_ref, func_Q_ref, pqr_vl, t_val);

% Diagnóstico: Forças Fx, Fy, Fz na janela de validação
plot_forces(P_final, time_vl, pwm_vl, att_vl, acc_vl, func_T_ref, constants_sim, t_val, res_final);

%% 10. RESUMO FINAL: P0 vs P_final (apenas validação)
fprintf('\n==========================================================\n');
fprintf('  RESUMO FINAL — Validação (%d-%ds)\n', t_val(1), t_val(2));
fprintf('==========================================================\n');
fprintf('  [P0] Chute Inicial:\n');
if res_P0.pqr_vl_ok
    fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
        R2_func(pqr_vl(:,1), res_P0.p_s_vl), R2_func(pqr_vl(:,2), res_P0.q_s_vl), R2_func(pqr_vl(:,3), res_P0.r_s_vl));
else, fprintf('    p,q,r -> DIVERGIU\n');
end
if res_P0.full_vl_ok
    fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
        R2_func(acc_vl(:,1), res_P0.accX_s_vl), R2_func(acc_vl(:,2), res_P0.accY_s_vl), R2_func(acc_vl(:,3), res_P0.accZ_s_vl));
end
fprintf('  [P_final] Após Otimização:\n');
if res_final.pqr_vl_ok
    fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
        R2_func(pqr_vl(:,1), res_final.p_s_vl), R2_func(pqr_vl(:,2), res_final.q_s_vl), R2_func(pqr_vl(:,3), res_final.r_s_vl));
else, fprintf('    p,q,r -> DIVERGIU\n');
end
if res_final.full_vl_ok
    fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
        R2_func(acc_vl(:,1), res_final.accX_s_vl), R2_func(acc_vl(:,2), res_final.accY_s_vl), R2_func(acc_vl(:,3), res_final.accZ_s_vl));
end
fprintf('==========================================================\n');

disp('Script finalizado.');

%% ========================================================================
%  FUNÇÕES LOCAIS
%  ========================================================================

function e = oem_multi_seg_cost(P, segs, win_sec, m, g, dt, weights_pqr, weights_acc, cost_mode)
% Custo OEM multi-segmento: concatena resíduos de cada segmento.
% cost_mode: 'full' (padrão), 'rotational', 'translational'
    if nargin < 9, cost_mode = 'full'; end
    e_all = [];
    for s = 1:length(segs)
        sg = segs{s};
        N = sg.N;

        if isinf(win_sec)
            win_len = N;
        else
            win_len = round(win_sec / dt);
        end

        win_starts = 1:win_len:N;
        if length(win_starts) > 1 && win_starts(end) == N
            win_starts(end) = [];
        end
        win_ends = [win_starts(2:end)-1, N];

        e_seg = oem_ms_cost_func(P, sg.pqr, sg.acc, sg.att_rad, ...
            sg.T_ref, sg.Q_ref, m, g, dt, N, ...
            win_starts, win_ends, weights_pqr, weights_acc, cost_mode);

        e_all = [e_all; e_seg]; %#ok<AGROW>
    end
    e = e_all;
end

function e = oem_ms_cost_func(P, pqr, acc, att_rad, T_ref, Q_ref, m, g, dt, N, ...
    win_starts, win_ends, weights_pqr, weights_acc, cost_mode)

    if nargin < 15, cost_mode = 'full'; end

    % Inércias → constantes G (corpo rígido, consistência garantida)
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
    gamma0 = Jx*Jz - Jxz^2;
    G1 = Jxz*(Jx - Jy + Jz) / gamma0;
    G2 = (Jz*(Jz - Jy) + Jxz^2) / gamma0;
    G3 = Jz / gamma0;
    G4 = Jxz / gamma0;
    G5 = (Jz - Jx) / Jy;
    G6 = Jxz / Jy;
    G7 = (Jx*(Jx - Jy) + Jxz^2) / gamma0;
    G8 = Jx / gamma0;
    invJy = 1 / Jy;

    k_T = P(5:8); k_Q = P(9:12);
    Dp = P(13); Dq = P(14); Dr = P(15);
    Bp = P(16); Bq = P(17); Br = P(18);
    dx_cg = P(19); dy_cg = P(20);
    Xu_m = P(21); Yv_m = P(22); Zw_m = P(23); Bz = P(24);

    % Braços efetivos com offset do CG
    Lx_r = 0.232 - dy_cg;   % direita (motores 1,4)
    Lx_l = 0.232 + dy_cg;   % esquerda (motores 2,3)
    Ly_f = 0.311185 - dx_cg; % frente (motores 1,3)
    Ly_r = 0.342865 + dx_cg; % traseira (motores 2,4)

    n_sub = 5;
    dt_sub = dt / n_sub;
    h2 = dt_sub / 2;
    MAX_VAL = 50;

    p_sim = zeros(N, 1);
    q_sim = zeros(N, 1);
    r_sim = zeros(N, 1);

    for w = 1:length(win_starts)
        i_s = win_starts(w);
        i_e = win_ends(w);
        ps = pqr(i_s,1); qs = pqr(i_s,2); rs = pqr(i_s,3);
        p_sim(i_s) = ps; q_sim(i_s) = qs; r_sim(i_s) = rs;

        for k = i_s:i_e-1
            Tmr = k_T(:)' .* T_ref(k,:);
            Qmr = k_Q(:)' .* Q_ref(k,:);
            Mx = -(Lx_r*Tmr(1) - Lx_l*Tmr(2) - Lx_l*Tmr(3) + Lx_r*Tmr(4));
            My = Ly_f*Tmr(1) - Ly_r*Tmr(2) + Ly_f*Tmr(3) - Ly_r*Tmr(4);
            Mz = Qmr(1) + Qmr(2) - Qmr(3) - Qmr(4);

            for si = 1:n_sub
                % RK4 sub-stepping (Mx, My, Mz constantes no sub-step)
                % k1
                pd1 = G1*ps*qs - G2*qs*rs + G3*Mx + G4*Mz - Dp*ps + Bp;
                qd1 = G5*ps*rs - G6*(ps^2 - rs^2) + invJy*My - Dq*qs + Bq;
                rd1 = G7*ps*qs - G1*qs*rs + G4*Mx + G8*Mz - Dr*rs + Br;
                % k2
                p2 = ps+h2*pd1; q2 = qs+h2*qd1; r2 = rs+h2*rd1;
                pd2 = G1*p2*q2 - G2*q2*r2 + G3*Mx + G4*Mz - Dp*p2 + Bp;
                qd2 = G5*p2*r2 - G6*(p2^2 - r2^2) + invJy*My - Dq*q2 + Bq;
                rd2 = G7*p2*q2 - G1*q2*r2 + G4*Mx + G8*Mz - Dr*r2 + Br;
                % k3
                p3 = ps+h2*pd2; q3 = qs+h2*qd2; r3 = rs+h2*rd2;
                pd3 = G1*p3*q3 - G2*q3*r3 + G3*Mx + G4*Mz - Dp*p3 + Bp;
                qd3 = G5*p3*r3 - G6*(p3^2 - r3^2) + invJy*My - Dq*q3 + Bq;
                rd3 = G7*p3*q3 - G1*q3*r3 + G4*Mx + G8*Mz - Dr*r3 + Br;
                % k4
                p4 = ps+dt_sub*pd3; q4 = qs+dt_sub*qd3; r4 = rs+dt_sub*rd3;
                pd4 = G1*p4*q4 - G2*q4*r4 + G3*Mx + G4*Mz - Dp*p4 + Bp;
                qd4 = G5*p4*r4 - G6*(p4^2 - r4^2) + invJy*My - Dq*q4 + Bq;
                rd4 = G7*p4*q4 - G1*q4*r4 + G4*Mx + G8*Mz - Dr*r4 + Br;
                % update
                ps = ps + dt_sub/6*(pd1 + 2*pd2 + 2*pd3 + pd4);
                qs = qs + dt_sub/6*(qd1 + 2*qd2 + 2*qd3 + qd4);
                rs = rs + dt_sub/6*(rd1 + 2*rd2 + 2*rd3 + rd4);
            end

            if ~isfinite(ps), ps = MAX_VAL; end
            if ~isfinite(qs), qs = MAX_VAL; end
            if ~isfinite(rs), rs = MAX_VAL; end
            ps = max(min(ps, MAX_VAL), -MAX_VAL);
            qs = max(min(qs, MAX_VAL), -MAX_VAL);
            rs = max(min(rs, MAX_VAL), -MAX_VAL);

            p_sim(k+1) = ps; q_sim(k+1) = qs; r_sim(k+1) = rs;
        end
    end

    sw_r = sqrt(weights_pqr(:));
    e_rot = [sw_r(1)*(pqr(:,1) - p_sim); ...
             sw_r(2)*(pqr(:,2) - q_sim); ...
             sw_r(3)*(pqr(:,3) - r_sim)];

    T_total = k_T(1)*T_ref(:,1) + k_T(2)*T_ref(:,2) + ...
              k_T(3)*T_ref(:,3) + k_T(4)*T_ref(:,4);

    phi_r = att_rad(:,1); theta_r = att_rad(:,2);
    gx = -g * sin(theta_r);
    gy = g * cos(theta_r) .* sin(phi_r);
    gz = g * cos(theta_r) .* cos(phi_r);

    u_int = zeros(N,1); v_int = zeros(N,1); w_int = zeros(N,1);
    for k = 1:N-1
        pk = pqr(k,1); qk = pqr(k,2); rk = pqr(k,3);
        ud = rk*v_int(k) - qk*w_int(k) + gx(k) + Xu_m*u_int(k);
        vd = pk*w_int(k) - rk*u_int(k) + gy(k) + Yv_m*v_int(k);
        wd = qk*u_int(k) - pk*v_int(k) - T_total(k)/m + gz(k) + Zw_m*w_int(k) + Bz;
        u_int(k+1) = u_int(k) + dt*ud;
        v_int(k+1) = v_int(k) + dt*vd;
        w_int(k+1) = w_int(k) + dt*wd;
    end

    % Força específica = aceleração inercial - gravidade_body
    % IMU mede f_esp, não aceleração inercial
    accX_m = pqr(:,3).*v_int - pqr(:,2).*w_int + Xu_m*u_int;
    accY_m = pqr(:,1).*w_int - pqr(:,3).*u_int + Yv_m*v_int;
    accZ_m = pqr(:,2).*u_int - pqr(:,1).*v_int - T_total/m + Zw_m*w_int + Bz;

    sw_a = sqrt(weights_acc(:));
    e_acc = [sw_a(1)*(acc(:,1) - accX_m); ...
             sw_a(2)*(acc(:,2) - accY_m); ...
             sw_a(3)*(acc(:,3) - accZ_m)];

    if strcmp(cost_mode, 'rotational')
        e = e_rot;
    elseif strcmp(cost_mode, 'translational')
        e = e_acc;
    else
        e = [e_rot; e_acc];
    end
end

function res = simulate_full(P, time_tr, pwm_tr, pqr_tr, acc_tr, att_tr, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, ode_opts)
% Simulação completa model-driven (replica estrutura do Simulink):
%   Etapa 1: Rotacional  (3 estados) → p, q, r
%   Etapa 2: Cinemática  (3 estados) → phi, theta, psi  (com p,q,r simulados)
%   Etapa 3: Translacional (3 estados) → u, v, w        (com p,q,r e att simulados)
% Cada etapa usa ode45 independente — evita contaminação numérica entre subsistemas.
% Totalmente model-driven: nenhum dado medido alimenta a simulação.

    N_tr = length(time_tr);
    N_vl = length(time_vl);
    g = constants_sim.g;
    m = constants_sim.m;

    res = simulate_one_window(P, time_tr, pwm_tr, pqr_tr, att_tr, N_tr, ...
        func_T_ref, func_Q_ref, constants_sim, ode_opts, g, 'tr');
    res2 = simulate_one_window(P, time_vl, pwm_vl, pqr_vl, att_vl, N_vl, ...
        func_T_ref, func_Q_ref, constants_sim, ode_opts, g, 'vl');

    % Mesclar resultados validação no struct
    fnames = fieldnames(res2);
    for i = 1:length(fnames)
        res.(fnames{i}) = res2.(fnames{i});
    end
end

function res = simulate_one_window(P, time, pwm, pqr, att, N, ...
    func_T_ref, func_Q_ref, constants_sim, ode_opts, g, suffix)
% Integração sequencial de uma janela (treino ou validação)

    nan3 = NaN(N,1);

    % === Etapa 1: Rotacional (3 estados) ===
    ode_rot = @(t,y) vtol_dynamics(t, y, P, time, pwm, func_T_ref, func_Q_ref);
    try
        [t_s, y_s] = ode45(ode_rot, time, pqr(1,:)', ode_opts);
        pqr_sim = interp1(t_s, y_s, time, 'linear', 'extrap');
        res.(['p_s_' suffix]) = pqr_sim(:,1);
        res.(['q_s_' suffix]) = pqr_sim(:,2);
        res.(['r_s_' suffix]) = pqr_sim(:,3);
        res.(['pqr_' suffix '_ok']) = true;
    catch
        res.(['p_s_' suffix]) = nan3; res.(['q_s_' suffix]) = nan3; res.(['r_s_' suffix]) = nan3;
        res.(['phi_s_' suffix]) = nan3; res.(['theta_s_' suffix]) = nan3; res.(['psi_s_' suffix]) = nan3;
        res.(['accX_s_' suffix]) = nan3; res.(['accY_s_' suffix]) = nan3; res.(['accZ_s_' suffix]) = nan3;
        res.(['pqr_' suffix '_ok']) = false;
        res.(['full_' suffix '_ok']) = false;
        return;
    end

    % === Etapa 2: Cinemática Euler (3 estados, usa p,q,r SIMULADOS) ===
    att_rad0 = deg2rad(att(1,:))';
    ode_kin = @(t,y) kin_ode(t, y, time, pqr_sim);
    try
        [t_s, y_s] = ode45(ode_kin, time, att_rad0, ode_opts);
        att_sim_rad = interp1(t_s, y_s, time, 'linear', 'extrap');
        res.(['phi_s_' suffix])   = rad2deg(att_sim_rad(:,1));
        res.(['theta_s_' suffix]) = rad2deg(att_sim_rad(:,2));
        res.(['psi_s_' suffix])   = rad2deg(att_sim_rad(:,3));
    catch
        att_sim_rad = repmat(att_rad0', N, 1);  % fallback: atitude constante
        res.(['phi_s_' suffix]) = nan3; res.(['theta_s_' suffix]) = nan3; res.(['psi_s_' suffix]) = nan3;
    end

    % === Etapa 3: Translacional (3 estados, usa p,q,r e att SIMULADOS) ===
    ode_trans = @(t,y) trans_ode(t, y, P, time, pwm, pqr_sim, att_sim_rad, func_T_ref, constants_sim);
    try
        [t_s, y_s] = ode45(ode_trans, time, [0;0;0], ode_opts);
        uvw_sim = interp1(t_s, y_s, time, 'linear', 'extrap');

        % Força específica: acc_modelo = u_dot - g_body
        accX = zeros(N,1); accY = zeros(N,1); accZ = zeros(N,1);
        for k = 1:N
            dy = trans_ode(time(k), uvw_sim(k,:)', P, time, pwm, pqr_sim, att_sim_rad, func_T_ref, constants_sim);
            phi_k   = att_sim_rad(k,1);
            theta_k = att_sim_rad(k,2);
            gx_k = -g * sin(theta_k);
            gy_k =  g * cos(theta_k) * sin(phi_k);
            gz_k =  g * cos(theta_k) * cos(phi_k);
            accX(k) = dy(1) - gx_k;
            accY(k) = dy(2) - gy_k;
            accZ(k) = dy(3) - gz_k;
        end
        res.(['accX_s_' suffix]) = accX;
        res.(['accY_s_' suffix]) = accY;
        res.(['accZ_s_' suffix]) = accZ;
        res.(['full_' suffix '_ok']) = true;
    catch
        res.(['accX_s_' suffix]) = nan3; res.(['accY_s_' suffix]) = nan3; res.(['accZ_s_' suffix]) = nan3;
        res.(['full_' suffix '_ok']) = false;
    end
end

% =========================================================================
%  ODE auxiliar: Cinemática de Euler  y = [phi; theta; psi]
%  Usa p,q,r simulados interpolados
% =========================================================================
function dydt = kin_ode(t, y, time_data, pqr_data)
    phi = y(1); theta = y(2);
    p = interp1(time_data, pqr_data(:,1), t, 'linear', 'extrap');
    q = interp1(time_data, pqr_data(:,2), t, 'linear', 'extrap');
    r = interp1(time_data, pqr_data(:,3), t, 'linear', 'extrap');

    ct = cos(theta);
    if abs(ct) < 1e-7, ct = 1e-7 * sign(ct); end
    sp = sin(phi); cp = cos(phi);
    tt = sin(theta) / ct;

    dydt = [p + (q*sp + r*cp)*tt;
            q*cp - r*sp;
            (q*sp + r*cp)/ct];
end

% =========================================================================
%  ODE auxiliar: Translacional  y = [u; v; w]
%  Usa p,q,r e atitude SIMULADOS interpolados
% =========================================================================
function dydt = trans_ode(t, y, P, time_data, pwm_data, pqr_data, att_data, func_T_ref, constants)
    u = y(1); v = y(2); w = y(3);

    p = interp1(time_data, pqr_data(:,1), t, 'linear', 'extrap');
    q = interp1(time_data, pqr_data(:,2), t, 'linear', 'extrap');
    r = interp1(time_data, pqr_data(:,3), t, 'linear', 'extrap');
    phi   = interp1(time_data, att_data(:,1), t, 'linear', 'extrap');
    theta = interp1(time_data, att_data(:,2), t, 'linear', 'extrap');
    psi   = interp1(time_data, att_data(:,3), t, 'linear', 'extrap');

    k_T = P(5:8);
    cpwm = zeros(1,4);
    for i = 1:4
        cpwm(i) = interp1(time_data, pwm_data(:,i), t, 'linear', 'extrap');
    end
    T_total = k_T(1)*func_T_ref(cpwm(1)) + k_T(2)*func_T_ref(cpwm(2)) + ...
              k_T(3)*func_T_ref(cpwm(3)) + k_T(4)*func_T_ref(cpwm(4));

    m_body = constants.m;  g_acc = constants.g;
    if length(P) >= 24
        Xu_m = P(21); Yv_m = P(22); Zw_m = P(23); Bz = P(24);
    else
        Xu_m = -4.0; Yv_m = -4.0; Zw_m = -0.1; Bz = -0.5;
    end

    sp = sin(phi); cp = cos(phi); ct = cos(theta);
    R_nb = [ cos(theta)*cos(psi), sp*sin(theta)*cos(psi)-cp*sin(psi), cp*sin(theta)*cos(psi)+sp*sin(psi);
             cos(theta)*sin(psi), sp*sin(theta)*sin(psi)+cp*cos(psi), cp*sin(theta)*sin(psi)-sp*cos(psi);
            -sin(theta),          sp*ct,                               cp*ct];
    G_body = R_nb' * [0; 0; m_body*g_acc];

    Fx = G_body(1);
    Fy = G_body(2);
    Fz = -T_total + G_body(3);

    dydt = [r*v - q*w + Fx/m_body + Xu_m*u;
            p*w - r*u + Fy/m_body + Yv_m*v;
            q*u - p*v + Fz/m_body + Zw_m*w + Bz];
end

function print_R2(label, res, ~, pqr_vl, ~, acc_vl, ~, t_val, R2_func)
    fprintf('  [%s] Validação (%d-%ds):\n', label, t_val(1), t_val(2));
    if res.pqr_vl_ok
        fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
            R2_func(pqr_vl(:,1), res.p_s_vl), R2_func(pqr_vl(:,2), res.q_s_vl), R2_func(pqr_vl(:,3), res.r_s_vl));
    else, fprintf('    p,q,r -> DIVERGIU\n');
    end
    if res.full_vl_ok
        fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
            R2_func(acc_vl(:,1), res.accX_s_vl), R2_func(acc_vl(:,2), res.accY_s_vl), R2_func(acc_vl(:,3), res.accZ_s_vl));
    end
end

function plot_all_results(titulo, res, ~, ~, ~, ...
    time_vl, pqr_vl, acc_vl, R2_func, ~, t_val_range, att_vl)

    lbl_vl = sprintf('VALIDAÇÃO [%d-%ds]', t_val_range(1), t_val_range(2));
    safe_titulo = strrep(titulo, ' ', '_');
    safe_titulo = regexprep(safe_titulo, '[^a-zA-Z0-9_]', '');

    fig1 = figure('Name', ['pqr - ' titulo], 'Position', [80 50 900 700], 'Visible', 'off');
    subplot(3,1,1);
    plot(time_vl, pqr_vl(:,1), 'b-', time_vl, res.p_s_vl, 'r--', 'LineWidth', 1.3);
    legend('Exp','Sim'); ylabel('p (rad/s)');
    title(sprintf('%s — p (R^2=%.3f)', lbl_vl, R2_func(pqr_vl(:,1), res.p_s_vl))); grid on;
    subplot(3,1,2);
    plot(time_vl, pqr_vl(:,2), 'b-', time_vl, res.q_s_vl, 'r--', 'LineWidth', 1.3);
    legend('Exp','Sim'); ylabel('q (rad/s)');
    title(sprintf('%s — q (R^2=%.3f)', lbl_vl, R2_func(pqr_vl(:,2), res.q_s_vl))); grid on;
    subplot(3,1,3);
    plot(time_vl, pqr_vl(:,3), 'b-', time_vl, res.r_s_vl, 'r--', 'LineWidth', 1.3);
    legend('Exp','Sim'); xlabel('Tempo (s)'); ylabel('r (rad/s)');
    title(sprintf('%s — r (R^2=%.3f)', lbl_vl, R2_func(pqr_vl(:,3), res.r_s_vl))); grid on;
    sgtitle(['Vel. Angulares (Val) — ' titulo]);
    saveas(fig1, fullfile('C:/Users/Henrique/ARTIGO/identification/images', ['pqr_val_' safe_titulo '.png']));

    fig2 = figure('Name', ['Acc - ' titulo], 'Position', [120 90 900 700], 'Visible', 'off');
    subplot(3,1,1);
    plot(time_vl, acc_vl(:,1), 'b-', time_vl, res.accX_s_vl, 'r--', 'LineWidth', 1.3);
    legend('IMU','Sim'); ylabel('AccX (m/s²)');
    title(sprintf('%s — AccX (R^2=%.3f)', lbl_vl, R2_func(acc_vl(:,1), res.accX_s_vl))); grid on;
    subplot(3,1,2);
    plot(time_vl, acc_vl(:,2), 'b-', time_vl, res.accY_s_vl, 'r--', 'LineWidth', 1.3);
    legend('IMU','Sim'); ylabel('AccY (m/s²)');
    title(sprintf('%s — AccY (R^2=%.3f)', lbl_vl, R2_func(acc_vl(:,2), res.accY_s_vl))); grid on;
    subplot(3,1,3);
    plot(time_vl, acc_vl(:,3), 'b-', time_vl, res.accZ_s_vl, 'r--', 'LineWidth', 1.3);
    legend('IMU','Sim'); xlabel('Tempo (s)'); ylabel('AccZ (m/s²)');
    title(sprintf('%s — AccZ (R^2=%.3f)', lbl_vl, R2_func(acc_vl(:,3), res.accZ_s_vl))); grid on;
    sgtitle(['Acelerações (Val) — ' titulo]);
    saveas(fig2, fullfile('C:/Users/Henrique/ARTIGO/identification/images', ['acc_val_' safe_titulo '.png']));

    % --- Gráfico de Atitude (phi, theta, psi) ---
    if nargin >= 12 && ~isempty(att_vl) && isfield(res, 'phi_s_vl') && res.pqr_vl_ok
        fig3 = figure('Name', ['Atitude - ' titulo], 'Position', [160 130 900 700], 'Visible', 'off');
        subplot(3,1,1);
        plot(time_vl, att_vl(:,1), 'b-', time_vl, res.phi_s_vl, 'r--', 'LineWidth', 1.3);
        legend('EKF','Sim'); ylabel('\phi (°)');
        title(sprintf('%s — \\phi (R^2=%.3f)', lbl_vl, R2_func(att_vl(:,1), res.phi_s_vl))); grid on;
        subplot(3,1,2);
        plot(time_vl, att_vl(:,2), 'b-', time_vl, res.theta_s_vl, 'r--', 'LineWidth', 1.3);
        legend('EKF','Sim'); ylabel('\theta (°)');
        title(sprintf('%s — \\theta (R^2=%.3f)', lbl_vl, R2_func(att_vl(:,2), res.theta_s_vl))); grid on;
        subplot(3,1,3);
        plot(time_vl, att_vl(:,3), 'b-', time_vl, res.psi_s_vl, 'r--', 'LineWidth', 1.3);
        legend('EKF','Sim'); xlabel('Tempo (s)'); ylabel('\psi (°)');
        title(sprintf('%s — \\psi (R^2=%.3f)', lbl_vl, R2_func(att_vl(:,3), res.psi_s_vl))); grid on;
        sgtitle(['Atitude (Val) — ' titulo]);
        saveas(fig3, fullfile('C:/Users/Henrique/ARTIGO/identification/images', ['att_val_' safe_titulo '.png']));
    end
end

function plot_torques(P, time_vl, pwm_vl, func_T_ref, func_Q_ref, pqr_vl, t_val)
    k_T = P(5:8);
    k_Q = P(9:12);
    dx_cg = P(19); dy_cg = P(20);
    N_vl = length(time_vl);

    T_ref_vl = zeros(N_vl, 4);
    Q_ref_vl = zeros(N_vl, 4);
    for j = 1:4
        T_ref_vl(:,j) = func_T_ref(pwm_vl(:,j));
        Q_ref_vl(:,j) = func_Q_ref(pwm_vl(:,j));
    end

    Tmr = T_ref_vl .* k_T';
    Qmr = Q_ref_vl .* k_Q';

    % Braços efetivos com offset do CG
    Lx_r = 0.232 - dy_cg;   % direita (motores 1,4)
    Lx_l = 0.232 + dy_cg;   % esquerda (motores 2,3)
    Ly_f = 0.311185 - dx_cg; % frente (motores 1,3)
    Ly_r = 0.342865 + dx_cg; % traseira (motores 2,4)

    Mx = -(Lx_r*Tmr(:,1) - Lx_l*Tmr(:,2) - Lx_l*Tmr(:,3) + Lx_r*Tmr(:,4));
    My = Ly_f*Tmr(:,1) - Ly_r*Tmr(:,2) + Ly_f*Tmr(:,3) - Ly_r*Tmr(:,4);
    Mz = Qmr(:,1) + Qmr(:,2) - Qmr(:,3) - Qmr(:,4);

    fig3 = figure('Name', 'Torques Validação', 'Position', [80 50 1000 900], 'Visible', 'off');

    subplot(4,1,1);
    plot(time_vl, Tmr(:,1), time_vl, Tmr(:,2), time_vl, Tmr(:,3), time_vl, Tmr(:,4), 'LineWidth', 1.2);
    legend('T_1','T_2','T_3','T_4'); ylabel('Empuxo (N)'); grid on;
    title(sprintf('Empuxo individual (k_T=[%.3f, %.3f, %.3f, %.3f])', k_T(1), k_T(2), k_T(3), k_T(4)));

    subplot(4,1,2);
    plot(time_vl, Mx, 'r-', 'LineWidth', 1.3);
    hold on; plot(time_vl, pqr_vl(:,1)*0.1, 'b--', 'LineWidth', 1); hold off;
    legend('Mx (N·m)', 'p×0.1 (ref)'); ylabel('Mx (N·m)'); grid on;
    title(sprintf('Roll: Mx  (Lx_r=%.3f, Lx_l=%.3f, dy_{cg}=%.4f)', Lx_r, Lx_l, dy_cg));

    subplot(4,1,3);
    plot(time_vl, My, 'r-', 'LineWidth', 1.3);
    hold on; plot(time_vl, pqr_vl(:,2)*0.1, 'b--', 'LineWidth', 1); hold off;
    legend('My (N·m)', 'q×0.1 (ref)'); ylabel('My (N·m)'); grid on;
    title(sprintf('Pitch: My  (Ly_f=%.3f, Ly_r=%.3f, dx_{cg}=%.4f)', Ly_f, Ly_r, dx_cg));

    subplot(4,1,4);
    plot(time_vl, Mz, 'r-', 'LineWidth', 1.3);
    hold on; plot(time_vl, pqr_vl(:,3)*0.1, 'b--', 'LineWidth', 1); hold off;
    legend('Mz (N·m)', 'r×0.1 (ref)'); ylabel('Mz (N·m)'); xlabel('Tempo (s)'); grid on;
    title(sprintf('Yaw: Mz  (k_Q=[%.3f, %.3f, %.3f, %.3f])', k_Q(1), k_Q(2), k_Q(3), k_Q(4)));

    sgtitle(sprintf('Torques na Validação [%d-%ds]  (dx_{cg}=%.4f, dy_{cg}=%.4f)', t_val(1), t_val(2), dx_cg, dy_cg));
    saveas(fig3, fullfile('C:/Users/Henrique/ARTIGO/identification/images', 'torques_validacao.png'));

    fprintf('\n  --- Torques (validação) | dx_cg=%.4f dy_cg=%.4f ---\n', dx_cg, dy_cg);
    fprintf('    Braços: Lx_r=%.4f Lx_l=%.4f Ly_f=%.4f Ly_r=%.4f\n', Lx_r, Lx_l, Ly_f, Ly_r);
    fprintf('    Mx: min=%.4f  max=%.4f  std=%.4f\n', min(Mx), max(Mx), std(Mx));
    fprintf('    My: min=%.4f  max=%.4f  std=%.4f\n', min(My), max(My), std(My));
    fprintf('    Mz: min=%.4f  max=%.4f  std=%.4f\n', min(Mz), max(Mz), std(Mz));
end

function plot_forces(P, time_vl, pwm_vl, att_vl, acc_vl, func_T_ref, constants_sim, t_val, res)
    k_T = P(5:8);
    Xu_m = P(21); Yv_m = P(22); Zw_m = P(23); Bz = P(24);
    m = constants_sim.m;
    g = constants_sim.g;
    N_vl = length(time_vl);

    % Atitude em radianos
    phi   = deg2rad(att_vl(:,1));
    theta = deg2rad(att_vl(:,2));

    % Empuxo individual e total
    Ti = zeros(N_vl, 4);
    for j = 1:4
        Ti(:,j) = k_T(j) * func_T_ref(pwm_vl(:,j));
    end
    T_total = sum(Ti, 2);

    % Componentes da gravidade no body frame (por unidade de massa)
    gx = -g * sin(theta);
    gy =  g * cos(theta) .* sin(phi);
    gz =  g * cos(theta) .* cos(phi);

    % Forças por unidade de massa (acelerações)
    Fz_m = -T_total/m + gz;       % thrust + gravidade z

    % =====================================================================
    %  Figura 1: Forças translacionais com modelo simulado
    % =====================================================================
    fig4 = figure('Name', 'Forças Translacionais', 'Position', [80 50 1000 900], 'Visible', 'off');

    % --- Subplot 1: Empuxo total vs peso ---
    subplot(4,1,1);
    plot(time_vl, T_total, 'r-', 'LineWidth', 1.3);
    hold on;
    yline(m*g, 'k--', 'LineWidth', 1.2);
    hold off;
    legend('T_{total}', sprintf('mg = %.2f N', m*g));
    ylabel('Força (N)'); grid on;
    title(sprintf('Empuxo total vs Peso  (T_{med}=%.2f N, mg=%.2f N, ratio=%.3f)', ...
        mean(T_total), m*g, mean(T_total)/(m*g)));

    % --- Subplot 2: AccX medido vs simulado ---
    subplot(4,1,2);
    plot(time_vl, acc_vl(:,1), 'b-', 'LineWidth', 1.0);
    hold on;
    if res.full_vl_ok
        plot(time_vl, res.accX_s_vl, 'r--', 'LineWidth', 1.3);
    end
    plot(time_vl, gx, 'k:', 'LineWidth', 0.8);
    hold off;
    if res.full_vl_ok
        legend('AccX_{IMU} (medido)', 'AccX_{modelo}', 'gx = -g sin\theta');
    else
        legend('AccX_{IMU} (medido)', 'gx = -g sin\theta');
    end
    ylabel('m/s^2'); grid on;
    title('Eixo X: força específica medida vs modelo');

    % --- Subplot 3: AccY medido vs simulado ---
    subplot(4,1,3);
    plot(time_vl, acc_vl(:,2), 'b-', 'LineWidth', 1.0);
    hold on;
    if res.full_vl_ok
        plot(time_vl, res.accY_s_vl, 'r--', 'LineWidth', 1.3);
    end
    plot(time_vl, gy, 'k:', 'LineWidth', 0.8);
    hold off;
    if res.full_vl_ok
        legend('AccY_{IMU} (medido)', 'AccY_{modelo}', 'gy = g cos\theta sin\phi');
    else
        legend('AccY_{IMU} (medido)', 'gy = g cos\theta sin\phi');
    end
    ylabel('m/s^2'); grid on;
    title('Eixo Y: força específica medida vs modelo');

    % --- Subplot 4: AccZ medido vs simulado ---
    subplot(4,1,4);
    plot(time_vl, acc_vl(:,3), 'b-', 'LineWidth', 1.0);
    hold on;
    if res.full_vl_ok
        plot(time_vl, res.accZ_s_vl, 'r--', 'LineWidth', 1.3);
    end
    plot(time_vl, Fz_m, 'k:', 'LineWidth', 0.8);
    yline(-g, 'k--', 'LineWidth', 0.6);
    hold off;
    if res.full_vl_ok
        legend('AccZ_{IMU} (medido)', 'AccZ_{modelo}', 'Fz/m = -T/m+gz', '-g');
    else
        legend('AccZ_{IMU} (medido)', 'Fz/m = -T/m+gz', '-g');
    end
    ylabel('m/s^2'); xlabel('Tempo (s)'); grid on;
    title(sprintf('Eixo Z: força específica medida vs modelo  (Bz=%.3f, Zw=%.3f)', Bz, Zw_m));

    sgtitle(sprintf('Forças Translacionais — Validação [%d-%ds]  (m=%.1f kg)', t_val(1), t_val(2), m));
    saveas(fig4, fullfile('C:/Users/Henrique/ARTIGO/identification/images', 'forcas_validacao.png'));

    % =====================================================================
    %  Figura 2: Análise dos PWMs por motor
    % =====================================================================
    fig5 = figure('Name', 'Análise PWM', 'Position', [100 50 1000 700], 'Visible', 'off');

    % Subplot 1: PWM de cada motor ao longo do tempo
    subplot(3,1,1);
    plot(time_vl, pwm_vl(:,1), 'r-', 'LineWidth', 1.0); hold on;
    plot(time_vl, pwm_vl(:,2), 'b-', 'LineWidth', 1.0);
    plot(time_vl, pwm_vl(:,3), 'g-', 'LineWidth', 1.0);
    plot(time_vl, pwm_vl(:,4), 'm-', 'LineWidth', 1.0);
    hold off;
    legend('M1 (FR,CW)', 'M2 (RL,CW)', 'M3 (FL,CCW)', 'M4 (RR,CCW)');
    ylabel('PWM (\mus)'); grid on;
    title('PWM por motor');

    % Subplot 2: Empuxo individual (com k_T)
    subplot(3,1,2);
    plot(time_vl, Ti(:,1), 'r-', 'LineWidth', 1.0); hold on;
    plot(time_vl, Ti(:,2), 'b-', 'LineWidth', 1.0);
    plot(time_vl, Ti(:,3), 'g-', 'LineWidth', 1.0);
    plot(time_vl, Ti(:,4), 'm-', 'LineWidth', 1.0);
    yline(m*g/4, 'k--', 'LineWidth', 1.0);
    hold off;
    legend(sprintf('T1 (k_T=%.3f)', k_T(1)), sprintf('T2 (k_T=%.3f)', k_T(2)), ...
           sprintf('T3 (k_T=%.3f)', k_T(3)), sprintf('T4 (k_T=%.3f)', k_T(4)), 'mg/4');
    ylabel('Empuxo (N)'); grid on;
    title('Empuxo por motor (T_i = k_{Ti} \cdot T_{ref}(PWM_i))');

    % Subplot 3: T_ref cru (sem k_T) — para ver o que o FC realmente comanda
    T_ref_raw = zeros(N_vl, 4);
    for j = 1:4
        T_ref_raw(:,j) = func_T_ref(pwm_vl(:,j));
    end
    subplot(3,1,3);
    plot(time_vl, T_ref_raw(:,1), 'r-', 'LineWidth', 1.0); hold on;
    plot(time_vl, T_ref_raw(:,2), 'b-', 'LineWidth', 1.0);
    plot(time_vl, T_ref_raw(:,3), 'g-', 'LineWidth', 1.0);
    plot(time_vl, T_ref_raw(:,4), 'm-', 'LineWidth', 1.0);
    yline(m*g/4, 'k--', 'LineWidth', 1.0);
    hold off;
    T_ref_total = sum(T_ref_raw, 2);
    legend(sprintf('T_{ref1} (med=%.2f)', mean(T_ref_raw(:,1))), ...
           sprintf('T_{ref2} (med=%.2f)', mean(T_ref_raw(:,2))), ...
           sprintf('T_{ref3} (med=%.2f)', mean(T_ref_raw(:,3))), ...
           sprintf('T_{ref4} (med=%.2f)', mean(T_ref_raw(:,4))), 'mg/4');
    ylabel('Empuxo ref (N)'); xlabel('Tempo (s)'); grid on;
    title(sprintf('T_{ref}(PWM) cru (sem k_T) — Total medio=%.2f N vs mg=%.2f N, ratio=%.3f', ...
        mean(T_ref_total), m*g, mean(T_ref_total)/(m*g)));

    sgtitle(sprintf('Análise PWM e Empuxo — Validação [%d-%ds]', t_val(1), t_val(2)));
    saveas(fig5, fullfile('C:/Users/Henrique/ARTIGO/identification/images', 'pwm_analise.png'));

    % Estatísticas
    fprintf('\n  --- Forças (validação) ---\n');
    fprintf('    T_total: média=%.3f N  (mg=%.3f N, ratio=%.4f)\n', mean(T_total), m*g, mean(T_total)/(m*g));
    fprintf('    Ti médios: [%.2f, %.2f, %.2f, %.2f] N\n', mean(Ti(:,1)), mean(Ti(:,2)), mean(Ti(:,3)), mean(Ti(:,4)));
    fprintf('    T_ref cru (sem k_T): [%.2f, %.2f, %.2f, %.2f] N  (total=%.2f)\n', ...
        mean(T_ref_raw(:,1)), mean(T_ref_raw(:,2)), mean(T_ref_raw(:,3)), mean(T_ref_raw(:,4)), mean(T_ref_total));
    fprintf('    PWM médios: [%.0f, %.0f, %.0f, %.0f] us\n', mean(pwm_vl(:,1)), mean(pwm_vl(:,2)), mean(pwm_vl(:,3)), mean(pwm_vl(:,4)));
    fprintf('    PWM std:    [%.1f, %.1f, %.1f, %.1f] us\n', std(pwm_vl(:,1)), std(pwm_vl(:,2)), std(pwm_vl(:,3)), std(pwm_vl(:,4)));
    fprintf('    Fz/m média=%.4f  (esperado ~0 em hover: -T/m+gz ~0)\n', mean(Fz_m));
    fprintf('    Xu_m=%.4f  Yv_m=%.4f  Zw_m=%.4f  Bz=%.4f\n', Xu_m, Yv_m, Zw_m, Bz);

    % Análise: assimetria constante (CG/motor) vs variável (vento)
    fprintf('\n  --- Análise de assimetria PWM ---\n');
    pwm_mean = mean(pwm_vl);
    pwm_std  = std(pwm_vl);
    cv = pwm_std ./ pwm_mean * 100;  % coeficiente de variação (%)
    fprintf('    Coeficiente de variação PWM: [%.1f%%, %.1f%%, %.1f%%, %.1f%%]\n', cv(1), cv(2), cv(3), cv(4));
    fprintf('    (CV baixo = patamar constante → diferença CG/motor)\n');
    fprintf('    (CV alto  = variação temporal → vento/perturbação)\n');

    % Assimetria roll: (M1+M4) vs (M2+M3)
    pwm_right = pwm_vl(:,1) + pwm_vl(:,4);  % motores direita
    pwm_left  = pwm_vl(:,2) + pwm_vl(:,3);  % motores esquerda
    roll_asym = pwm_right - pwm_left;
    fprintf('    Assimetria roll (dir-esq): média=%.1f us  std=%.1f us\n', mean(roll_asym), std(roll_asym));

    % Assimetria pitch: (M1+M3) vs (M2+M4)
    pwm_front = pwm_vl(:,1) + pwm_vl(:,3);  % motores frente
    pwm_rear  = pwm_vl(:,2) + pwm_vl(:,4);  % motores traseira
    pitch_asym = pwm_front - pwm_rear;
    fprintf('    Assimetria pitch (frente-trás): média=%.1f us  std=%.1f us\n', mean(pitch_asym), std(pitch_asym));

    % Assimetria yaw: (M1+M2) CW vs (M3+M4) CCW
    pwm_cw  = pwm_vl(:,1) + pwm_vl(:,2);
    pwm_ccw = pwm_vl(:,3) + pwm_vl(:,4);
    yaw_asym = pwm_cw - pwm_ccw;
    fprintf('    Assimetria yaw (CW-CCW): média=%.1f us  std=%.1f us\n', mean(yaw_asym), std(yaw_asym));

    if mean(abs(std(roll_asym))) < 5 && mean(abs(std(pitch_asym))) < 5
        fprintf('    → Assimetrias ESTÁVEIS: predomina CG offset / diferença de motores\n');
    else
        fprintf('    → Assimetrias VARIÁVEIS: predomina vento / perturbação externa\n');
    end
end
