%% run_identification.m — Identificacao do modelo VTOL (analysis/)
%
% Pipeline autocontido: EEM (Fase A) + OEM (Fase B)
%
% Vetor de parametros (24):
%   P(1:6)   = [Jx, Jy, Jz, Jxy, Jxz, Jyz]  (momentos de inercia)
%   P(7:10)  = k_T1..k_T4   (fatores de escala empuxo)
%   P(11:14) = k_Q1..k_Q4   (fatores de escala torque)
%   P(15:17) = Dp, Dq, Dr   (amortecimento rotacional)
%   P(18:20) = Bp, Bq, Br   (bias rotacional)
%   P(21:24) = Xu, Yv, Zw, Bz (arrasto translacional + bias vertical)
%
% Fase A: EEM rotacional (algebraico) -> P(1:20)
% Fase B: OEM progressivo (RK4 windowed) -> P(1:24)

clear; close all; clc;

%% ========================================================================
%  1. CONFIGURACAO
% =========================================================================
t_trains = {[157, 167]; [167, 177]; [177, 187]};
t_val    = [147, 157];

dt = 0.1;
smooth_win = 5;

R2_func = @(y_e, y_s) 1 - sum((y_e - y_s).^2) / max(sum((y_e - mean(y_e)).^2), 1e-12);

img_dir = fullfile(fileparts(mfilename('fullpath')), 'images');
if ~exist(img_dir, 'dir'), mkdir(img_dir); end

%% ========================================================================
%  2. CARREGAR DADOS E INTERPOLAR
% =========================================================================
fprintf('=== Carregando dados ===\n');
load(fullfile('..', 'identification', 'log_data.mat'));

IMU.TimeS  = double(IMU.TimeUS) / 1e6;
ATT.TimeS  = double(ATT.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;

idx_imu = IMU.I == 0;
time_IMU  = IMU.TimeS(idx_imu);
time_ATT  = ATT.TimeS;
time_RCOU = RCOU.TimeS;

t_start = max([min(time_IMU), min(time_ATT), min(time_RCOU)]);
t_end   = min([max(time_IMU), max(time_ATT), max(time_RCOU)]);
t_common = (t_start:dt:t_end)';

gyrX = interp1(time_IMU, IMU.GyrX(idx_imu), t_common, 'linear');
gyrY = interp1(time_IMU, IMU.GyrY(idx_imu), t_common, 'linear');
gyrZ = interp1(time_IMU, IMU.GyrZ(idx_imu), t_common, 'linear');
accX = interp1(time_IMU, IMU.AccX(idx_imu), t_common, 'linear');
accY = interp1(time_IMU, IMU.AccY(idx_imu), t_common, 'linear');
accZ = interp1(time_IMU, IMU.AccZ(idx_imu), t_common, 'linear');
roll_d  = interp1(time_ATT, ATT.Roll,  t_common, 'linear');
pitch_d = interp1(time_ATT, ATT.Pitch, t_common, 'linear');
yaw_d   = interp1(time_ATT, ATT.Yaw,   t_common, 'linear');
pwm1 = interp1(time_RCOU, double(RCOU.C1), t_common, 'linear');
pwm2 = interp1(time_RCOU, double(RCOU.C2), t_common, 'linear');
pwm3 = interp1(time_RCOU, double(RCOU.C3), t_common, 'linear');
pwm4 = interp1(time_RCOU, double(RCOU.C4), t_common, 'linear');

fprintf('  Dados interpolados: %d amostras, dt=%.1fs\n', length(t_common), dt);

%% ========================================================================
%  3. PARAMETROS DO MODELO (CAD)
% =========================================================================
params = setup_params();
fprintf('  Massa: %.5f kg\n', params.mass);
fprintf('  J (CAD): Jx=%.3e  Jy=%.3e  Jz=%.3e\n', params.Jx, params.Jy, params.Jz);
fprintf('           Jxy=%.3e Jxz=%.3e Jyz=%.3e\n', params.Jxy, params.Jxz, params.Jyz);

%% ========================================================================
%  4. EXTRAIR SEGMENTOS DE TREINO
% =========================================================================
n_seg = length(t_trains);
segs = cell(n_seg, 1);

fprintf('\n  Segmentos de treino: %d\n', n_seg);
for s = 1:n_seg
    idx_s = (t_common >= t_trains{s}(1)) & (t_common <= t_trains{s}(2));
    seg.time = t_common(idx_s);
    seg.N    = sum(idx_s);
    seg.pwm  = [pwm1(idx_s), pwm2(idx_s), pwm3(idx_s), pwm4(idx_s)];
    seg.pqr  = [gyrX(idx_s), gyrY(idx_s), gyrZ(idx_s)];
    seg.acc  = [accX(idx_s), accY(idx_s), accZ(idx_s)];
    seg.att  = [roll_d(idx_s), pitch_d(idx_s), yaw_d(idx_s)];
    seg.att_rad = deg2rad(seg.att);
    seg.T_ref = zeros(seg.N, 4);
    seg.Q_ref = zeros(seg.N, 4);
    for j = 1:4
        seg.T_ref(:,j) = params.T_ref(seg.pwm(:,j));
        seg.Q_ref(:,j) = params.Q_ref(seg.pwm(:,j));
    end
    segs{s} = seg;
    fprintf('    [%d] %d-%ds  (%d pontos)\n', s, t_trains{s}(1), t_trains{s}(2), seg.N);
end

% Dados de validacao
idx_val = (t_common >= t_val(1)) & (t_common <= t_val(2));
time_vl = t_common(idx_val);
pwm_vl  = [pwm1(idx_val), pwm2(idx_val), pwm3(idx_val), pwm4(idx_val)];
pqr_vl  = [gyrX(idx_val), gyrY(idx_val), gyrZ(idx_val)];
acc_vl  = [accX(idx_val), accY(idx_val), accZ(idx_val)];
att_vl  = [roll_d(idx_val), pitch_d(idx_val), yaw_d(idx_val)];

%% ========================================================================
%  5. CHUTE INICIAL E LIMITES (24 parametros)
% =========================================================================
n_params = 24;
n_rot = 20;  % EEM otimiza P(1:20) — inercias + rotacional

% Chute inicial: CAD para inercias, 1.0 para k, 0 para damping/bias
P0 = [params.Jx; params.Jy; params.Jz;        ... % Jx, Jy, Jz
      params.Jxy; params.Jxz; params.Jyz;      ... % Jxy, Jxz, Jyz
      1.0; 1.0; 1.0; 1.0;                      ... % k_T
      1.0; 1.0; 1.0; 1.0;                      ... % k_Q
      0.0; 0.0; 0.0;                           ... % Dp, Dq, Dr
      0.0; 0.0; 0.0;                           ... % Bp, Bq, Br
      0.0; 0.0; 0.0; 0.0];                          % Xu, Yv, Zw, Bz

% Limites: inercias +-200% do CAD, k_T/k_Q [0.1, 5], etc.
lb = [params.Jx*0.2;  params.Jy*0.2;  params.Jz*0.2;   ... % Jx, Jy, Jz
      -0.010; -0.005; -0.005;                            ... % Jxy, Jxz, Jyz
      0.1; 0.1; 0.1; 0.1;                               ... % k_T
      0.1; 0.1; 0.1; 0.1;                               ... % k_Q
        0;   0;   0;                                     ... % Dp >= 0
      -10; -10; -10;                                     ... % Bp
      -30; -30; -5; -5];                                      % Xu, Yv, Zw, Bz

ub = [params.Jx*5.0;  params.Jy*5.0;  params.Jz*5.0;   ... % Jx, Jy, Jz
       0.005;  0.010;  0.005;                            ... % Jxy, Jxz, Jyz
      5.0; 5.0; 5.0; 5.0;                               ... % k_T
      5.0; 5.0; 5.0; 5.0;                               ... % k_Q
       50;  50;  50;                                     ... % Dp
       10;  10;  10;                                     ... % Bp
        0;   0;   0;  5];                                     % Xu, Yv, Zw, Bz

param_names = {'Jx','Jy','Jz','Jxy','Jxz','Jyz', ...
    'k_T1','k_T2','k_T3','k_T4','k_Q1','k_Q2','k_Q3','k_Q4', ...
    'Dp','Dq','Dr','Bp','Bq','Br', ...
    'Xu','Yv','Zw','Bz'};

%% ========================================================================
%  6. FASE A: EEM ROTACIONAL (20 params: inercias + k + D + B)
% =========================================================================
fprintf('\n==========================================================\n');
fprintf('  FASE A: EEM rotacional (%d segmentos, %d params)\n', n_seg, n_rot);
fprintf('==========================================================\n');

p_all = []; q_all = []; r_all = [];
pd_all = []; qd_all = []; rd_all = [];
Tr_all = []; Qr_all = [];

for s = 1:n_seg
    sg = segs{s};
    p_s = sg.pqr(:,1); q_s = sg.pqr(:,2); r_s = sg.pqr(:,3);
    p_all  = [p_all; p_s]; %#ok<AGROW>
    q_all  = [q_all; q_s]; %#ok<AGROW>
    r_all  = [r_all; r_s]; %#ok<AGROW>
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

cost_eem = @(Prot) eem_cost(Prot, params, p_all, q_all, r_all, ...
    pd_all, qd_all, rd_all, Tr_all, Qr_all, weights_eem);

opts_eem = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter', ...
    'MaxIterations', 2000, ...
    'MaxFunctionEvaluations', 80000, ...
    'StepTolerance', 1e-14, ...
    'FunctionTolerance', 1e-14);

fprintf('  Otimizando %d parametros (inercias + rotacional)...\n', n_rot);
[P_eem_rot, rn_eem, ~, ef_eem] = lsqnonlin(cost_eem, P0(1:n_rot), lb(1:n_rot), ub(1:n_rot), opts_eem);

P_eem = [P_eem_rot; P0(n_rot+1:n_params)];

fprintf('\n  EEM Resnorm: %.4f  |  Exit flag: %d\n', rn_eem, ef_eem);
fprintf('\n  --- Parametros EEM ---\n');
for i = 1:n_rot
    fprintf('    %-6s: %12.6f  (CAD/chute: %12.6f)\n', param_names{i}, P_eem(i), P0(i));
end

% Comparar inercias
fprintf('\n  Inercias identificadas vs CAD:\n');
fprintf('    Jx:  %.6f -> %.6f  (ratio: %.2f)\n', params.Jx, P_eem(1), P_eem(1)/params.Jx);
fprintf('    Jy:  %.6f -> %.6f  (ratio: %.2f)\n', params.Jy, P_eem(2), P_eem(2)/params.Jy);
fprintf('    Jz:  %.6f -> %.6f  (ratio: %.2f)\n', params.Jz, P_eem(3), P_eem(3)/params.Jz);

% R2 do EEM
e_eem = cost_eem(P_eem_rot);
N_eem = length(p_all);
pd_pred = pd_all - e_eem(1:N_eem) / sqrt(weights_eem(1));
qd_pred = qd_all - e_eem(N_eem+1:2*N_eem) / sqrt(weights_eem(2));
rd_pred = rd_all - e_eem(2*N_eem+1:end) / sqrt(weights_eem(3));
fprintf('\n  R2 EEM (treino concatenado):\n');
fprintf('    p_dot: %.4f\n', R2_func(pd_all, pd_pred));
fprintf('    q_dot: %.4f\n', R2_func(qd_all, qd_pred));
fprintf('    r_dot: %.4f\n', R2_func(rd_all, rd_pred));

%% ========================================================================
%  7. FASE B: OEM PROGRESSIVO (janelas 1s -> 2s -> 3s)
% =========================================================================
fprintf('\n==========================================================\n');
fprintf('  FASE B: OEM Progressivo (24 params)\n');
fprintf('==========================================================\n');

pqr_cat = cell2mat(cellfun(@(s) s.pqr, segs, 'UniformOutput', false));
acc_cat = cell2mat(cellfun(@(s) s.acc, segs, 'UniformOutput', false));

weights_pqr = [1/max(var(pqr_cat(:,1)),1e-12); ...
               1/max(var(pqr_cat(:,2)),1e-12); ...
               1/max(var(pqr_cat(:,3)),1e-12)];

weights_acc = [1/max(var(acc_cat(:,1)),1e-12); ...
               1/max(var(acc_cat(:,2)),1e-12); ...
               1/max(var(acc_cat(:,3)),1e-12)];

win_durations = [1.0, 2.0, 3.0];
P_current = P_eem;
best_R2_mean = -Inf;
P_best = P_eem;
best_stage_name = 'EEM';

for stage = 1:length(win_durations)
    win_sec = win_durations(stage);
    stage_name = sprintf('%.0fs', win_sec);

    fprintf('\n  --- OEM Stage [%s] (%d segmentos) ---\n', stage_name, n_seg);

    cost_oem = @(P) oem_multi_seg_cost(P, segs, win_sec, params, dt, ...
        weights_pqr, weights_acc);

    opts_oem = optimoptions('lsqnonlin', ...
        'Algorithm', 'trust-region-reflective', ...
        'Display', 'iter', ...
        'MaxIterations', 500, ...
        'MaxFunctionEvaluations', 80000, ...
        'StepTolerance', 1e-14, ...
        'FunctionTolerance', 1e-14);

    [P_current, rn, ~, ef] = lsqnonlin(cost_oem, P_current, lb, ub, opts_oem);
    fprintf('  [%s] Resnorm: %.4f | Exit: %d\n', stage_name, rn, ef);

    [R2_val_p, R2_val_q, R2_val_r] = evaluate_rotational(P_current, params, ...
        time_vl, pwm_vl, pqr_vl, dt);

    R2_mean_stage = mean([R2_val_p, R2_val_q, R2_val_r]);
    fprintf('  [%s] Val R2: p=%.4f | q=%.4f | r=%.4f | media=%.4f\n', ...
        stage_name, R2_val_p, R2_val_q, R2_val_r, R2_mean_stage);

    if R2_mean_stage > best_R2_mean
        best_R2_mean = R2_mean_stage;
        P_best = P_current;
        best_stage_name = stage_name;
    end
end

fprintf('\n  >>> Melhor estagio: [%s] com R2 medio = %.4f\n', best_stage_name, best_R2_mean);
P_final = P_best;

%% ========================================================================
%  8. RESULTADOS FINAIS
% =========================================================================
fprintf('\n==========================================================\n');
fprintf('  PARAMETROS FINAIS\n');
fprintf('==========================================================\n');
for i = 1:n_params
    fprintf('  %-6s: %12.6f  (EEM: %12.6f)  (P0: %12.6f)\n', ...
        param_names{i}, P_final(i), P_eem(i), P0(i));
end

% Montar params identificado
params_id = params;
params_id.Jx = P_final(1); params_id.Jy = P_final(2); params_id.Jz = P_final(3);
params_id.Jxy = P_final(4); params_id.Jxz = P_final(5); params_id.Jyz = P_final(6);
params_id.J = [P_final(1), P_final(4), P_final(5);
               P_final(4), P_final(2), P_final(6);
               P_final(5), P_final(6), P_final(3)];
params_id.J_inv = inv(params_id.J);
params_id.k_T = P_final(7:10);
params_id.k_Q = P_final(11:14);
params_id.Dp = P_final(15); params_id.Dq = P_final(16); params_id.Dr = P_final(17);
params_id.Bp = P_final(18); params_id.Bq = P_final(19); params_id.Br = P_final(20);
params_id.Xu = P_final(21); params_id.Yv = P_final(22);
params_id.Zw = P_final(23); params_id.Bz = P_final(24);

save(fullfile(fileparts(mfilename('fullpath')), 'P_identified.mat'), ...
    'P_final', 'P_eem', 'P0', 'param_names', 'params_id');
fprintf('\n  Parametros salvos em P_identified.mat\n');

%% ========================================================================
%  9. VALIDACAO COM ODE45 (simulacao completa de 9 estados)
% =========================================================================
fprintf('\n==========================================================\n');
fprintf('  VALIDACAO COMPLETA (ode45, 9 estados)\n');
fprintf('==========================================================\n');

ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% --- Treino (primeiro segmento) ---
sg = segs{1};
att_rad0 = deg2rad(sg.att(1,:));
x0_tr = [sg.pqr(1,:)'; att_rad0(:); 0; 0; 0];

try
    [t_s, x_s] = ode45(@(tt,xx) vtol_dynamics_legacy(tt, xx, sg.time, sg.pwm, params_id), ...
        sg.time, x0_tr, ode_opts);
    x_tr = interp1(t_s, x_s, sg.time, 'linear', 'extrap');
    tr_ok = true;
    fprintf('  Treino [%d-%ds]: OK (%d passos)\n', t_trains{1}(1), t_trains{1}(2), length(t_s));
    fprintf('    R2 p=%.4f | q=%.4f | r=%.4f\n', ...
        R2_func(sg.pqr(:,1), x_tr(:,1)), ...
        R2_func(sg.pqr(:,2), x_tr(:,2)), ...
        R2_func(sg.pqr(:,3), x_tr(:,3)));
catch ME
    fprintf('  Treino: DIVERGIU (%s)\n', ME.message);
    tr_ok = false;
end

% --- Validacao ---
att_rad0_vl = deg2rad(att_vl(1,:));
x0_vl = [pqr_vl(1,:)'; att_rad0_vl(:); 0; 0; 0];

try
    [t_sv, x_sv] = ode45(@(tt,xx) vtol_dynamics_legacy(tt, xx, time_vl, pwm_vl, params_id), ...
        time_vl, x0_vl, ode_opts);
    x_vl = interp1(t_sv, x_sv, time_vl, 'linear', 'extrap');
    vl_ok = true;
    fprintf('  Validacao [%d-%ds]: OK (%d passos)\n', t_val(1), t_val(2), length(t_sv));
    fprintf('    R2 p=%.4f | q=%.4f | r=%.4f\n', ...
        R2_func(pqr_vl(:,1), x_vl(:,1)), ...
        R2_func(pqr_vl(:,2), x_vl(:,2)), ...
        R2_func(pqr_vl(:,3), x_vl(:,3)));
    fprintf('    R2 phi=%.4f | theta=%.4f | psi=%.4f\n', ...
        R2_func(att_vl(:,1), rad2deg(x_vl(:,4))), ...
        R2_func(att_vl(:,2), rad2deg(x_vl(:,5))), ...
        R2_func(att_vl(:,3), rad2deg(x_vl(:,6))));
catch ME
    fprintf('  Validacao: DIVERGIU (%s)\n', ME.message);
    vl_ok = false;
end

%% ========================================================================
%  10. ANALISE DE ACELERACOES (forca especifica semi-acoplada)
% =========================================================================
fprintf('\n==========================================================\n');
fprintf('  ANALISE DE ACELERACOES (semi-acoplada)\n');
fprintf('==========================================================\n');

% Semi-acoplado: usa pqr e atitude MEDIDOS, integra u,v,w, compara com IMU
datasets = {{'treino', sg.time, sg.pwm, sg.pqr, sg.acc, sg.att, t_trains{1}}, ...
            {'validacao', time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, t_val}};

for d = 1:2
    ds_name = datasets{d}{1};
    ds_time = datasets{d}{2};
    ds_pwm  = datasets{d}{3};
    ds_pqr  = datasets{d}{4};
    ds_acc  = datasets{d}{5};
    ds_att  = datasets{d}{6};
    ds_range = datasets{d}{7};
    Nd = length(ds_time);

    % Empuxo total
    T_tot_d = zeros(Nd, 1);
    for k = 1:Nd
        for j = 1:4
            T_tot_d(k) = T_tot_d(k) + P_final(6+j) * params.T_ref(ds_pwm(k,j));
        end
    end

    % Gravidade no corpo
    att_rad_d = deg2rad(ds_att);
    gx_d = -params.g * sin(att_rad_d(:,2));
    gy_d =  params.g * cos(att_rad_d(:,2)) .* sin(att_rad_d(:,1));
    gz_d =  params.g * cos(att_rad_d(:,2)) .* cos(att_rad_d(:,1));

    Xu_f = P_final(21); Yv_f = P_final(22); Zw_f = P_final(23); Bz_f = P_final(24);

    % Integrar u,v,w com RK4
    n_sub = 5; dt_sub = dt/n_sub; h2 = dt_sub/2;
    u_d = zeros(Nd,1); v_d = zeros(Nd,1); w_d = zeros(Nd,1);
    for k = 1:Nd-1
        pk=ds_pqr(k,1); qk=ds_pqr(k,2); rk=ds_pqr(k,3);
        gxk=gx_d(k); gyk=gy_d(k); gzk=gz_d(k);
        Tk_m = T_tot_d(k)/params.mass;
        us=u_d(k); vs=v_d(k); ws=w_d(k);
        for si=1:n_sub
            ud1=rk*vs-qk*ws+gxk+Xu_f*us; vd1=pk*ws-rk*us+gyk+Yv_f*vs; wd1=qk*us-pk*vs-Tk_m+gzk+Zw_f*ws+Bz_f;
            u2=us+h2*ud1; v2=vs+h2*vd1; w2=ws+h2*wd1;
            ud2=rk*v2-qk*w2+gxk+Xu_f*u2; vd2=pk*w2-rk*u2+gyk+Yv_f*v2; wd2=qk*u2-pk*v2-Tk_m+gzk+Zw_f*w2+Bz_f;
            u3=us+h2*ud2; v3=vs+h2*vd2; w3=ws+h2*wd2;
            ud3=rk*v3-qk*w3+gxk+Xu_f*u3; vd3=pk*w3-rk*u3+gyk+Yv_f*v3; wd3=qk*u3-pk*v3-Tk_m+gzk+Zw_f*w3+Bz_f;
            u4=us+dt_sub*ud3; v4=vs+dt_sub*vd3; w4=ws+dt_sub*wd3;
            ud4=rk*v4-qk*w4+gxk+Xu_f*u4; vd4=pk*w4-rk*u4+gyk+Yv_f*v4; wd4=qk*u4-pk*v4-Tk_m+gzk+Zw_f*w4+Bz_f;
            us=us+dt_sub/6*(ud1+2*ud2+2*ud3+ud4);
            vs=vs+dt_sub/6*(vd1+2*vd2+2*vd3+vd4);
            ws=ws+dt_sub/6*(wd1+2*wd2+2*wd3+wd4);
        end
        u_d(k+1)=us; v_d(k+1)=vs; w_d(k+1)=ws;
    end

    % Forca especifica modelo (u_dot, v_dot incluem gravidade -> compativeis com IMU)
    accX_mod = ds_pqr(:,3).*v_d - ds_pqr(:,2).*w_d + gx_d + Xu_f*u_d;
    accY_mod = ds_pqr(:,1).*w_d - ds_pqr(:,3).*u_d + gy_d + Yv_f*v_d;
    accZ_mod = -T_tot_d/params.mass + Zw_f*w_d + Bz_f;

    R2_ax = R2_func(ds_acc(:,1), accX_mod);
    R2_ay = R2_func(ds_acc(:,2), accY_mod);
    R2_az = R2_func(ds_acc(:,3), accZ_mod);

    fprintf('  %s [%d-%ds]:\n', ds_name, ds_range(1), ds_range(2));
    fprintf('    R2 AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', R2_ax, R2_ay, R2_az);

    % Salvar para plot
    if d == 2  % validacao
        accX_mod_vl = accX_mod; accY_mod_vl = accY_mod; accZ_mod_vl = accZ_mod;
        u_vl = u_d; v_vl = v_d; w_vl = w_d;
    end
end

%% ========================================================================
%  11. FIGURAS
% =========================================================================

% --- Fig 12: EEM derivadas (treino) ---
fig12 = figure('Position', [50 50 1100 700], 'Color', 'w', 'Visible', 'off');
subplot(3,1,1);
plot(1:N_eem, pd_all, 'b-', 1:N_eem, pd_pred, 'r--', 'LineWidth', 0.8);
ylabel('p_{dot} [rad/s^2]'); grid on;
legend('Medido', sprintf('EEM (R^2=%.3f)', R2_func(pd_all, pd_pred)));
title('Fase A - EEM: Derivadas preditas vs medidas (treino concatenado)');
subplot(3,1,2);
plot(1:N_eem, qd_all, 'b-', 1:N_eem, qd_pred, 'r--', 'LineWidth', 0.8);
ylabel('q_{dot} [rad/s^2]'); grid on;
legend('Medido', sprintf('EEM (R^2=%.3f)', R2_func(qd_all, qd_pred)));
subplot(3,1,3);
plot(1:N_eem, rd_all, 'b-', 1:N_eem, rd_pred, 'r--', 'LineWidth', 0.8);
ylabel('r_{dot} [rad/s^2]'); xlabel('Amostra'); grid on;
legend('Medido', sprintf('EEM (R^2=%.3f)', R2_func(rd_all, rd_pred)));
saveas(fig12, fullfile(img_dir, 'fig12_eem_derivadas.png'));

% --- Fig 13: Validacao p,q,r ---
if vl_ok
    fig13 = figure('Position', [80 50 1000 700], 'Color', 'w', 'Visible', 'off');
    subplot(3,1,1);
    plot(time_vl, pqr_vl(:,1), 'b-', time_vl, x_vl(:,1), 'r--', 'LineWidth', 1.2);
    ylabel('p [rad/s]'); grid on;
    legend('Medido', sprintf('Sim (R^2=%.3f)', R2_func(pqr_vl(:,1), x_vl(:,1))));
    title(sprintf('Validacao [%d-%ds] - Vel. Angulares (9 estados, ode45)', t_val(1), t_val(2)));
    subplot(3,1,2);
    plot(time_vl, pqr_vl(:,2), 'b-', time_vl, x_vl(:,2), 'r--', 'LineWidth', 1.2);
    ylabel('q [rad/s]'); grid on;
    legend('Medido', sprintf('Sim (R^2=%.3f)', R2_func(pqr_vl(:,2), x_vl(:,2))));
    subplot(3,1,3);
    plot(time_vl, pqr_vl(:,3), 'b-', time_vl, x_vl(:,3), 'r--', 'LineWidth', 1.2);
    ylabel('r [rad/s]'); xlabel('Tempo [s]'); grid on;
    legend('Medido', sprintf('Sim (R^2=%.3f)', R2_func(pqr_vl(:,3), x_vl(:,3))));
    saveas(fig13, fullfile(img_dir, 'fig13_val_pqr.png'));

    % --- Fig 14: Validacao atitude ---
    fig14 = figure('Position', [110 50 1000 700], 'Color', 'w', 'Visible', 'off');
    subplot(3,1,1);
    plot(time_vl, att_vl(:,1), 'b-', time_vl, rad2deg(x_vl(:,4)), 'r--', 'LineWidth', 1.2);
    ylabel('Roll [deg]'); grid on;
    legend('EKF', sprintf('Sim (R^2=%.3f)', R2_func(att_vl(:,1), rad2deg(x_vl(:,4)))));
    title(sprintf('Validacao [%d-%ds] - Atitude (9 estados, ode45)', t_val(1), t_val(2)));
    subplot(3,1,2);
    plot(time_vl, att_vl(:,2), 'b-', time_vl, rad2deg(x_vl(:,5)), 'r--', 'LineWidth', 1.2);
    ylabel('Pitch [deg]'); grid on;
    legend('EKF', sprintf('Sim (R^2=%.3f)', R2_func(att_vl(:,2), rad2deg(x_vl(:,5)))));
    subplot(3,1,3);
    plot(time_vl, att_vl(:,3), 'b-', time_vl, rad2deg(x_vl(:,6)), 'r--', 'LineWidth', 1.2);
    ylabel('Yaw [deg]'); xlabel('Tempo [s]'); grid on;
    legend('EKF', sprintf('Sim (R^2=%.3f)', R2_func(att_vl(:,3), rad2deg(x_vl(:,6)))));
    saveas(fig14, fullfile(img_dir, 'fig14_val_atitude.png'));
end

% --- Fig 15: Aceleracoes (validacao semi-acoplada) ---
fig15 = figure('Position', [50 50 1100 800], 'Color', 'w', 'Visible', 'off');
subplot(3,1,1);
plot(time_vl, acc_vl(:,1), 'b-', time_vl, accX_mod_vl, 'r--', 'LineWidth', 1.0);
ylabel('AccX [m/s^2]'); grid on;
legend('IMU', sprintf('Modelo (R^2=%.3f)', R2_func(acc_vl(:,1), accX_mod_vl)));
title(sprintf('Validacao [%d-%ds] - Aceleracoes (semi-acoplado)', t_val(1), t_val(2)));
subplot(3,1,2);
plot(time_vl, acc_vl(:,2), 'b-', time_vl, accY_mod_vl, 'r--', 'LineWidth', 1.0);
ylabel('AccY [m/s^2]'); grid on;
legend('IMU', sprintf('Modelo (R^2=%.3f)', R2_func(acc_vl(:,2), accY_mod_vl)));
subplot(3,1,3);
plot(time_vl, acc_vl(:,3), 'b-', time_vl, accZ_mod_vl, 'r--', 'LineWidth', 1.0);
ylabel('AccZ [m/s^2]'); xlabel('Tempo [s]'); grid on;
legend('IMU', sprintf('Modelo (R^2=%.3f)', R2_func(acc_vl(:,3), accZ_mod_vl)));
saveas(fig15, fullfile(img_dir, 'fig15_val_aceleracoes.png'));

% --- Fig 16: Parametros bar chart ---
fig16 = figure('Position', [50 50 1200 500], 'Color', 'w', 'Visible', 'off');
subplot(1,3,1);
bar_J = [P0(1:3)*1000, P_eem(1:3)*1000, P_final(1:3)*1000];
b = bar(bar_J);
b(1).FaceColor = [0.7 0.7 0.7]; b(2).FaceColor = [0.3 0.5 0.8]; b(3).FaceColor = [0.2 0.7 0.3];
set(gca, 'XTickLabel', {'Jx','Jy','Jz'});
ylabel('J [g.m^2]'); grid on; legend('CAD','EEM','Final');
title('Inercias');

subplot(1,3,2);
bar_kT = [P0(7:10), P_eem(7:10), P_final(7:10)];
b = bar(bar_kT);
b(1).FaceColor = [0.7 0.7 0.7]; b(2).FaceColor = [0.3 0.5 0.8]; b(3).FaceColor = [0.2 0.7 0.3];
set(gca, 'XTickLabel', {'C1','C2','C3','C4'});
ylabel('k_T'); grid on; legend('P0','EEM','Final');
title('Empuxo k_T'); yline(1, 'k--');

subplot(1,3,3);
bar_kQ = [P0(11:14), P_eem(11:14), P_final(11:14)];
b = bar(bar_kQ);
b(1).FaceColor = [0.7 0.7 0.7]; b(2).FaceColor = [0.3 0.5 0.8]; b(3).FaceColor = [0.2 0.7 0.3];
set(gca, 'XTickLabel', {'C1','C2','C3','C4'});
ylabel('k_Q'); grid on; legend('P0','EEM','Final');
title('Torque k_Q'); yline(1, 'k--');

saveas(fig16, fullfile(img_dir, 'fig16_parametros.png'));

% --- Fig 17: Velocidades u,v,w integradas (validacao) ---
fig17 = figure('Position', [80 50 1000 600], 'Color', 'w', 'Visible', 'off');
subplot(3,1,1);
plot(time_vl, u_vl, 'r-', 'LineWidth', 1.2);
ylabel('u [m/s]'); grid on; title(sprintf('Velocidades lineares integradas (validacao [%d-%ds])', t_val(1), t_val(2)));
subplot(3,1,2);
plot(time_vl, v_vl, 'g-', 'LineWidth', 1.2);
ylabel('v [m/s]'); grid on;
subplot(3,1,3);
plot(time_vl, w_vl, 'b-', 'LineWidth', 1.2);
ylabel('w [m/s]'); xlabel('Tempo [s]'); grid on;
saveas(fig17, fullfile(img_dir, 'fig17_velocidades_uvw.png'));

fprintf('\n=== Identificacao concluida! ===\n');
fprintf('  Figuras salvas em: %s\n', img_dir);

%% ========================================================================
%  FUNCOES LOCAIS
% =========================================================================

function e = oem_multi_seg_cost(P, segs, win_sec, params, dt, weights_pqr, weights_acc)
    e_all = [];
    for s = 1:length(segs)
        sg = segs{s};
        e_seg = oem_seg_cost(P, sg, win_sec, params, dt, weights_pqr, weights_acc);
        e_all = [e_all; e_seg]; %#ok<AGROW>
    end
    e = e_all;
end

function e = oem_seg_cost(P, sg, win_sec, params, dt, weights_pqr, weights_acc)
% OEM: inercias variaveis + k + D + B + translacional
    N = sg.N;
    x_m = params.x_m; y_m = params.y_m; d_m = params.d_m;
    m_kg = params.mass; g = params.g;

    % Inercias do vetor P
    J = [P(1), P(4), P(5); P(4), P(2), P(6); P(5), P(6), P(3)];
    J_inv = inv(J);

    k_T = P(7:10); k_Q = P(11:14);
    Dp = P(15); Dq = P(16); Dr = P(17);
    Bp = P(18); Bq = P(19); Br = P(20);
    Xu = P(21); Yv = P(22); Zw = P(23); Bz = P(24);

    Ti_mat = sg.T_ref .* k_T'; Qi_mat = sg.Q_ref .* k_Q';
    Mx_motor = sum(-y_m' .* Ti_mat, 2) + Bp;
    My_motor = sum(+x_m' .* Ti_mat, 2) + Bq;
    Mz_motor = sum( d_m' .* Qi_mat, 2) + Br;
    T_total = sum(Ti_mat, 2);

    win_len = round(win_sec / dt);
    win_starts = 1:win_len:N;
    if length(win_starts) > 1 && win_starts(end) == N, win_starts(end) = []; end
    win_ends = min([win_starts(2:end)-1, N], N);

    n_sub = 5; dt_sub = dt/n_sub; h2 = dt_sub/2; MAX_VAL = 50;

    rot_dot = @(omega_v, Mm) J_inv * (-cross(omega_v, J*omega_v) + ...
        Mm - [Dp*omega_v(1); Dq*omega_v(2); Dr*omega_v(3)]);

    p_sim = zeros(N,1); q_sim = zeros(N,1); r_sim = zeros(N,1);

    for w = 1:length(win_starts)
        i_s = win_starts(w); i_e = win_ends(w);
        ps = sg.pqr(i_s,1); qs = sg.pqr(i_s,2); rs = sg.pqr(i_s,3);
        p_sim(i_s) = ps; q_sim(i_s) = qs; r_sim(i_s) = rs;

        for k = i_s:i_e-1
            Mm_k = [Mx_motor(k); My_motor(k); Mz_motor(k)];
            for si = 1:n_sub
                omega = [ps; qs; rs];
                od1 = rot_dot(omega, Mm_k);
                od2 = rot_dot(omega + h2*od1, Mm_k);
                od3 = rot_dot(omega + h2*od2, Mm_k);
                od4 = rot_dot(omega + dt_sub*od3, Mm_k);
                omega = omega + dt_sub/6*(od1 + 2*od2 + 2*od3 + od4);
                ps = omega(1); qs = omega(2); rs = omega(3);
            end
            if ~isfinite(ps), ps=MAX_VAL; end
            if ~isfinite(qs), qs=MAX_VAL; end
            if ~isfinite(rs), rs=MAX_VAL; end
            ps=max(min(ps,MAX_VAL),-MAX_VAL);
            qs=max(min(qs,MAX_VAL),-MAX_VAL);
            rs=max(min(rs,MAX_VAL),-MAX_VAL);
            p_sim(k+1)=ps; q_sim(k+1)=qs; r_sim(k+1)=rs;
        end
    end

    sw_r = sqrt(weights_pqr(:));
    e_rot = [sw_r(1)*(sg.pqr(:,1)-p_sim); sw_r(2)*(sg.pqr(:,2)-q_sim); sw_r(3)*(sg.pqr(:,3)-r_sim)];

    % Translacional semi-acoplado
    phi_m = sg.att_rad(:,1); theta_m = sg.att_rad(:,2);
    gx = -g*sin(theta_m); gy = g*cos(theta_m).*sin(phi_m); gz = g*cos(theta_m).*cos(phi_m);

    u_int=zeros(N,1); v_int=zeros(N,1); w_int=zeros(N,1);
    for k = 1:N-1
        pk=sg.pqr(k,1); qk=sg.pqr(k,2); rk=sg.pqr(k,3);
        gxk=gx(k); gyk=gy(k); gzk=gz(k); Tk_m=T_total(k)/m_kg;
        us=u_int(k); vs=v_int(k); ws=w_int(k);
        for si=1:n_sub
            ud1=rk*vs-qk*ws+gxk+Xu*us; vd1=pk*ws-rk*us+gyk+Yv*vs; wd1=qk*us-pk*vs-Tk_m+gzk+Zw*ws+Bz;
            u2=us+h2*ud1; v2=vs+h2*vd1; w2=ws+h2*wd1;
            ud2=rk*v2-qk*w2+gxk+Xu*u2; vd2=pk*w2-rk*u2+gyk+Yv*v2; wd2=qk*u2-pk*v2-Tk_m+gzk+Zw*w2+Bz;
            u3=us+h2*ud2; v3=vs+h2*vd2; w3=ws+h2*wd2;
            ud3=rk*v3-qk*w3+gxk+Xu*u3; vd3=pk*w3-rk*u3+gyk+Yv*v3; wd3=qk*u3-pk*v3-Tk_m+gzk+Zw*w3+Bz;
            u4=us+dt_sub*ud3; v4=vs+dt_sub*vd3; w4=ws+dt_sub*wd3;
            ud4=rk*v4-qk*w4+gxk+Xu*u4; vd4=pk*w4-rk*u4+gyk+Yv*v4; wd4=qk*u4-pk*v4-Tk_m+gzk+Zw*w4+Bz;
            us=us+dt_sub/6*(ud1+2*ud2+2*ud3+ud4);
            vs=vs+dt_sub/6*(vd1+2*vd2+2*vd3+vd4);
            ws=ws+dt_sub/6*(wd1+2*wd2+2*wd3+wd4);
        end
        if ~isfinite(us), us=MAX_VAL; end; if ~isfinite(vs), vs=MAX_VAL; end; if ~isfinite(ws), ws=MAX_VAL; end
        us=max(min(us,MAX_VAL),-MAX_VAL); vs=max(min(vs,MAX_VAL),-MAX_VAL); ws=max(min(ws,MAX_VAL),-MAX_VAL);
        u_int(k+1)=us; v_int(k+1)=vs; w_int(k+1)=ws;
    end

    accX_m = sg.pqr(:,3).*v_int - sg.pqr(:,2).*w_int + gx + Xu*u_int;
    accY_m = sg.pqr(:,1).*w_int - sg.pqr(:,3).*u_int + gy + Yv*v_int;
    accZ_m = -T_total/m_kg + Zw*w_int + Bz;

    sw_a = sqrt(weights_acc(:));
    e_acc = [sw_a(1)*(sg.acc(:,1)-accX_m); sw_a(2)*(sg.acc(:,2)-accY_m); sw_a(3)*(sg.acc(:,3)-accZ_m)];

    e = [e_rot; e_acc];
end

function [R2_p, R2_q, R2_r] = evaluate_rotational(P, params, ~, pwm_vl, pqr_vl, dt)
    R2f = @(y_e,y_s) 1-sum((y_e-y_s).^2)/max(sum((y_e-mean(y_e)).^2),1e-12);

    N_vl = length(pqr_vl(:,1));
    x_m = params.x_m; y_m = params.y_m; d_m = params.d_m;

    J = [P(1),P(4),P(5); P(4),P(2),P(6); P(5),P(6),P(3)];
    J_inv = inv(J);

    k_T=P(7:10); k_Q=P(11:14);
    Dp=P(15); Dq=P(16); Dr=P(17);
    Bp=P(18); Bq=P(19); Br=P(20);

    T_ref_vl=zeros(N_vl,4); Q_ref_vl=zeros(N_vl,4);
    for j=1:4
        T_ref_vl(:,j)=params.T_ref(pwm_vl(:,j));
        Q_ref_vl(:,j)=params.Q_ref(pwm_vl(:,j));
    end
    Ti=T_ref_vl.*k_T'; Qi=Q_ref_vl.*k_Q';
    Mx_motor=sum(-y_m'.*Ti,2)+Bp; My_motor=sum(+x_m'.*Ti,2)+Bq; Mz_motor=sum(d_m'.*Qi,2)+Br;

    rot_dot=@(o,Mm) J_inv*(-cross(o,J*o)+Mm-[Dp*o(1);Dq*o(2);Dr*o(3)]);

    n_sub=5; dt_sub=dt/n_sub; h2=dt_sub/2;
    ps=pqr_vl(1,1); qs=pqr_vl(1,2); rs=pqr_vl(1,3);
    p_sim=zeros(N_vl,1); q_sim=zeros(N_vl,1); r_sim=zeros(N_vl,1);
    p_sim(1)=ps; q_sim(1)=qs; r_sim(1)=rs;

    for k=1:N_vl-1
        Mm_k=[Mx_motor(k);My_motor(k);Mz_motor(k)];
        for si=1:n_sub
            omega=[ps;qs;rs];
            od1=rot_dot(omega,Mm_k); od2=rot_dot(omega+h2*od1,Mm_k);
            od3=rot_dot(omega+h2*od2,Mm_k); od4=rot_dot(omega+dt_sub*od3,Mm_k);
            omega=omega+dt_sub/6*(od1+2*od2+2*od3+od4);
            ps=omega(1); qs=omega(2); rs=omega(3);
        end
        if ~isfinite(ps), ps=0; end; if ~isfinite(qs), qs=0; end; if ~isfinite(rs), rs=0; end
        p_sim(k+1)=ps; q_sim(k+1)=qs; r_sim(k+1)=rs;
    end

    R2_p=R2f(pqr_vl(:,1),p_sim); R2_q=R2f(pqr_vl(:,2),q_sim); R2_r=R2f(pqr_vl(:,3),r_sim);
end
