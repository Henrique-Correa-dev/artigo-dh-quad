% simulate.m - Simulador standalone do modelo VTOL (equivalente ao Simulink)
%
% Carrega dados de voo, integra os 9 estados via vtol_dynamics.m (ode45),
% e gera graficos comparativos (pqr, att, acc) na pasta output/.
%
% Uso:
%   >> simulate              % usa trecho padrao 147-157s e P0
%   >> simulate              % edite t_range e P abaixo conforme necessario
%
% Equivale ao modelo Simulink quad_model_v3: todos os subsistemas
% (rotacional, cinematica de Euler, translacional) integrados simultaneamente
% pelo mesmo solver, sem interpolacao entre etapas.

%% ========================================================================
%  1. CONFIGURACAO
%  ========================================================================
% Trecho de tempo a simular (segundos do log)
t_range = [147, 157];

% ---------- Selecao de parametros (comente/descomente) ----------
% Opcao 1: Chute inicial (P0)
P = [0.063244; 0.250554; 0.116192; 0.001571;   ... % Jx, Jy, Jz, Jxz
     0.55; 0.45; 1.0; 0.75;                    ... % k_T1..k_T4
     0.55; 0.45; 1.0; 0.75;                    ... % k_Q1..k_Q4
     10; 5; 0.5;                                ... % Dp, Dq, Dr
     0.7; 1.4; 0.3;                            ... % Bp, Bq, Br
     0.0; 0.0;                                  ... % dx_cg, dy_cg
     -4.0; -4.0; -0.1; -0.5];                     % Xu_m, Yv_m, Zw_m, Bz
param_source = 'P0 (chute inicial)';

%Opcao 2: Parametros identificados (descomente as 3 linhas abaixo)
%id = load(fullfile(fileparts(mfilename('fullpath')), 'P_identified.mat'));
%P = id.P_final;
%param_source = 'P_final (identificado)';

% Constantes fisicas
constants.m = 1.6011;   % massa (kg)
constants.g = 9.81;     % gravidade (m/s^2)

% Solver
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Pasta de saida
output_dir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ========================================================================
%  2. CARREGAMENTO E INTERPOLACAO DOS DADOS
%  ========================================================================
fprintf('Carregando dados...\n');
load("log_data.mat")

ATT.TimeS  = double(ATT.TimeUS) / 1e6;
IMU.TimeS  = double(IMU.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;

idx = IMU.I == 0;
gyrX_raw = IMU.GyrX(idx); gyrY_raw = IMU.GyrY(idx); gyrZ_raw = IMU.GyrZ(idx);
accX_raw = IMU.AccX(idx); accY_raw = IMU.AccY(idx); accZ_raw = IMU.AccZ(idx);
time_IMU = IMU.TimeS(idx);

time_ATT  = ATT.TimeS;
time_RCOU = RCOU.TimeS;

pwm1_raw = double(RCOU.C1); pwm2_raw = double(RCOU.C2);
pwm3_raw = double(RCOU.C3); pwm4_raw = double(RCOU.C4);

t_start = max([min(time_IMU), min(time_ATT), min(time_RCOU)]);
t_end   = min([max(time_IMU), max(time_ATT), max(time_RCOU)]);
dt = 0.1;
t_common = t_start:dt:t_end;

gyrX_interp  = interp1(time_IMU, gyrX_raw, t_common, 'linear');
gyrY_interp  = interp1(time_IMU, gyrY_raw, t_common, 'linear');
gyrZ_interp  = interp1(time_IMU, gyrZ_raw, t_common, 'linear');
accX_interp  = interp1(time_IMU, accX_raw, t_common, 'linear');
accY_interp  = interp1(time_IMU, accY_raw, t_common, 'linear');
accZ_interp  = interp1(time_IMU, accZ_raw, t_common, 'linear');
roll_interp  = interp1(time_ATT, ATT.Roll, t_common, 'linear');
pitch_interp = interp1(time_ATT, ATT.Pitch, t_common, 'linear');
yaw_interp   = interp1(time_ATT, ATT.Yaw, t_common, 'linear');
pwm1_interp  = interp1(time_RCOU, pwm1_raw, t_common, 'linear');
pwm2_interp  = interp1(time_RCOU, pwm2_raw, t_common, 'linear');
pwm3_interp  = interp1(time_RCOU, pwm3_raw, t_common, 'linear');
pwm4_interp  = interp1(time_RCOU, pwm4_raw, t_common, 'linear');

fprintf('Dados interpolados a 10 Hz. Base de tempo: %.1f a %.1f s (%d pontos)\n', ...
    t_common(1), t_common(end), length(t_common));

%% ========================================================================
%  3. MODELOS DE REFERENCIA DOS MOTORES
%  ========================================================================
pwm_values_exp   = [1000; 1200; 1400; 1600; 1800; 2000];
thrust_grams_exp = [0; 143; 328; 532; 784; 843];
torque_Nm_exp    = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];
poly_degree = 3;

func_T_ref = create_thrust_model(pwm_values_exp, thrust_grams_exp, poly_degree);
func_Q_ref = create_torque_model(pwm_values_exp, torque_Nm_exp, poly_degree);

%% ========================================================================
%  4. EXTRACAO DO TRECHO DE SIMULACAO
%  ========================================================================
idx_sim = (t_common >= t_range(1)) & (t_common <= t_range(2));
time = t_common(idx_sim)';
N = sum(idx_sim);

pwm = [pwm1_interp(idx_sim)', pwm2_interp(idx_sim)', ...
       pwm3_interp(idx_sim)', pwm4_interp(idx_sim)'];
pqr = [gyrX_interp(idx_sim)', gyrY_interp(idx_sim)', gyrZ_interp(idx_sim)'];
acc = [accX_interp(idx_sim)', accY_interp(idx_sim)', accZ_interp(idx_sim)'];
att = [roll_interp(idx_sim)',  pitch_interp(idx_sim)', yaw_interp(idx_sim)'];

fprintf('\nTrecho: %d a %d s  (%d pontos, dt=%.2f s)\n', t_range(1), t_range(2), N, dt);
fprintf('Parametros: %s\n', param_source);

%% ========================================================================
%  5. INTEGRACAO 9 ESTADOS (vtol_dynamics via ode45)
%  ========================================================================
fprintf('Integrando 9 estados via ode45...\n');

% Condicao inicial: [p,q,r] e [phi,theta,psi] medidos, [u,v,w] = 0
att_rad0 = deg2rad(att(1,:));
y0 = [pqr(1,:)'; att_rad0(:); 0; 0; 0];

fprintf('  y0 = [p=%.4f, q=%.4f, r=%.4f, phi=%.4f, theta=%.4f, psi=%.4f, u=0, v=0, w=0]\n', ...
    y0(1), y0(2), y0(3), y0(4), y0(5), y0(6));

ode_func = @(t,y) vtol_dynamics(t, y, P, time, pwm, func_T_ref, func_Q_ref, constants);
tic;
[t_s, y_s] = ode45(ode_func, time, y0, ode_opts);
tempo_sim = toc;
fprintf('Integracao concluida em %.2f s  (%d passos adaptativos)\n', tempo_sim, length(t_s));

% Interpolar para a grade 10 Hz
y_out = interp1(t_s, y_s, time, 'linear', 'extrap');

% Extrair estados
p_sim     = y_out(:,1);
q_sim     = y_out(:,2);
r_sim     = y_out(:,3);
phi_sim   = rad2deg(y_out(:,4));
theta_sim = rad2deg(y_out(:,5));
psi_sim   = rad2deg(y_out(:,6));
u_sim     = y_out(:,7);
v_sim     = y_out(:,8);
w_sim     = y_out(:,9);

%% ========================================================================
%  6. CALCULO DA FORCA ESPECIFICA (AccX, AccY, AccZ)
%  ========================================================================
% Utiliza subfunções centralizadas do vtol_dynamics.m.
% Para testar modelos diferentes, edite APENAS vtol_dynamics.m.

dyn_h = vtol_dynamics('get_handles');
trans_dot_fn  = dyn_h.trans_dot;
spec_force_fn = dyn_h.specific_force;

g = constants.g;
m = constants.m;
Xu_m = P(21); Yv_m = P(22); Zw_m = P(23); Bz = P(24);
k_T = P(5:8);

% Empuxo total vetorizado
T_vec = zeros(N,1);
for k = 1:N
    for j = 1:4
        T_vec(k) = T_vec(k) + k_T(j) * func_T_ref(pwm(k,j));
    end
end

% Projecao da gravidade no corpo
gx = -g * sin(y_out(:,5));
gy =  g * sin(y_out(:,4)) .* cos(y_out(:,5));
gz =  g * cos(y_out(:,4)) .* cos(y_out(:,5));

% Derivadas translacionais (u_dot, v_dot, w_dot) — via subfunção centralizada
[udot_sim, vdot_sim, wdot_sim] = trans_dot_fn( ...
    y_out(:,1), y_out(:,2), y_out(:,3), ...
    y_out(:,7), y_out(:,8), y_out(:,9), ...
    gx, gy, gz, T_vec/m, Xu_m, Yv_m, Zw_m, Bz);

% Saída do modelo (derivada completa, gz subtraído em z) — via subfunção centralizada
[accX_sim, accY_sim, accZ_sim] = spec_force_fn( ...
    y_out(:,1), y_out(:,2), y_out(:,3), ...
    y_out(:,7), y_out(:,8), y_out(:,9), ...
    gx, gy, gz, T_vec/m, Xu_m, Yv_m, Zw_m, Bz);

%% ========================================================================
%  6b. DIAGNOSTICO: TRANSLACIONAL COM ATITUDE MEDIDA (isola erros de att)
%  ========================================================================
% Integra APENAS u,v,w via RK4, usando p,q,r e phi,theta MEDIDOS.
% Serve para testar o modelo translacional sem cascata de erros da atitude.

fprintf('\n--- Diagnostico: translacional com atitude medida ---\n');

att_rad = deg2rad(att);       % [phi, theta, psi] em rad
phi_m   = att_rad(:,1);
theta_m = att_rad(:,2);
dt = time(2) - time(1);

% Gravidade projetada a partir de atitude MEDIDA
gx_m = -g * sin(theta_m);
gy_m =  g * sin(phi_m) .* cos(theta_m);
gz_m =  g * cos(phi_m) .* cos(theta_m);

% RK4 para u,v,w usando p,q,r e atitude medidos
n_sub = 5;          % sub-passos por amostra
dt_sub = dt / n_sub;
h2 = dt_sub / 2;

u_diag = zeros(N,1); v_diag = zeros(N,1); w_diag = zeros(N,1);
for k = 1:N-1
    pk = pqr(k,1); qk = pqr(k,2); rk = pqr(k,3);
    gxk = gx_m(k); gyk = gy_m(k); gzk = gz_m(k);
    Tk_m = T_vec(k)/m;

    us = u_diag(k); vs = v_diag(k); ws = w_diag(k);
    for si = 1:n_sub
        [ud1,vd1,wd1] = trans_dot_fn(pk,qk,rk, us,vs,ws, gxk,gyk,gzk, Tk_m, Xu_m,Yv_m,Zw_m,Bz);
        u2=us+h2*ud1; v2=vs+h2*vd1; w2=ws+h2*wd1;
        [ud2,vd2,wd2] = trans_dot_fn(pk,qk,rk, u2,v2,w2, gxk,gyk,gzk, Tk_m, Xu_m,Yv_m,Zw_m,Bz);
        u3=us+h2*ud2; v3=vs+h2*vd2; w3=ws+h2*wd2;
        [ud3,vd3,wd3] = trans_dot_fn(pk,qk,rk, u3,v3,w3, gxk,gyk,gzk, Tk_m, Xu_m,Yv_m,Zw_m,Bz);
        u4=us+dt_sub*ud3; v4=vs+dt_sub*vd3; w4=ws+dt_sub*wd3;
        [ud4,vd4,wd4] = trans_dot_fn(pk,qk,rk, u4,v4,w4, gxk,gyk,gzk, Tk_m, Xu_m,Yv_m,Zw_m,Bz);
        us = us + dt_sub/6*(ud1 + 2*ud2 + 2*ud3 + ud4);
        vs = vs + dt_sub/6*(vd1 + 2*vd2 + 2*vd3 + vd4);
        ws = ws + dt_sub/6*(wd1 + 2*wd2 + 2*wd3 + wd4);
    end
    u_diag(k+1) = us; v_diag(k+1) = vs; w_diag(k+1) = ws;
end

% Saída do modelo com u,v,w do diagnostico (atitude medida)
[accX_diag, accY_diag, accZ_diag] = spec_force_fn( ...
    pqr(:,1), pqr(:,2), pqr(:,3), ...
    u_diag, v_diag, w_diag, ...
    gx_m, gy_m, gz_m, T_vec/m, Xu_m, Yv_m, Zw_m, Bz);

%% ========================================================================
%  7. METRICAS R^2
%  ========================================================================
R2 = @(y_e, y_s) 1 - sum((y_e - y_s).^2) / max(sum((y_e - mean(y_e)).^2), 1e-12);

R2_p = R2(pqr(:,1), p_sim);
R2_q = R2(pqr(:,2), q_sim);
R2_r = R2(pqr(:,3), r_sim);

R2_phi   = R2(att(:,1), phi_sim);
R2_theta = R2(att(:,2), theta_sim);
R2_psi   = R2(att(:,3), psi_sim);

R2_ax = R2(acc(:,1), accX_sim);
R2_ay = R2(acc(:,2), accY_sim);
R2_az = R2(acc(:,3), accZ_sim);

% Diagnostico (atitude medida)
R2_ax_diag = R2(acc(:,1), accX_diag);
R2_ay_diag = R2(acc(:,2), accY_diag);
R2_az_diag = R2(acc(:,3), accZ_diag);

fprintf('\n==========================================================\n');
fprintf('  RESULTADOS — Trecho [%d-%ds]\n', t_range(1), t_range(2));
fprintf('==========================================================\n');
fprintf('  Vel. Angulares:\n');
fprintf('    R2 p = %.4f\n', R2_p);
fprintf('    R2 q = %.4f\n', R2_q);
fprintf('    R2 r = %.4f\n', R2_r);
fprintf('  Atitude:\n');
fprintf('    R2 phi   = %.4f\n', R2_phi);
fprintf('    R2 theta = %.4f\n', R2_theta);
fprintf('    R2 psi   = %.4f\n', R2_psi);
fprintf('  Aceleracoes (9-estados, atitude integrada):\n');
fprintf('    R2 AccX = %.4f\n', R2_ax);
fprintf('    R2 AccY = %.4f\n', R2_ay);
fprintf('    R2 AccZ = %.4f\n', R2_az);
fprintf('  Aceleracoes (DIAGNOSTICO, atitude medida):\n');
fprintf('    R2 AccX = %.4f\n', R2_ax_diag);
fprintf('    R2 AccY = %.4f\n', R2_ay_diag);
fprintf('    R2 AccZ = %.4f\n', R2_az_diag);
fprintf('==========================================================\n');

%% ========================================================================
%  8. GRAFICOS
%  ========================================================================
lbl = sprintf('[%d-%ds]', t_range(1), t_range(2));

% --- Figura 1: Velocidades angulares (p, q, r) ---
fig1 = figure('Name', 'pqr', 'Position', [80 50 900 500], 'Visible', 'off');

ax1=subplot(3,1,1);
plot(time, pqr(:,1), 'b-', time, p_sim, 'r--', 'LineWidth', 1.3);
legend('Experimental', 'Simulado', 'Location', 'best'); ylabel('p (rad/s)');
title(sprintf('%s — p  (R^2 = %.3f)', lbl, R2_p)); grid on;
yd = [pqr(:,1); p_sim]; yr = max(yd)-min(yd); if yr<1e-6, yr=1; end
ylim([min(yd)-0.5*yr, max(yd)+0.5*yr]);

ax2=subplot(3,1,2);
plot(time, pqr(:,2), 'b-', time, q_sim, 'r--', 'LineWidth', 1.3);
legend('Experimental', 'Simulado', 'Location', 'best'); ylabel('q (rad/s)');
title(sprintf('%s — q  (R^2 = %.3f)', lbl, R2_q)); grid on;
yd = [pqr(:,2); q_sim]; yr = max(yd)-min(yd); if yr<1e-6, yr=1; end
ylim([min(yd)-0.5*yr, max(yd)+0.5*yr]);

ax3=subplot(3,1,3);
plot(time, pqr(:,3), 'b-', time, r_sim, 'r--', 'LineWidth', 1.3);
legend('Experimental', 'Simulado', 'Location', 'best'); ylabel('r (rad/s)'); xlabel('Tempo (s)');
title(sprintf('%s — r  (R^2 = %.3f)', lbl, R2_r)); grid on;
yd = [pqr(:,3); r_sim]; yr = max(yd)-min(yd); if yr<1e-6, yr=1; end
ylim([min(yd)-0.5*yr, max(yd)+0.5*yr]);

sgtitle(sprintf('Velocidades Angulares — %s', lbl));
for a = [ax1 ax2 ax3], set(a, 'LooseInset', max(get(a,'TightInset'), 0.02)); end
saveas(fig1, fullfile(output_dir, 'pqr.png'));
fprintf('Salvo: %s\n', fullfile(output_dir, 'pqr.png'));

% --- Figura 2: Atitude (phi, theta, psi) ---
fig2 = figure('Name', 'att', 'Position', [120 90 900 700], 'Visible', 'off');

subplot(3,1,1);
plot(time, att(:,1), 'b-', time, phi_sim, 'r--', 'LineWidth', 1.3);
legend('EKF', 'Simulado'); ylabel('\phi (graus)');
title(sprintf('%s — \\phi  (R^2 = %.3f)', lbl, R2_phi)); grid on;

subplot(3,1,2);
plot(time, att(:,2), 'b-', time, theta_sim, 'r--', 'LineWidth', 1.3);
legend('EKF', 'Simulado'); ylabel('\theta (graus)');
title(sprintf('%s — \\theta  (R^2 = %.3f)', lbl, R2_theta)); grid on;

subplot(3,1,3);
plot(time, att(:,3), 'b-', time, psi_sim, 'r--', 'LineWidth', 1.3);
legend('EKF', 'Simulado'); ylabel('\psi (graus)'); xlabel('Tempo (s)');
title(sprintf('%s — \\psi  (R^2 = %.3f)', lbl, R2_psi)); grid on;

sgtitle(sprintf('Atitude — %s', lbl));
saveas(fig2, fullfile(output_dir, 'att.png'));
fprintf('Salvo: %s\n', fullfile(output_dir, 'att.png'));

% --- Figura 3: Aceleracoes / Forca especifica (AccX, AccY, AccZ) ---
fig3 = figure('Name', 'acc', 'Position', [160 130 900 500], 'Visible', 'off');

ax1=subplot(3,1,1);
plot(time, acc(:,1), 'b-', time, accX_sim, 'r--', 'LineWidth', 1.3);
legend('IMU', 'Simulado', 'Location', 'best'); ylabel('AccX (m/s^2)');
title(sprintf('%s — AccX  (R^2 = %.3f)', lbl, R2_ax)); grid on;
yd = acc(:,1); yr = max(yd)-min(yd); if yr<1e-6, yr=1; end
ylim([min(yd)-0.5*yr, max(yd)+0.5*yr]);

ax2=subplot(3,1,2);
plot(time, acc(:,2), 'b-', time, accY_sim, 'r--', 'LineWidth', 1.3);
legend('IMU', 'Simulado', 'Location', 'best'); ylabel('AccY (m/s^2)');
title(sprintf('%s — AccY  (R^2 = %.3f)', lbl, R2_ay)); grid on;
yd = acc(:,2); yr = max(yd)-min(yd); if yr<1e-6, yr=1; end
ylim([min(yd)-0.5*yr, max(yd)+0.5*yr]);

ax3=subplot(3,1,3);
plot(time, acc(:,3), 'b-', time, accZ_sim, 'r--', 'LineWidth', 1.3);
legend('IMU', 'Simulado', 'Location', 'best'); ylabel('AccZ (m/s^2)'); xlabel('Tempo (s)');
title(sprintf('%s — AccZ  (R^2 = %.3f)', lbl, R2_az)); grid on;
yd = acc(:,3); yr = max(yd)-min(yd); if yr<1e-6, yr=1; end
ylim([min(yd)-0.5*yr, max(yd)+0.5*yr]);

sgtitle(sprintf('Forca Especifica (Aceleracao IMU) — %s', lbl));
for a = [ax1 ax2 ax3], set(a, 'LooseInset', max(get(a,'TightInset'), 0.02)); end
saveas(fig3, fullfile(output_dir, 'acc.png'));
fprintf('Salvo: %s\n', fullfile(output_dir, 'acc.png'));

% --- Figura 4: Estados translacionais (u, v, w) ---
fig4 = figure('Name', 'uvw', 'Position', [200 170 900 700], 'Visible', 'off');

subplot(3,1,1);
plot(time, u_sim, 'r-', 'LineWidth', 1.3);
ylabel('u (m/s)'); title(sprintf('%s — u (body-x vel)', lbl)); grid on;

subplot(3,1,2);
plot(time, v_sim, 'r-', 'LineWidth', 1.3);
ylabel('v (m/s)'); title(sprintf('%s — v (body-y vel)', lbl)); grid on;

subplot(3,1,3);
plot(time, w_sim, 'r-', 'LineWidth', 1.3);
ylabel('w (m/s)'); xlabel('Tempo (s)');
title(sprintf('%s — w (body-z vel)', lbl)); grid on;

sgtitle(sprintf('Velocidades Translacionais — %s', lbl));
saveas(fig4, fullfile(output_dir, 'uvw.png'));
fprintf('Salvo: %s\n', fullfile(output_dir, 'uvw.png'));

% --- Figura 5: Comparacao u_dot vs f_especifica vs IMU ---
fig5 = figure('Name', 'acc_debug', 'Position', [240 210 900 700], 'Visible', 'off');

subplot(3,1,1);
plot(time, acc(:,1), 'b-', time, accX_sim, 'r--', time, udot_sim, 'g:', 'LineWidth', 1.3);
legend('IMU (f_{espec})', 'Xu_m \cdot u', 'u_{dot} (deriv)');
ylabel('m/s^2'); title(sprintf('%s — AccX: IMU vs f_{espec} vs u_{dot}', lbl)); grid on;

subplot(3,1,2);
plot(time, acc(:,2), 'b-', time, accY_sim, 'r--', time, vdot_sim, 'g:', 'LineWidth', 1.3);
legend('IMU (f_{espec})', 'Yv_m \cdot v', 'v_{dot} (deriv)');
ylabel('m/s^2'); title(sprintf('%s — AccY: IMU vs f_{espec} vs v_{dot}', lbl)); grid on;

subplot(3,1,3);
plot(time, acc(:,3), 'b-', time, accZ_sim, 'r--', time, wdot_sim, 'g:', 'LineWidth', 1.3);
legend('IMU (f_{espec})', '-T/m+Zw*w+Bz', 'w_{dot} (deriv)');
ylabel('m/s^2'); xlabel('Tempo (s)');
title(sprintf('%s — AccZ: IMU vs f_{espec} vs w_{dot}', lbl)); grid on;

sgtitle(sprintf('Debug: Forca Especifica vs Derivada de Velocidade — %s', lbl));
saveas(fig5, fullfile(output_dir, 'acc_debug.png'));
fprintf('Salvo: %s\n', fullfile(output_dir, 'acc_debug.png'));

% --- Figura 6: Diagnostico — Acc com atitude medida vs integrada ---
fig6 = figure('Name', 'acc_diag', 'Position', [280 250 900 700], 'Visible', 'off');

subplot(3,1,1);
plot(time, acc(:,1), 'b-', time, accX_sim, 'r--', time, accX_diag, 'g-.', 'LineWidth', 1.3);
legend('IMU', '9-estados', 'Diag (att medida)'); ylabel('AccX (m/s^2)');
title(sprintf('%s — AccX  9-est R^2=%.3f | Diag R^2=%.3f', lbl, R2_ax, R2_ax_diag)); grid on;

subplot(3,1,2);
plot(time, acc(:,2), 'b-', time, accY_sim, 'r--', time, accY_diag, 'g-.', 'LineWidth', 1.3);
legend('IMU', '9-estados', 'Diag (att medida)'); ylabel('AccY (m/s^2)');
title(sprintf('%s — AccY  9-est R^2=%.3f | Diag R^2=%.3f', lbl, R2_ay, R2_ay_diag)); grid on;

subplot(3,1,3);
plot(time, acc(:,3), 'b-', time, accZ_sim, 'r--', time, accZ_diag, 'g-.', 'LineWidth', 1.3);
legend('IMU', '9-estados', 'Diag (att medida)'); ylabel('AccZ (m/s^2)'); xlabel('Tempo (s)');
title(sprintf('%s — AccZ  9-est R^2=%.3f | Diag R^2=%.3f', lbl, R2_az, R2_az_diag)); grid on;

sgtitle(sprintf('Diagnostico Translacional: Atitude Medida vs Integrada — %s', lbl));
saveas(fig6, fullfile(output_dir, 'acc_diag.png'));
fprintf('Salvo: %s\n', fullfile(output_dir, 'acc_diag.png'));

fprintf('\nSimulacao concluida. Graficos salvos em: %s\n', output_dir);
