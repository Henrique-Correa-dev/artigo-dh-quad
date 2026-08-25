%SETUP_QUAD_LINEAR  Prepara o workspace pra simular quad_model_linear.slx.
%
% Modelo linear de HOVER, entrada em FORÇAS:  δẋ = A·δx + B·δv  (c=0 no equilíbrio),
% v = [T, Mx, My, Mz].
%
% NÃO recalcula a linearização — CARREGA de outputs/linear_model.mat (gerado
% por linearize.m), via load_linear_model (que avisa se estiver desatualizado).
% Assim o .slx usa EXATAMENTE o mesmo A,B do resto do projeto (fonte única).
% >>> Se você re-identificou (mudou P), rode linearize.m ANTES. <<<
%
% Pra validar contra voo, o modelo é alimentado por [T,M](t) DERIVADO do PWM
% logado (a ALOCAÇÃO, via vtol_dynamics forces handle + tabelas fT,fQ). Pra
% controle, basta desconectar os From Workspace e ligar o controlador no [T,M].
% VALE SÓ PERTO DO HOVER: valide num trecho calmo; em manobras grandes diverge.
%
% Popula no base workspace (lidos pelo bloco State-Space):
%   A_lin (9x9) B_lin (9x4) x0_lin (9x1) u0_lin (4x1=[mg;0;0;0]) dx0_lin (9x1)
%   T_in, Mx_in, My_in, Mz_in (timeseries)  |  t_sim  |  ref (struct p/ comparação)
%
% Uso:
%   >> setup_quad_linear
%   >> simOut = sim('quad_model_linear', 'StopTime', num2str(t_sim));

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ===================================================================
%  CONFIG
%  ===================================================================
LOG_FILE = 'logs_concat.mat';
t_window = [67, 87];            % janela de PWM que alimenta o modelo

%% ===================================================================
%  1) Carregar log e recortar janela (dt=0.1)
%  ===================================================================
L = load_log_data(fullfile(paths.data, LOG_FILE));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
t_common = (t_lo:0.1:t_hi)';
if t_window(1) < t_lo || t_window(2) > t_hi
    error('Janela [%g, %g] fora do log [%.1f, %.1f].', t_window, t_lo, t_hi);
end
idx = (t_common >= t_window(1)) & (t_common <= t_window(2));
time_w   = t_common(idx);
time_rel = time_w - time_w(1);
t_sim    = time_rel(end);

rs = @(tt, yy) interp1(tt, yy, t_common, 'linear');
pwm = [rs(L.time_RCOU,L.pwm1_raw), rs(L.time_RCOU,L.pwm2_raw), ...
       rs(L.time_RCOU,L.pwm3_raw), rs(L.time_RCOU,L.pwm4_raw)];
pqr = [rs(L.time_IMU,L.gyrX_raw), rs(L.time_IMU,L.gyrY_raw), rs(L.time_IMU,L.gyrZ_raw)];
acc = [rs(L.time_IMU,L.accX_raw), rs(L.time_IMU,L.accY_raw), rs(L.time_IMU,L.accZ_raw)];
att = [rs(L.time_ATT,L.roll_deg), rs(L.time_ATT,L.pitch_deg), rs(L.time_ATT,L.yaw_deg)];
pwm = pwm(idx,:);  pqr = pqr(idx,:);  acc = acc(idx,:);  att = att(idx,:);

% IC dos estados = medido em t=0 (u,v,w não medidos → 0)
x_ic = [pqr(1,:)'; deg2rad(att(1,:))'; 0; 0; 0];

%% ===================================================================
%  2) Carregar o modelo linear de hover (fonte única: linear_model.mat)
%  ===================================================================
lm = load_linear_model(paths);     % avisa se estiver stale vs P_identified.mat
A_lin  = lm.A;    B_lin  = lm.B;
x0_lin = lm.x0;   u0_lin = lm.u0;   % u0 = [m·g; 0; 0; 0] (trim de força)
dx0_lin = x_ic - x0_lin;            % δx(0) = IC medida - ponto de operação (hover)

%% ===================================================================
%  2b) ALOCAÇÃO: PWM logado → [T, Mx, My, Mz](t) (entradas do modelo)
%      Usa o handle forces (tabelas fT,fQ + k) — a mesma alocação do modelo NL.
%  ===================================================================
[func_T_ref, func_Q_ref] = motor_models();
dyn = vtol_dynamics('get_handles');
[Tt, Mxt, Myt, Mzt] = dyn.forces(pwm, lm.P, func_T_ref, func_Q_ref);   % Nx1 cada

T_in  = timeseries(Tt,  time_rel, 'Name','T_in');
Mx_in = timeseries(Mxt, time_rel, 'Name','Mx_in');
My_in = timeseries(Myt, time_rel, 'Name','My_in');
Mz_in = timeseries(Mzt, time_rel, 'Name','Mz_in');

fprintf('Modelo linear carregado de linear_model.mat\n');
fprintf('  v0 (trim de força) = [T=%.2f N, Mx=0, My=0, Mz=0]\n', u0_lin(1));
fprintf('  [T,M] da janela: T méd=%.2f N | Mx,My,Mz std=[%.3f %.3f %.3f] N·m\n', ...
        mean(Tt), std(Mxt), std(Myt), std(Mzt));

%% ===================================================================
%  3) ref auxiliar (pra comparação posterior em script)
%  ===================================================================
ref = struct('time', time_rel, 'time_abs', time_w, 'pqr', pqr, 'acc', acc, ...
             'att_deg', att, 'pwm', pwm, 'P_J', lm.P);

%% ===================================================================
%  Resumo
%  ===================================================================
fprintf('\nsetup_quad_linear pronto (hover):\n');
fprintf('  janela de entrada [%g, %g] s  (sim t=0 até %.2f s, %d amostras)\n', ...
        t_window(1), t_window(2), t_sim, sum(idx));
fprintf('  autovalores de A_lin: %s\n', mat2str(round(eig(A_lin),3)'));
fprintf('\nPra simular:\n');
fprintf('  simOut = sim(''quad_model_linear'', ''StopTime'', num2str(t_sim));\n');
