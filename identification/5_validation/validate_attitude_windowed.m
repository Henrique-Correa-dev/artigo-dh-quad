% validate_attitude_windowed.m — Validação da atitude por JANELA CURTA
%
% Prova que o modelo de atitude está CORRETO em horizonte curto, e que a
% deriva open-loop é só acúmulo de integração (não defeito de modelo).
%
% A atitude é φ,θ,ψ = ∫(cinemática de Euler dos pqr). Como é integrador puro
% (polo em 0), qualquer resíduo de pqr acumula → deriva ilimitada em janela
% longa. Mas em janela CURTA, ∫(resíduo) é pequeno → atitude bate.
%
% Compara 3 curvas por ângulo:
%   1. Medido (ATT do log)                      — referência
%   2. Modelo, JANELA CURTA (reset a cada N s)  — bate bem (kinematics OK)
%   3. Modelo, OPEN-LOOP (sem reset)            — deriva (artefato de integração)
%
% Tanto (2) quanto (3) integram o MESMO pqr do modelo (gerado do PWM).
% A ÚNICA diferença é o reset periódico. Se (2) bate e (3) deriva, está
% provado: a física está certa, a deriva é integração acumulada.

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ===================== CONFIG =====================
LOG_FILE = 'logs_concat.mat';
t_window = [75, 110];      % [s] janela de validação
N_RESET  = 3;              % [s] intervalo de reset da atitude (janela curta)
P_SOURCE = 'p_final';      % 'p0' | 'p_final' | 'manual'

%% ===================== CARREGAR =====================
proj = parameters();
switch lower(P_SOURCE)
    case 'p0',      P_J = proj.P0_J;
    case 'p_final', d = load(fullfile(paths.outputs,'P_identified.mat')); P_J = d.P_final;
    case 'manual',  P_J = proj.P0_J;   % edite se quiser
end
constants = struct('m', proj.m, 'g', proj.g);

L = load_log_data(fullfile(paths.data, LOG_FILE));
t_lo = max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]);
t_grid = (t_lo:0.1:t_hi)';
idx = (t_grid>=t_window(1)) & (t_grid<=t_window(2));
time = t_grid(idx);  Nt = numel(time);

pwm = [interp1(L.time_RCOU,L.pwm1_raw,time,'linear'), ...
       interp1(L.time_RCOU,L.pwm2_raw,time,'linear'), ...
       interp1(L.time_RCOU,L.pwm3_raw,time,'linear'), ...
       interp1(L.time_RCOU,L.pwm4_raw,time,'linear')];
pqr_meas = [interp1(L.time_IMU,L.gyrX_raw,time,'linear'), ...
            interp1(L.time_IMU,L.gyrY_raw,time,'linear'), ...
            interp1(L.time_IMU,L.gyrZ_raw,time,'linear')];
att_meas = [interp1(L.time_ATT,L.roll_deg, time,'linear'), ...
            interp1(L.time_ATT,L.pitch_deg,time,'linear'), ...
            interp1(L.time_ATT,L.yaw_deg,  time,'linear')];

%% ===================== pqr DO MODELO (do PWM) =====================
% Modo hybrid integra pqr da dinâmica (entrada = PWM). É o pqr do modelo.
res = sim_window('hybrid', P_J, time, pwm, pqr_meas, att_meas, constants);
pqr_model = [res.p, res.q, res.r];

%% ===================== INTEGRAÇÃO DA ATITUDE =====================
% (2) modelo janela curta   (3) modelo open-loop   (4) MEDIDO janela curta
seg = floor((time - time(1)) / N_RESET);     % índice do segmento por amostra

[phi_w, th_w, psi_w] = integ_euler(time, pqr_model, att_meas, seg);   % modelo, reset
[phi_o, th_o, psi_o] = integ_euler(time, pqr_model, att_meas, []);    % modelo, sem reset
[phi_m, th_m, psi_m] = integ_euler(time, pqr_meas,  att_meas, seg);   % MEDIDO, reset (cinemática pura)

%% ===================== MÉTRICAS =====================
rmse = @(a,b) sqrt(mean((a-b).^2));
wrap = @(e) mod(e+180,360)-180;
fprintf('\n=== RMSE atitude vs ATT do log (janela [%g,%g]s, reset %gs) ===\n', ...
    t_window(1), t_window(2), N_RESET);
fprintf('  %-22s  %7s  %7s  %7s\n','','phi','theta','psi');
fprintf('  %-22s  %7.2f  %7.2f  %7.2f\n','janela curta (modelo)', ...
    rmse(phi_w,att_meas(:,1)), rmse(th_w,att_meas(:,2)), rmse(wrap(psi_w-att_meas(:,3)),0));
fprintf('  %-22s  %7.2f  %7.2f  %7.2f\n','open-loop (modelo)', ...
    rmse(phi_o,att_meas(:,1)), rmse(th_o,att_meas(:,2)), rmse(wrap(psi_o-att_meas(:,3)),0));
fprintf('  %-22s  %7.2f  %7.2f  %7.2f\n','janela curta (MEDIDO)', ...
    rmse(phi_m,att_meas(:,1)), rmse(th_m,att_meas(:,2)), rmse(wrap(psi_m-att_meas(:,3)),0));
fprintf('\n  MEDIDO janela curta ~0  → cinemática EXATA (sem nada a identificar).\n');
fprintf('  MODELO janela curta     → erro = resíduo do pqr do modelo na janela.\n');
fprintf('  MODELO open-loop >>      → deriva = integração acumulada (não é bug).\n');

%% ===================== PLOT =====================
fig = figure('Name','validate_attitude_windowed','Position',[80 60 1300 800]);
labels = {'\phi (°)','\theta (°)','\psi (°)'};
meas  = {att_meas(:,1), att_meas(:,2), att_meas(:,3)};
win   = {phi_w, th_w, psi_w};
ol    = {phi_o, th_o, psi_o};
winm  = {phi_m, th_m, psi_m};
for i = 1:3
    subplot(3,1,i); hold on; grid on;
    plot(time, meas{i}, 'b-',  'LineWidth',1.4, 'DisplayName','medido (ATT log)');
    plot(time, winm{i}, 'g-',  'LineWidth',1.2, 'DisplayName',sprintf('MEDIDO janela %gs (cinemática)',N_RESET));
    plot(time, win{i},  'r--', 'LineWidth',1.4, 'DisplayName',sprintf('modelo janela %gs',N_RESET));
    plot(time, ol{i},   'm:',  'LineWidth',1.1, 'DisplayName','modelo open-loop');
    ylabel(labels{i});
    if i==1, legend('Location','best'); end
    if i==3, xlabel('t [s]'); end
end
sgtitle(sprintf('Atitude: janela curta (reset %gs) vs open-loop — [%g,%g]s', ...
    N_RESET, t_window(1), t_window(2)));
saveas(fig, fullfile(paths.images,'validate_attitude_windowed.png'));
fprintf('\nFigura: %s\n', fullfile(paths.images,'validate_attitude_windowed.png'));


%% ===================== HELPER: integração de Euler =====================
function [phi, theta, psi] = integ_euler(time, pqr, att_meas, seg)
% Integra a cinemática de Euler. Se 'seg' não vazio, reseta a atitude pra
% att_meas no início de cada segmento (janela curta). Saída em graus.
    Nt = numel(time);
    phi = zeros(Nt,1); theta = zeros(Nt,1); psi = zeros(Nt,1);
    ph = deg2rad(att_meas(1,1)); th = deg2rad(att_meas(1,2)); ps = deg2rad(att_meas(1,3));
    phi(1)=rad2deg(ph); theta(1)=rad2deg(th); psi(1)=rad2deg(ps);
    for k = 1:Nt-1
        dt = time(k+1)-time(k); if dt<=0||~isfinite(dt), dt=0.1; end
        p = pqr(k,1); q = pqr(k,2); r = pqr(k,3);
        ct = cos(th); if abs(ct)<1e-6, ct = sign(ct+1e-12)*1e-6; end
        sp = sin(ph); cp = cos(ph); tt = tan(th);
        ph = ph + dt*(p + (q*sp + r*cp)*tt);
        th = th + dt*(q*cp - r*sp);
        ps = ps + dt*((q*sp + r*cp)/ct);
        % reset no início de novo segmento (janela curta)
        if ~isempty(seg) && seg(k+1) ~= seg(k)
            ph = deg2rad(att_meas(k+1,1));
            th = deg2rad(att_meas(k+1,2));
            ps = deg2rad(att_meas(k+1,3));
        end
        phi(k+1)=rad2deg(ph); theta(k+1)=rad2deg(th); psi(k+1)=rad2deg(ps);
    end
end
