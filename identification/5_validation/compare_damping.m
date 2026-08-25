% compare_damping.m — Figura do teste de amortecimento: constante × ∝V × híbrido
% =========================================================================
% Simula (modo full, mesma cadeia de motor) os três modelos identificados em
% outputs/diag/{unified,damp_aero,damp_hybrid} na janela de validação e sobrepõe
% p, q, r medidos. Para os modos ∝V, o gancho de vtol_dynamics (appdata
% 'diag_damp') é ativado com a mesma V(t) usada na identificação.
%
% Uso:  >> compare_damping
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
T_WINDOW = [610 630];  MODE = 'full';

L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';  idx = tg>=T_WINDOW(1) & tg<=T_WINDOW(2);  time = tg(idx);
ip = @(tt,xx) interp1(tt, xx, time, 'linear');
pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
pqr_meas = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
att_meas = [ip(L.time_ATT,L.roll_deg), ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
constants = struct('m', proj.m, 'g', proj.g);

% velocidade estimada (mesmo método da identificação ∝V)
VEst = estimate_velocity(L, tg);  Vgrid = VEst.V;  Vgrid(~isfinite(Vgrid)) = 0;

C = struct('tag',{'unified','damp_aero','damp_hybrid'}, 'mode',{'const','aero','hybrid'}, ...
           'c0',{[0;0;0],[0;0;0],[1.19;1.20;0.08]}, ...
           'label',{'constante (c_p, c_q, c_r)','só \propto V (C_{lp}, C_{mq}, C_{nr})','híbrido (rotor fixo + \propto V)'}, ...
           'cor',{[0 0.45 0.7],[0.85 0.37 0.01],[0.3 0.6 0.3]});
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
for k = 1:numel(C)
    P = load(fullfile(paths.outputs,'diag',C(k).tag,'P_identified.mat')).P_final(:);
    if strcmp(C(k).mode,'const'), if isappdata(0,'diag_damp'), rmappdata(0,'diag_damp'); end
    else, setappdata(0,'diag_damp', struct('mode',C(k).mode,'t',tg,'V',Vgrid,'c0',C(k).c0)); end
    r = sim_window(MODE, P, time, pwm, pqr_meas, att_meas, constants);
    C(k).res = r;  C(k).R2 = [R2(pqr_meas(:,1),r.p), R2(pqr_meas(:,2),r.q), R2(pqr_meas(:,3),r.r)];
    fprintf('%-12s R² p %.3f | q %.3f | r %.3f\n', C(k).tag, C(k).R2);
end
if isappdata(0,'diag_damp'), rmappdata(0,'diag_damp'); end

f = figure('Position',[60 60 1150 760],'Color','w'); try, f.Theme='light'; catch, end
lab = {'p','q','r'};
for i = 1:3
    subplot(4,1,i); hold on; grid on;
    h = plot(time, pqr_meas(:,i), 'k-', 'LineWidth', 1.8);
    for k = 1:numel(C)
        Y = [C(k).res.p, C(k).res.q, C(k).res.r];
        h(end+1) = plot(time, Y(:,i), '-', 'Color', C(k).cor, 'LineWidth', 1.3); %#ok<SAGROW>
    end
    ylabel([lab{i} ' [rad/s]']);
    if i == 1
        lg = arrayfun(@(c) sprintf('%s  (R² %.2f/%.2f/%.2f)', c.label, c.R2), C, 'UniformOutput', false);
        legend(h, [{'medido'}, lg], 'Location','northoutside', 'NumColumns', 2, 'FontSize', 9);
    end
end
subplot(4,1,4); hold on; grid on;
plot(time, Vgrid(idx), 'k-', 'LineWidth', 1.3); ylabel('V estimada [m/s]'); xlabel('tempo [s]');
text(0.005, 1.08, 'velocidade estimada (atitude + barômetro) que escala os termos \propto V', 'Units','normalized');
exportgraphics(f, fullfile(paths.images,'compare_damping.png'), 'BackgroundColor','white','Resolution',150);
fprintf('  Figura: %s\n', fullfile(paths.images,'compare_damping.png'));
