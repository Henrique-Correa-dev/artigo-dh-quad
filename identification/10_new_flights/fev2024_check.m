% fev2024_check.m — por que a identificação de fev/2024 não fecha
% =========================================================================
% Simula (modo full, mesma cadeia de motor da validação) a janela de treino e a
% de validação do voo de 17/02/2024 com os parâmetros identificados, e mostra
% ao lado o comando de PWM com o limite de saturação do ArduPilot
% (MOT_SPIN_MAX = 0,95 → 1950 µs). Serve para separar erro de modelo de dado
% inutilizável por saturação.
%
% Uso:  >> fev2024_check
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();
parameters('campaign','fev2024');  clear aero_gains
proj = parameters();
PWM_SAT = 1950;

FEV = fullfile(fileparts(paths.root), 'babyshark_vtol_model', 'examples', ...
               'example_inputs', 'Logs do Voo Vertical', 'VooVert_001.mat');
P = load(fullfile(paths.outputs,'runs','fev2024_cg','P_identified.mat')).P_final(:);
L = load_log_data(FEV);
tg = (max([L.time_IMU(1) L.time_ATT(1) L.time_RCOU(1)]):0.1:min([L.time_IMU(end) L.time_ATT(end) L.time_RCOU(end)]))';
ip = @(tt,xx,t) interp1(tt, xx, t, 'linear');
constants = struct('m', proj.m, 'g', proj.g);
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

WINS = [161 171; 191 211];  nomes = {'treino (161–171 s)','validação (191–211 s)'};
f = figure('Position',[40 40 1350 900],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(4, 2, 'TileSpacing','compact','Padding','compact');
lab = {'p','q','r'};
for c = 1:2
    tw = WINS(c,:);  time = tg(tg>=tw(1) & tg<=tw(2));
    pwm = [ip(L.time_RCOU,L.pwm1_raw,time), ip(L.time_RCOU,L.pwm2_raw,time), ...
           ip(L.time_RCOU,L.pwm3_raw,time), ip(L.time_RCOU,L.pwm4_raw,time)];
    pqr = [ip(L.time_IMU,L.gyrX_raw,time), ip(L.time_IMU,L.gyrY_raw,time), ip(L.time_IMU,L.gyrZ_raw,time)];
    att = [ip(L.time_ATT,L.roll_deg,time), ip(L.time_ATT,L.pitch_deg,time), ip(L.time_ATT,L.yaw_deg,time)];
    r = sim_window('full', P, time, pwm, pqr, att, constants);
    Y = [r.p, r.q, r.r];
    for i = 1:3
        nexttile((i-1)*2 + c); hold on; grid on;
        plot(time, pqr(:,i), 'k-', 'LineWidth', 1.6);
        plot(time, Y(:,i), '-', 'Color', [0.85 0.37 0.01], 'LineWidth', 1.2);
        ylabel([lab{i} ' [rad/s]']); xlim(tw);
        text(0.01, 0.94, sprintf('R^2 = %.2f', R2(pqr(:,i), Y(:,i))), 'Units','normalized', ...
            'VerticalAlignment','top', 'FontSize',9, 'BackgroundColor',[1 1 1 0.75]);
        if i == 1
            text(0.0, 1.10, nomes{c}, 'Units','normalized','FontWeight','bold','FontSize',10);
            if c == 1, legend({'medido','simulado'},'Location','southwest','FontSize',8); end
        end
    end
    nexttile(6 + c); hold on; grid on;
    plot(time, pwm, 'LineWidth', 0.9);
    yline(PWM_SAT, 'r--', 'MOT\_SPIN\_MAX = 1950 \mus', 'LineWidth',1.2, 'LabelHorizontalAlignment','left');
    sat = 100*mean(any(pwm >= PWM_SAT-2, 2));
    text(0.01, 0.94, sprintf('%.0f%% do tempo com algum motor saturado', sat), 'Units','normalized', ...
        'VerticalAlignment','top','FontSize',9,'BackgroundColor',[1 1 1 0.75]);
    ylabel('PWM [\mus]'); xlabel('t [s]'); xlim(tw); ylim([1300 2050]);
    if c == 1, legend({'M1','M2','M3','M4'},'Location','southwest','Orientation','horizontal','FontSize',8); end
end
title(tl, 'Voo de 17/02/2024 — modelo identificado nesta campanha (CG corrigido pelo equilíbrio)');
fn = fullfile(paths.images,'new_flights','fev2024_check.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',130);
fprintf('  Figura: %s\n', fn);

% quanto tempo cada motor passa saturado, no voo inteiro
ii = L.time_RCOU >= 141 & L.time_RCOU <= 211;
W = [L.pwm1_raw(ii), L.pwm2_raw(ii), L.pwm3_raw(ii), L.pwm4_raw(ii)];
fprintf('\n  Saturação por motor (141–211 s, limite %d µs):\n', PWM_SAT);
for j = 1:4
    fprintf('    M%d: média %.0f µs | %.1f%% das amostras no limite\n', j, mean(W(:,j)), 100*mean(W(:,j) >= PWM_SAT-2));
end
fprintf('    algum motor no limite: %.1f%% das amostras\n', 100*mean(any(W >= PWM_SAT-2, 2)));
parameters('campaign','');
