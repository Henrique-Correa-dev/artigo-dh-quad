% compare_dissertacao.m — Os valores da Tabela 5.3 rodam? (só simula, não otimiza)
% =========================================================================
% MOTIVO
%   A Tabela 5.3 da dissertação e o P_identified.mat do programa DIVERGEM nas
%   inércias: a tabela traz Jy = 0,0984, Jz = 0,1453 e Jxz = 0,00261, enquanto o
%   programa usa Jy = 0,0934 com Jz e Jxz TRAVADOS no CAD (0,1262 e 0,00157).
%   cp, cq, cr, c_T e c_M batem entre os dois. Este script simula os dois vetores
%   na MESMA janela de validação (605–625 s) e mede R², RMSE e TIC, para decidir
%   qual é a versão oficial com base no ajuste, e não na memória.
%
%   Observação física: na tabela, Jx + Jy = 0,0469 + 0,0984 = 0,1453 = Jz, ou
%   seja, o valor está EXATAMENTE sobre a desigualdade triangular (corpo planar).
%   O script reporta a folga triangular de cada conjunto.
%
% CONJUNTOS COMPARADOS
%   diss_T0     Θ0  da Tabela 5.3 (CAD + k=1 + cp=cq=cr=0,5)
%   prior_T0    Θ0  novo, com o amortecimento a priori (prior_damping.m)
%   diss_EEM    Θ_EEM da Tabela 5.3
%   diss_OEM    Θ_OEM da Tabela 5.3  ← o que está escrito na dissertação
%   prog_OEM    P_final de outputs/P_identified.mat  ← o que roda no programa
%
%   Escalas da tabela: c_T em 1e-7 N/RPM² e c_M em 1e-9 N·m/RPM². No código,
%   k_T e k_Q são multiplicadores sobre kT_rpm = 1,764e-7 e kQ_rpm = 5,04e-9
%   (motor_models('physical')), então converte-se dividindo pelas referências.
%
% Uso:  >> compare_dissertacao
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

LOG_FILE  = 'logs_concat.mat';
T_WINDOW  = [610, 630];      % mesma janela de validação do identify_plant
DELAY_PWM = 0;               % o atraso está DENTRO de sim_window (motor_chain); 0 aqui
MODE      = 'full';          % 'full' = referência (Tabela 5.5); 'hybrid' isola o rotacional

kT_ref = 1.764e-7;  kQ_ref = 5.04e-9;   % referências de motor_models('physical')
sc = @(cT, cM) [cT(:)*1e-7/kT_ref; cM(:)*1e-9/kQ_ref];   % tabela → k_T, k_Q

%% ---------------- CONJUNTOS DE PARÂMETROS ----------------
proj = parameters();
prog = load(fullfile(paths.outputs,'P_identified.mat'));

S(1).tag = 'diss_T0';   S(1).P = [0.0432; 0.0844; 0.1262; 0.00157; sc([1.76 1.76 1.76 1.76],[5.04 5.04 5.04 5.04]); 0.5; 0.5; 0.5];
S(2).tag = 'prior_T0';  S(2).P = proj.P0_J(:);
S(3).tag = 'diss_EEM';  S(3).P = [0.0526; 0.0804; 0.1131; 0.00047; sc([0.89 0.88 0.86 0.88],[3.33 2.62 4.28 2.48]); 1.63; 1.50; 0.58];
S(4).tag = 'diss_OEM';  S(4).P = [0.0469; 0.0984; 0.1453; 0.00261; sc([1.87 1.86 1.86 1.89],[3.38 3.16 3.92 3.47]); 6.85; 4.68; 0.81];
S(5).tag = 'prog_OEM';  S(5).P = prog.P_final(:);

%% ---------------- DADOS DA JANELA ----------------
L = load_log_data(fullfile(paths.data, LOG_FILE));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg   = (t_lo:0.1:t_hi)';
idx  = (tg >= T_WINDOW(1)) & (tg <= T_WINDOW(2));
time = tg(idx);
ip   = @(tt,xx) interp1(tt, xx, time, 'linear');

pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ...
       ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
if DELAY_PWM > 0
    dl = @(x) [repmat(x(1),DELAY_PWM,1); x(1:end-DELAY_PWM)];
    for ci = 1:4, pwm(:,ci) = dl(pwm(:,ci)); end
end
pqr_meas = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
acc_meas = [ip(L.time_IMU,L.accX_raw), ip(L.time_IMU,L.accY_raw), ip(L.time_IMU,L.accZ_raw)];
att_meas = [ip(L.time_ATT,L.roll_deg),  ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
constants = struct('m', proj.m, 'g', proj.g);

fprintf('\n=== Janela %g–%g s | modo %s | %d amostras ===\n', T_WINDOW, MODE, numel(time));

%% ---------------- SIMULAÇÃO E MÉTRICAS ----------------
sig = {'p','q','r','accX','accY','accZ'};
for k = 1:numel(S)
    P = S(k).P;
    S(k).tri = (P(1)+P(2)) - P(3);          % folga triangular Jz <= Jx+Jy
    r = sim_window(MODE, P, time, pwm, pqr_meas, att_meas, constants);
    ms = [pqr_meas, acc_meas];
    sm = [r.p, r.q, r.r, r.accX, r.accY, r.accZ];
    for i = 1:6
        m = fit_metrics(ms(:,i), sm(:,i));
        S(k).R2(i) = m.R2;  S(k).TIC(i) = m.TIC;  S(k).RMSE(i) = m.RMSE;
    end
    S(k).res = r;
end

fprintf('\n%-10s | %7s %7s %7s | %7s %7s %7s | %9s\n', 'conjunto', 'R2 p','R2 q','R2 r','TIC p','TIC q','TIC r','Jx+Jy-Jz');
fprintf('%s\n', repmat('-',1,72));
for k = 1:numel(S)
    fprintf('%-10s | %7.3f %7.3f %7.3f | %7.3f %7.3f %7.3f | %9.5f\n', ...
        S(k).tag, S(k).R2(1:3), S(k).TIC(1:3), S(k).tri);
end
fprintf('\n%-10s | %7s %7s %7s | %7s %7s %7s\n', 'conjunto', 'R2 aX','R2 aY','R2 aZ','TIC aX','TIC aY','TIC aZ');
fprintf('%s\n', repmat('-',1,62));
for k = 1:numel(S)
    fprintf('%-10s | %7.3f %7.3f %7.3f | %7.3f %7.3f %7.3f\n', S(k).tag, S(k).R2(4:6), S(k).TIC(4:6));
end

fprintf('\n  Inércias por conjunto (kg·m²):\n');
fprintf('  %-10s %8s %8s %8s %9s   %s\n','conjunto','Jx','Jy','Jz','Jxz','folga triangular');
for k = 1:numel(S)
    P = S(k).P;
    if S(k).tri < 1e-6, nota = '<< NO LIMITE (corpo planar)'; else, nota = ''; end
    fprintf('  %-10s %8.4f %8.4f %8.4f %9.5f   %+8.5f %s\n', S(k).tag, P(1), P(2), P(3), P(4), S(k).tri, nota);
end

%% ---------------- FIGURA ----------------
f = figure('Position',[60 60 1150 680],'Color','w');
try, f.Theme = 'light'; catch, end
% ordem de plotagem: Θ0 (os dois) primeiro, atrás; depois os dois OEM
PLOT_K  = [4 5];            % só os dois OEM (dissertação × programa); [1 2 4 5] inclui os Θ0
cores   = [0.70 0.70 0.70; 0.45 0.45 0.45; 0 0 0; 0.85 0.37 0.01; 0.00 0.45 0.70];
estilos = {':','--','-','-','-'};
larg    = [1.2, 1.2, 1.6, 1.6, 1.4];
rotulos = {'\Theta_0 dissertação (c = 0,5)', '\Theta_0 a priori (influxo)', '', ...
           '\Theta_{OEM} dissertação', '\Theta_{OEM} programa'};
for i = 1:3
    subplot(3,1,i); hold on; grid on;
    h = plot(time, pqr_meas(:,i), 'k-', 'LineWidth', 1.8);
    for k = PLOT_K
        Y = [S(k).res.p, S(k).res.q, S(k).res.r];
        h(end+1) = plot(time, Y(:,i), estilos{k}, 'Color', cores(k,:), 'LineWidth', larg(k)); %#ok<SAGROW>
    end
    % escala pelo sinal medido: as curvas de Θ0 saem do eixo, e é esse o ponto
    if any(PLOT_K <= 2)   % com Θ0 na figura, trava a escala no medido (Θ0 sai do eixo)
        yr = [min(pqr_meas(:,i)), max(pqr_meas(:,i))];  dy = 0.35*diff(yr);
        ylim([yr(1)-dy, yr(2)+dy]);
    end
    ylabel([sig{i} ' [rad/s]']);
    if i == 1
        legend(h, [{'medido'}, rotulos(PLOT_K)], 'Orientation','horizontal', ...
               'Location','northoutside', 'NumColumns',3, 'FontSize',9);
        title(sprintf('Validação %g–%g s, modo %s: \\Theta_{OEM} da Tabela 5.3 × \\Theta_{OEM} do programa', T_WINDOW, MODE));
    end
    if i == 3, xlabel('tempo [s]'); end
end
exportgraphics(f, fullfile(paths.images,sprintf('compare_dissertacao_%s.png',MODE)), 'BackgroundColor','white','Resolution',150);
fprintf('\n  Figura: %s\n', fullfile(paths.images,sprintf('compare_dissertacao_%s.png',MODE)));

%% ---------------- FIGURA 2: acelerômetros ----------------
f2 = figure('Position',[60 60 1150 680],'Color','w');
try, f2.Theme = 'light'; catch, end
acc_lab = {'a_x','a_y','a_z'};
for i = 1:3
    subplot(3,1,i); hold on; grid on;
    h = plot(time, acc_meas(:,i), 'k-', 'LineWidth', 1.8);
    for k = PLOT_K
        Y = [S(k).res.accX, S(k).res.accY, S(k).res.accZ];
        h(end+1) = plot(time, Y(:,i), estilos{k}, 'Color', cores(k,:), 'LineWidth', larg(k)); %#ok<SAGROW>
    end
    if any(PLOT_K <= 2)
        yr = [min(acc_meas(:,i)), max(acc_meas(:,i))];  dy = 0.35*diff(yr);
        ylim([yr(1)-dy, yr(2)+dy]);
    end
    ylabel([acc_lab{i} ' [m/s^2]']);
    if i == 1
        legend(h, [{'medido'}, rotulos(PLOT_K)], 'Orientation','horizontal', ...
               'Location','northoutside', 'NumColumns',3, 'FontSize',9);
        title(sprintf('Validação %g–%g s, modo %s: acelerômetros (força específica na IMU)', T_WINDOW, MODE));
    end
    if i == 3, xlabel('tempo [s]'); end
end
exportgraphics(f2, fullfile(paths.images,sprintf('compare_dissertacao_%s_acc.png',MODE)), 'BackgroundColor','white','Resolution',150);
fprintf('  Figura: %s\n', fullfile(paths.images,sprintf('compare_dissertacao_%s_acc.png',MODE)));
