% plot_train_val_map.m — mapa dos trechos de treino e de validação
% =========================================================================
% Uma figura com os três voos usados no trabalho, mostrando onde cada trecho
% de identificação (treino) e de validação está no voo. Nenhum trecho de
% validação toca um de treino, em nenhum dos três.
%
%   painel 1  multirrotor, voo 25/05/2026 (logs_concat, 3 logs emendados)
%   painel 2  asa fixa longitudinal, voo 17/12/2024
%   painel 3  asa fixa látero-direcional, voo 17/12/2024
%
% Uso:  >> plot_train_val_map
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();
CT = [0.20 0.40 0.65];      % treino
CV = [0.85 0.37 0.01];      % validação

%% ------------------------------------------------------- multirrotor
L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
TR_q = [4 24; 25 41; 42 62; 63 99; 100 125];
VL_q = [605 625; 610 630; 550 600; 400 420; 420 440; 440 460; ...
        460 480; 525 545; 545 565; 565 585; 585 605; 630 650];
sig = {L.gyrX_raw, L.gyrY_raw, L.gyrZ_raw, L.accX_raw, L.accY_raw, L.accZ_raw};
rot = {'p [rad/s]','q [rad/s]','r [rad/s]','a_x [m/s²]','a_y [m/s²]','a_z [m/s²]'};
tq = L.time_IMU;

f1 = figure('Position',[30 30 1500 1150],'Color','w'); try, f1.Theme='light'; catch, end
tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
for j = 1:6
    ax = nexttile; hold on; grid on;
    faixas(TR_q, CT); faixas(VL_q, CV);
    plot(tq, sig{j}, '-', 'Color',[0.25 0.25 0.25], 'LineWidth',0.6);
    for k = 1:numel(L.boundaries), xline(L.boundaries(k), 'k:', 'LineWidth',1.2); end
    ylabel(rot{j}); xlim([tq(1) tq(end)]);
    if j <= 3, ylim([-4 4]); else, ylim(quantile(sig{j},[0.001 0.999]) + [-1 1]); end
    if j == 1, legenda(ax); end
    if j == 6, xlabel('t [s]'); end
end
title(tl, sprintf(['MULTIRROTOR, voo 25/05/2026 (3 logs emendados, %.0f s)  —  ' ...
    'treino %d trechos (%.0f s) em azul, validação %d trechos (%.0f s) em laranja'], ...
    tq(end)-tq(1), size(TR_q,1), sum(diff(TR_q,1,2)), size(VL_q,1), sum(diff(VL_q,1,2))), ...
    'FontWeight','bold','FontSize',13);
fn1 = fullfile(paths.images,'train_val_map.png');
exportgraphics(f1, fn1, 'BackgroundColor','white','Resolution',140);

%% ------------------------------------------------------- asa fixa
Fl = fw_load('long');
TR_l = [306.5 310.5; 326.5 330.5; 344.0 348.0; 367.5 371.5];
VL_l = [312 318; 334 340; 352 358; 360 366; 376 382; 396 402];
Ft = fw_load('lat');
TR_t = [303.5 307.5; 312.5 316.5; 322.5 326.5; 333.0 337.0];
VL_t = [340 346; 356 362; 380 386; 408 414; 432 438; 446 452];

f2 = figure('Position',[30 30 1500 1000],'Color','w'); try, f2.Theme='light'; catch, end
tl2 = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
LSET = {Fl, TR_l, VL_l, [250 520], {'de','q','az'}, {'\delta_e','q [rad/s]','a_z [m/s²]'}, 'LONGITUDINAL'; ...
        Ft, TR_t, VL_t, [240 570], {'da','p','ay'}, {'\delta_a','p [rad/s]','a_y [m/s²]'}, 'LÁTERO-DIRECIONAL'};
for c = 1:2
    F = LSET{c,1}; TR = LSET{c,2}; VL = LSET{c,3}; XL = LSET{c,4};
    for j = 1:3
        ax = nexttile; hold on; grid on;
        faixas(TR, CT); faixas(VL, CV);
        plot(F.t, F.(LSET{c,5}{j}), '-', 'Color',[0.25 0.25 0.25], 'LineWidth',0.7);
        ylabel(LSET{c,6}{j}); xlim(XL);
        if j == 1
            ylim([-1.1 1.1]);
            title(sprintf('ASA FIXA %s, voo 17/12/2024  —  treino %d doublets (%.0f s), validação %d trechos (%.0f s), asa fixa %.0f s', ...
                LSET{c,7}, size(TR,1), sum(diff(TR,1,2)), size(VL,1), sum(diff(VL,1,2)), sum(F.V>10)*F.dt), 'FontWeight','bold');
            if c == 1, legenda_fw(ax); end
        end
        if j == 3 && c == 2, xlabel('t [s]'); end
    end
end
title(tl2, 'Asa fixa: trechos de treino e de validação', 'FontWeight','bold','FontSize',13);
fn2 = fullfile(paths.images,'train_val_map_fw.png');
exportgraphics(f2, fn2, 'BackgroundColor','white','Resolution',140);

dd = fullfile(getenv('HOME'),'Desktop','DH_modelo_oficial');
if ~exist(dd,'dir'), mkdir(dd); end
copyfile(fn1, fullfile(dd,'train_val_map.png'));  copyfile(fn2, fullfile(dd,'train_val_map_fw.png'));
fprintf('\n  Figuras: %s\n           %s  (cópias na Mesa)\n', fn1, fn2);

% confere que nenhuma validação toca um treino
checa('multirrotor', TR_q, VL_q);  checa('longitudinal', TR_l, VL_l);  checa('lateral', TR_t, VL_t);

function faixas(W, cor)
    for k = 1:size(W,1)
        xregion(W(k,1), W(k,2), 'FaceColor', cor, 'FaceAlpha', 0.28, 'EdgeColor','none');
    end
end
function legenda_fw(ax)
    h1 = patch(ax, nan, nan, [0.20 0.40 0.65], 'FaceAlpha',0.28, 'EdgeColor','none');
    h2 = patch(ax, nan, nan, [0.85 0.37 0.01], 'FaceAlpha',0.28, 'EdgeColor','none');
    legend([h1 h2], {'treino','validação'}, 'Location','northwest','FontSize',9);
end
function legenda(ax)
    h1 = patch(ax, nan, nan, [0.20 0.40 0.65], 'FaceAlpha',0.28, 'EdgeColor','none');
    h2 = patch(ax, nan, nan, [0.85 0.37 0.01], 'FaceAlpha',0.28, 'EdgeColor','none');
    h3 = plot(ax, nan, nan, 'k:', 'LineWidth',1.2);
    legend([h1 h2 h3], {'treino (identificação)','validação','emenda entre logs'}, ...
        'Location','northwest','FontSize',9);
end
function checa(nome, TR, VL)
    n = 0;
    for i = 1:size(VL,1)
        if any(VL(i,1) < TR(:,2) & VL(i,2) > TR(:,1)), n = n + 1; end
    end
    if n == 0, fprintf('  %-14s: nenhuma janela de validação toca o treino\n', nome);
    else,      fprintf('  %-14s: ATENÇÃO, %d janelas de validação tocam o treino\n', nome, n); end
end
