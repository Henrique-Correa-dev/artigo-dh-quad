% fw_train_val_map.m — trechos de identificação e de validação nos voos de asa fixa
% =========================================================================
% Mostra os dois logs de asa fixa de 17/12/2024 por inteiro, marcando em azul
% os doublets usados na identificação e em laranja os reservados para
% validação, que não participam da estimação.
%
%   longitudinal   3 doublets de profundor no treino, 1 na validação
%   látero-direc.  4 doublets de aileron no treino, 1 na validação
%                  (o de validação é um dos doublets de amplitude cheia que
%                   sobraram do voo, nunca usado em nenhuma identificação)
%
% Uso:  >> fw_train_val_map
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();
CT = [0.20 0.40 0.65];   CV = [0.85 0.37 0.01];

TR_l = [306.5 310.5; 344.0 348.0; 367.5 371.5];    VL_l = [326.5 330.5];
TR_t = [303.5 307.5; 312.5 316.5; 322.5 326.5; 333.0 337.0];   VL_t = [343.0 346.0];

Fl = fw_load('long');   Ft = fw_load('lat');
f = figure('Position',[30 30 1500 900],'Color','w');  try, f.Theme='light'; catch, end
tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
SET = {Fl, TR_l, VL_l, [250 520], {'de','q','az'}, {'\delta_e','q [rad/s]','a_z [m/s²]'}, 'LONGITUDINAL', 'profundor'; ...
       Ft, TR_t, VL_t, [240 520], {'da','p','ay'}, {'\delta_a','p [rad/s]','a_y [m/s²]'}, 'LÁTERO-DIRECIONAL', 'aileron'};
for c = 1:2
    F = SET{c,1};  TR = SET{c,2};  VL = SET{c,3};  XL = SET{c,4};
    for j = 1:3
        ax = nexttile; hold on; grid on;
        for k = 1:size(TR,1)
            xregion(TR(k,1), TR(k,2), 'FaceColor',CT, 'FaceAlpha',0.30, 'EdgeColor','none');
        end
        for k = 1:size(VL,1)
            xregion(VL(k,1), VL(k,2), 'FaceColor',CV, 'FaceAlpha',0.30, 'EdgeColor','none');
        end
        plot(F.t, F.(SET{c,5}{j}), '-', 'Color',[0.25 0.25 0.25], 'LineWidth',0.8);
        ylabel(SET{c,6}{j});  xlim(XL);
        if j == 1
            ylim([-1.1 1.1]);
            title(sprintf('ASA FIXA %s, voo 17/12/2024  —  %d doublets de %s na identificação, %d na validação', ...
                SET{c,7}, size(TR,1), SET{c,8}, size(VL,1)), 'FontWeight','bold');
            h1 = patch(ax, nan, nan, CT, 'FaceAlpha',0.30, 'EdgeColor','none');
            h2 = patch(ax, nan, nan, CV, 'FaceAlpha',0.30, 'EdgeColor','none');
            legend([h1 h2], {'identificação','validação'}, 'Location','northeast','FontSize',9);
        end
        if j == 3, xlabel('t [s]'); end
    end
end
title(tl, 'Asa fixa: trechos de identificação e de validação', 'FontWeight','bold','FontSize',13);
fn = fullfile(paths.images,'fw_train_val_map.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',150);
dd = fullfile(getenv('HOME'),'Desktop','DH_asafixa_modelos');
if ~exist(dd,'dir'), mkdir(dd); end
copyfile(fn, fullfile(dd,'fw_train_val_map.png'));
fprintf('\n  Figura: %s  (cópia na Mesa)\n', fn);
% confere isolamento
for c = 1:2
    TR = SET{c,2};  VL = SET{c,3};  n = 0;
    for i = 1:size(VL,1)
        if any(VL(i,1) < TR(:,2) & VL(i,2) > TR(:,1)), n = n + 1; end
    end
    fprintf('  %-18s: %d janela(s) de validação tocando o treino\n', SET{c,7}, n);
end
