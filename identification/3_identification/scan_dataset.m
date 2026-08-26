% scan_dataset.m — Campanha de voo e seleção das janelas (treino / validação)
% =========================================================================
% Campanha COMPLETA em quatro sinais — RCIN | PWM | taxas | ACELERÔMETRO — com
% as janelas escolhidas destacadas, rotuladas e separadas por linhas pontilhadas:
%   verde  = TREINO    (eixos isolados):     DT1 0-40, DT2 40-80, DT3 80-125
%   laranja= VALIDAÇÃO (movimento composto): DT4 450-470, DT5 610-630 (janela principal)
% Só diagnóstico — não modifica nada.  Uso:  >> scan_dataset

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---------------- CONFIG ----------------
LOG_FILE = 'logs_concat.mat';
dt = 0.1;

% Todas as janelas em ordem cronológica; ISVAL marca as de validação.
WINS  = [  0  40;  40  80;  80 125; 450 470; 610 630];
ISVAL = logical([0 0 0 1 1]);
LBL   = {'\DeltaT_1','\DeltaT_2','\DeltaT_3','\DeltaT_4','\DeltaT_5'};

%% ---------------- CARREGAR + GRADE ----------------
L  = load_log_data(fullfile(paths.data, LOG_FILE));
t0 = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t1 = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t0:dt:t1)';
ip = @(tt, xx) interp1(tt, xx, tg, 'linear');

p  = ip(L.time_IMU, L.gyrX_raw);  q = ip(L.time_IMU, L.gyrY_raw);  r = ip(L.time_IMU, L.gyrZ_raw);
ax = ip(L.time_IMU, L.accX_raw);  ay= ip(L.time_IMU, L.accY_raw);  az= ip(L.time_IMU, L.accZ_raw);
W  = [ip(L.time_RCOU, L.pwm1_raw), ip(L.time_RCOU, L.pwm2_raw), ...
      ip(L.time_RCOU, L.pwm3_raw), ip(L.time_RCOU, L.pwm4_raw)];

hasRC = isfield(L,'time_RCIN') && ~isempty(L.time_RCIN);
if hasRC
    ipr = @(xx) interp1(L.time_RCIN, xx, tg, 'linear', 'extrap');
    RC  = [ipr(L.rcin_roll), ipr(L.rcin_pitch), ipr(L.rcin_yaw)];
end

gap = false(size(tg));
for j = find(diff(L.time_IMU) > 1.0)'
    gap = gap | (tg > L.time_IMU(j) & tg < L.time_IMU(j+1));
end

%% ---------------- CORTE 150–400 s (trecho sem manobras aproveitadas) ----------------
CUT   = [150 400];
T_END = 650;                                     % descarta o trecho final (sem manobras aproveitadas)
shift = diff(CUT);
mask  = (tg <= CUT(1)) | (tg >= CUT(2) & tg <= T_END);
k     = find(diff(tg(mask)) > 2*dt, 1);          % índice da emenda entre os blocos
mapt  = @(t) t - shift.*(t >= CUT(2));           % tempo real -> tempo plotado

tp = mapt(cut_insert(tg, mask, k));
p  = cut_insert(p,  mask, k);  q  = cut_insert(q,  mask, k);  r  = cut_insert(r,  mask, k);
ax = cut_insert(ax, mask, k);  ay = cut_insert(ay, mask, k);  az = cut_insert(az, mask, k);
W  = cut_insert(W,  mask, k);
if hasRC, RC = cut_insert(RC, mask, k); end
gp = cut_insert(double(gap), mask, k) > 0.5;
WINSp = mapt(WINS);

%% ---------------- FIGURA ----------------
FS_AX  = 14;   % eixos (ticks + ylabel/xlabel)
FS_LEG = 13;   % legendas
FS_TAG = 15;   % rótulos \DeltaT_i
FS_TIT = 16;   % título geral
LW     = 1.0;  % espessura das linhas de sinal

fig = figure('Name','scan_dataset', 'Position',[40 40 1150 870]);
set(fig, 'Color', 'w', 'DefaultAxesFontSize', FS_AX, 'DefaultLineLineWidth', LW);
try, fig.Theme = 'light'; catch, end
A = gobjects(4,1);

A(1) = subplot(4,1,1); hold on; grid on;
if hasRC
    plot(tp,RC(:,1),'b'); plot(tp,RC(:,2),'r'); plot(tp,RC(:,3),'Color',[0 .6 0]);
    legend('roll','pitch','yaw','Location','eastoutside','FontSize',FS_LEG); ylabel('RCIN (\mus)');
end
mark_wins(WINSp,ISVAL); shade_gaps(tp,gp); label_wins(WINSp,LBL,FS_TAG);

A(2) = subplot(4,1,2); hold on; grid on;
plot(tp,W(:,1)); plot(tp,W(:,2)); plot(tp,W(:,3)); plot(tp,W(:,4));
legend('M1','M2','M3','M4','Location','eastoutside','FontSize',FS_LEG); ylabel('PWM (\mus)');
mark_wins(WINSp,ISVAL); shade_gaps(tp,gp);

A(3) = subplot(4,1,3); hold on; grid on;
plot(tp,p,'b'); plot(tp,q,'r'); plot(tp,r,'Color',[0 .6 0]);
legend('p','q','r','Location','eastoutside','FontSize',FS_LEG); ylabel('taxas (rad/s)');
mark_wins(WINSp,ISVAL); shade_gaps(tp,gp);

A(4) = subplot(4,1,4); hold on; grid on;
plot(tp,ax,'b'); plot(tp,ay,'r'); plot(tp,az,'Color',[0 .6 0]);
legend('a_x','a_y','a_z','Location','eastoutside','FontSize',FS_LEG); ylabel('acel. (m/s^2)'); xlabel('t (s)');
mark_wins(WINSp,ISVAL); shade_gaps(tp,gp);

linkaxes(A,'x');  xlim([0 mapt(T_END)]);
XT  = [0 50 100 200 250 300 350 400];
XTL = {'0','50','100','450','500','550','600','650'};
for i = 1:4
    set(A(i), 'XTick', XT, 'XTickLabel', XTL);
    xline(A(i), CUT(1), '--', 'Color',[.3 .3 .3], 'LineWidth',1.2, 'HandleVisibility','off');
end
sgtitle('Campanha de voo — verde: treino (eixos isolados)  |  laranja: validação (movimento composto)', ...
        'FontSize', FS_TIT);
exportgraphics(fig, fullfile(paths.images,'scan_dataset.png'), 'BackgroundColor','white', 'Resolution',200);
fprintf('\n  Figura salva: %s\n', fullfile(paths.images,'scan_dataset.png'));

% Copia direto para a figura usada na dissertação (Fig. "Campanha de voo completa")
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tese_png = fullfile(root, 'dissertacao', 'Cap5', 'excitacao_dataset.png');
exportgraphics(fig, tese_png, 'BackgroundColor','white', 'Resolution',200);
fprintf('  Figura salva: %s\n', tese_png);


%% ======================= HELPERS =======================
function mark_wins(WINS, ISVAL)
    yl = ylim;
    for i = 1:size(WINS,1)
        if ISVAL(i), c = [1 .75 .2]; a = 0.32; else, c = [.55 1 .55]; a = 0.22; end
        patch([WINS(i,1) WINS(i,2) WINS(i,2) WINS(i,1)], [yl(1) yl(1) yl(2) yl(2)], ...
              c, 'FaceAlpha',a, 'EdgeColor','none', 'HandleVisibility','off');
    end
    for e = unique(WINS(:))'        % linhas pontilhadas separando os trechos
        xline(e, 'LineStyle',':', 'Color',[.25 .25 .25], 'LineWidth',0.9, 'HandleVisibility','off');
    end
    ylim(yl);
end

function label_wins(WINS, LBL, fs)
    % Rótulos ACIMA do eixo, com colchete horizontal marcando a extensão da janela.
    if nargin < 3, fs = 9; end
    yl = ylim;  h = yl(2) - yl(1);
    yb   = yl(2) + 0.03*h;    % altura do colchete
    yt   = yl(2) + 0.115*h;   % altura do texto
    tick = 0.04*h;            % perninhas do colchete
    for i = 1:size(WINS,1)
        x1 = WINS(i,1);  x2 = WINS(i,2);
        line([x1 x2], [yb yb], 'Color','k', 'LineWidth',1.1, ...
             'Clipping','off', 'HandleVisibility','off');
        line([x1 x1], [yb yb-tick], 'Color','k', 'LineWidth',1.1, ...
             'Clipping','off', 'HandleVisibility','off');
        line([x2 x2], [yb yb-tick], 'Color','k', 'LineWidth',1.1, ...
             'Clipping','off', 'HandleVisibility','off');
        text(mean(WINS(i,:)), yt, LBL{i}, 'HorizontalAlignment','center', ...
             'FontSize',fs, 'FontWeight','bold', 'Clipping','off');
    end
    ylim(yl);
end

function shade_gaps(tg, gap)
    yl = ylim;  d = diff([0; gap(:); 0]);
    s = tg(d(1:end-1) == 1);  e = tg(d(2:end) == -1);
    for i = 1:min(numel(s), numel(e))
        patch([s(i) e(i) e(i) s(i)], [yl(1) yl(1) yl(2) yl(2)], [.6 .6 .6], ...
              'FaceAlpha',0.5, 'EdgeColor','none', 'HandleVisibility','off');
    end
    ylim(yl);
end

function y = cut_insert(x, mask, k)
    % aplica a máscara e insere uma linha de NaN na emenda (quebra a linha do plot)
    y = x(mask, :);
    y = [y(1:k,:); nan(1, size(y,2)); y(k+1:end,:)];
end
