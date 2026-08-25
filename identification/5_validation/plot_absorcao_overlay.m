% plot_absorcao_overlay.m — com e sem aerodinâmica no MESMO gráfico
% =========================================================================
% Mesmos dois modelos do plot_absorcao.m (moment_win6_2026 sem aero e
% oficial_2026_final com aero), sobrepostos por janela de validação, com o
% medido em preto. R² dos dois em cada painel. Salva em
%   outputs/images/absorcao/sobrepostos/
% Uso:  >> plot_absorcao_overlay
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
WINS = [605 625; 610 630; 550 600; 400 420; 420 440; 440 460; ...
        460 480; 525 545; 545 565; 565 585; 585 605; 630 650];
W_ = getappdata(0,'absorcao_wins');  if ~isempty(W_), WINS = W_; end   % regenerar só algumas
RUNS = {'moment_win6_2026',  '';  'oficial_2026_final','measured'};
NOM  = {'sem aerodinâmica','com aerodinâmica'};
COR  = [0.20 0.40 0.65; 0.85 0.37 0.01];

L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';
set_aero_vel(L, tg);
setappdata(0,'damp_form_override','moment');  clear aero_gains
gapA = L.boundaries(end);  gapB = L.log_starts(end);
constants = struct('m', proj.m, 'g', proj.g);
R2  = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
TIC = @(y,yh) sqrt(mean((y-yh).^2)) / (sqrt(mean(y.^2)) + sqrt(mean(yh.^2)));
canais = {'p','q','r','accX','accY','accZ'};
rot = {'p [rad/s]','q [rad/s]','r [rad/s]','a_x [m/s²]','a_y [m/s²]','a_z [m/s²]'};
DIR = fullfile(paths.images,'absorcao','sobrepostos');
if ~exist(DIR,'dir'), mkdir(DIR); end
P1 = load(fullfile(paths.outputs,'runs',RUNS{1,1},'P_identified.mat')).P_final(:);
P2 = load(fullfile(paths.outputs,'runs',RUNS{2,1},'P_identified.mat')).P_final(:);

for w = 1:size(WINS,1)
    tw = WINS(w,:);
    if tw(1) < gapB && tw(2) > gapA, continue; end
    idx = tg>=tw(1) & tg<=tw(2);  time = tg(idx);
    ip = @(tt,xx) interp1(tt, xx, time, 'linear');
    pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ...
           ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
    pqr = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
    acc = [ip(L.time_IMU,L.accX_raw), ip(L.time_IMU,L.accY_raw), ip(L.time_IMU,L.accZ_raw)];
    att = [ip(L.time_ATT,L.roll_deg), ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
    med = [pqr, acc];
    if isappdata(0,'faero_on'), rmappdata(0,'faero_on'); end
    r1 = sim_window('full', P1, time, pwm, pqr, att, constants);
    setappdata(0,'faero_on','measured');
    r2 = sim_window('full', P2, time, pwm, pqr, att, constants);
    rmappdata(0,'faero_on');
    f = figure('Position',[30 30 1150 1250],'Color','w','Visible','off');
    try, f.Theme='light'; catch, end
    tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
    for j = 1:6
        nexttile; hold on; grid on;
        y1 = r1.(canais{j});  y2 = r2.(canais{j});
        plot(time, med(:,j), 'k-', 'LineWidth', 1.8);
        plot(time, y1, '-', 'Color',COR(1,:), 'LineWidth', 1.2);
        plot(time, y2, '-', 'Color',COR(2,:), 'LineWidth', 1.2);
        ylabel(rot{j}); xlim(tw);
        % folga de 25% no topo do eixo: as métricas moram ali, sem cobrir as curvas
        yl = ylim;  ylim([yl(1), yl(2) + 0.25*(yl(2)-yl(1))]);
        c1 = sprintf('\\color[rgb]{%.2f %.2f %.2f}', COR(1,:));
        c2 = sprintf('\\color[rgb]{%.2f %.2f %.2f}', COR(2,:));
        s1 = sprintf('{\\bfsem aero:}  {\\bfR^2} %s%.3f\\color{black}   {\\bfTIC} %s%.3f', ...
            c1, R2(med(:,j),y1), c1, TIC(med(:,j),y1));
        s2 = sprintf('{\\bfcom aero:}  {\\bfR^2} %s%.3f\\color{black}   {\\bfTIC} %s%.3f', ...
            c2, R2(med(:,j),y2), c2, TIC(med(:,j),y2));
        text(0.008, 0.98, s1, 'Units','normalized','VerticalAlignment','top','FontSize',9);
        text(0.40,  0.98, s2, 'Units','normalized','VerticalAlignment','top','FontSize',9);
        if j == 1
            legend([{'medido'}, NOM], 'Location','northoutside','Orientation','horizontal','FontSize',9);
        end
        if j == 6, xlabel('tempo [s]'); end
    end
    title(tl, sprintf(['Absorção dos efeitos aerodinâmicos: identificação sem e com os termos da estrutura\n' ...
        'validação em malha aberta, trecho %.0f a %.0f s'], tw), 'FontWeight','bold');
    exportgraphics(f, fullfile(DIR, sprintf('valid_%04.0f-%04.0f.png', tw)), 'BackgroundColor','white','Resolution',130);
    close(f);
    fprintf('  %4.0f-%-4.0f ok\n', tw);
end
rmappdata(0,'damp_form_override');
fprintf('\n  Figuras em %s\n', DIR);
