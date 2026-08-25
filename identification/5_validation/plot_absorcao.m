% plot_absorcao.m — comparação da absorção: modelo com e sem termos aerodinâmicos
% =========================================================================
% Duas identificações completas do mesmo conjunto de dados:
%   SEM aero   moment_win6_2026     L_p, M_q, N_r absorvem tudo (P(17:22) = 0)
%   COM aero   oficial_2026_final   asa medida em voo fixa + F_aero, mesmos dados
% Simula as duas em todas as janelas de validação, malha aberta, modo full, e
% salva uma figura de seis painéis por janela em duas pastas:
%   outputs/images/absorcao/sem_aero/   e   outputs/images/absorcao/com_aero/
% Uso:  >> plot_absorcao
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
WINS = [605 625; 610 630; 550 600; 400 420; 420 440; 440 460; ...
        460 480; 525 545; 545 565; 565 585; 585 605; 630 650];
RUNS = {'sem_aero', 'moment_win6_2026',  '', 'SEM termos aerodinâmicos (L_p, M_q, N_r absorvem a estrutura)'; ...
        'com_aero', 'oficial_2026_final','measured', 'COM termos aerodinâmicos (asa medida em voo, F_aero medida)'};

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
RES = nan(size(WINS,1), 2, 6);

for m = 1:2
    S = load(fullfile(paths.outputs,'runs',RUNS{m,2},'P_identified.mat'));
    P = S.P_final(:);
    DIR = fullfile(paths.images,'absorcao',RUNS{m,1});
    if ~exist(DIR,'dir'), mkdir(DIR); end
    if isempty(RUNS{m,3})
        if isappdata(0,'faero_on'), rmappdata(0,'faero_on'); end
    else
        setappdata(0,'faero_on',RUNS{m,3});
    end
    fprintf('\n===== %s (%s): L_p %.4f  M_q %.4f  N_r %.4f\n', RUNS{m,1}, RUNS{m,2}, P(13:15));
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
        r1 = sim_window('full', P, time, pwm, pqr, att, constants);
        f = figure('Position',[30 30 1150 1250],'Color','w','Visible','off');
        try, f.Theme='light'; catch, end
        tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
        for j = 1:6
            nexttile; hold on; grid on;
            y1 = r1.(canais{j});
            plot(time, med(:,j), 'k-', 'LineWidth', 1.7);
            plot(time, y1, '-', 'Color',[0.85 0.37 0.01], 'LineWidth', 1.3);
            RES(w,m,j) = R2(med(:,j), y1);
            ylabel(rot{j}); xlim(tw);
            text(0.006, 0.94, sprintf('R² %.3f   TIC %.3f', RES(w,m,j), TIC(med(:,j),y1)), ...
                'Units','normalized','VerticalAlignment','top','FontSize',9, ...
                'FontName','Menlo','BackgroundColor',[1 1 1 0.78]);
            if j == 1, legend({'medido','modelo'},'Location','northoutside','Orientation','horizontal','FontSize',9); end
            if j == 6, xlabel('tempo [s]'); end
        end
        title(tl, sprintf('%s\nvalidação em malha aberta, trecho %.0f a %.0f s', RUNS{m,4}, tw), 'FontWeight','bold');
        exportgraphics(f, fullfile(DIR, sprintf('valid_%04.0f-%04.0f.png', tw)), 'BackgroundColor','white','Resolution',130);
        close(f);
        fprintf('  %4.0f–%-4.0f  R² %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f\n', tw, squeeze(RES(w,m,:)));
    end
end
rmappdata(0,'damp_form_override');
if isappdata(0,'faero_on'), rmappdata(0,'faero_on'); end

fprintf('\n===== MÉDIA das janelas válidas\n');
fprintf('  %-10s %6s %6s %6s | %6s %6s %6s\n','modelo','p','q','r','a_x','a_y','a_z');
for m = 1:2
    mu = squeeze(mean(RES(:,m,:),1,'omitnan'));
    fprintf('  %-10s %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f\n', RUNS{m,1}, mu);
end
fprintf('\n  Figuras em %s\n', fullfile(paths.images,'absorcao'));
