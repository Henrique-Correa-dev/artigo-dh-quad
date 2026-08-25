% compare_damping_windows.m — constante × ∝V em janelas de 20 s de 400 a 630 s
% =========================================================================
% Só simula (modo full, mesma cadeia de motor). Janelas que cruzam o intervalo
% entre logs (512,1–517,1 s) são puladas. Reporta R² de p, q, r por janela e
% um gráfico de barras.
%
% Uso:  >> compare_damping_windows
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();  MODE = 'full';
% Trechos com voo efetivo: 400–480 (voo 6) e 525–650 (voo 7).
% O intervalo 480–525 é descartado: contém o fim do voo 6, o intervalo entre
% logs (512,1–517,1 s) e o início do voo 7 (aeronave parada / no solo).
WINS = [(400:20:460)', (420:20:480)'; (525:20:625)', (545:20:645)'; 630 650];

L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';
constants = struct('m', proj.m, 'g', proj.g);
VEst = estimate_velocity(L, tg);  Vgrid = VEst.V;  Vgrid(~isfinite(Vgrid)) = 0;

C = struct('tag',{'unified','damp_aero'}, 'mode',{'const','aero'}, 'c0',{[0;0;0],[0;0;0]}, ...
           'label',{'constante (c_p,c_q,c_r)','\propto V (C_{lp},C_{mq},C_{nr})'}, 'cor',{[0 0.45 0.7],[0.85 0.37 0.01]});
for k = 1:numel(C), C(k).P = load(fullfile(paths.outputs,'diag',C(k).tag,'P_identified.mat')).P_final(:); end
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

gapA = L.boundaries(end);  gapB = L.log_starts(end);        % 512,1 – 517,1
res = nan(size(WINS,1), numel(C), 3);  Vp95 = nan(size(WINS,1),1);  SIM = cell(size(WINS,1),1);
fprintf('\n%-9s | %6s | %-26s | %-26s\n', 'janela', 'V p95', 'constante  R² p / q / r', '∝V  R² p / q / r');
for w = 1:size(WINS,1)
    tw = WINS(w,:);
    if tw(1) < gapB && tw(2) > gapA, fprintf('%3d–%-3d   |  (cruza o intervalo entre logs — pulada)\n', tw); continue; end
    idx = tg>=tw(1) & tg<=tw(2);  time = tg(idx);
    ip = @(tt,xx) interp1(tt, xx, time, 'linear');
    pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
    pqr_meas = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
    att_meas = [ip(L.time_ATT,L.roll_deg), ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
    Vp95(w) = prctile(Vgrid(idx), 95);
    for k = 1:numel(C)
        if strcmp(C(k).mode,'const'), if isappdata(0,'diag_damp'), rmappdata(0,'diag_damp'); end
        else, setappdata(0,'diag_damp', struct('mode',C(k).mode,'t',tg,'V',Vgrid,'c0',C(k).c0)); end
        try
            r = sim_window(MODE, C(k).P, time, pwm, pqr_meas, att_meas, constants);
            res(w,k,:) = [R2(pqr_meas(:,1),r.p), R2(pqr_meas(:,2),r.q), R2(pqr_meas(:,3),r.r)];
            SIM{w}.time = time; SIM{w}.meas = pqr_meas; SIM{w}.V = Vgrid(idx); SIM{w}.sim{k} = [r.p, r.q, r.r];
        catch ME
            fprintf('  (janela %d, %s: %s)\n', w, C(k).tag, ME.message);
        end
    end
    fprintf('%3d–%-3d   | %6.2f | %7.3f %7.3f %7.3f     | %7.3f %7.3f %7.3f\n', tw, Vp95(w), squeeze(res(w,1,:)), squeeze(res(w,2,:)));
end
if isappdata(0,'diag_damp'), rmappdata(0,'diag_damp'); end
ok = all(isfinite(res(:,1,1)),2);
fprintf('%-9s | %6s | %7.3f %7.3f %7.3f     | %7.3f %7.3f %7.3f\n', 'MÉDIA', '', mean(res(ok,1,:),1), mean(res(ok,2,:),1));

%% figura: R² por janela
f = figure('Position',[60 60 1150 640],'Color','w'); try, f.Theme='light'; catch, end
lab = {'p','q','r'};  xt = mean(WINS,2);
for i = 1:3
    subplot(3,1,i); hold on; grid on;
    bar(xt, squeeze(res(:,:,i)), 0.8);
    colororder(vertcat(C.cor));
    ylim([-0.2 1]); ylabel(['R² ' lab{i}]); xlim([395 635]);
    if i == 1, legend({C.label}, 'Location','northoutside','Orientation','horizontal'); end
    if i == 3, xlabel('centro da janela de 20 s [s]'); end
end
exportgraphics(f, fullfile(paths.images,'compare_damping_windows.png'), 'BackgroundColor','white','Resolution',150);
fprintf('  Figura: %s\n', fullfile(paths.images,'compare_damping_windows.png'));

%% figuras de séries temporais por janela (2 painéis: voo 6 e voo 7)
grupos = {find(WINS(:,1) < 512), find(WINS(:,1) > 512)};  nomes = {'400–480 s (voo 6)', '525–650 s (voo 7)'};
for gI = 1:2
    ws = grupos{gI};  ws = ws(cellfun(@(c) ~isempty(c), SIM(ws)));  nR = numel(ws);
    f2 = figure('Position',[40 40 1250 230*nR+60],'Color','w'); try, f2.Theme='light'; catch, end
    tl = tiledlayout(nR, 3, 'TileSpacing','compact', 'Padding','compact');
    for ri = 1:nR
        w = ws(ri);  Sw = SIM{w};
        for i = 1:3
            nexttile; hold on; grid on;
            plot(Sw.time, Sw.meas(:,i), 'k-', 'LineWidth', 1.6);
            for k = 1:numel(C), plot(Sw.time, Sw.sim{k}(:,i), '-', 'Color', C(k).cor, 'LineWidth', 1.2); end
            xlim(WINS(w,:));
            text(0.01, 0.95, sprintf('%s   R^2 const %.2f | \\propto V %.2f   (V_{p95} %.1f m/s)', lab{i}, res(w,1,i), res(w,2,i), Vp95(w)), ...
                'Units','normalized', 'VerticalAlignment','top', 'FontSize', 8.5, 'BackgroundColor',[1 1 1 0.7]);
            if i == 1, ylabel(sprintf('%d–%d s', WINS(w,1), WINS(w,2)), 'FontWeight','bold'); end
            if ri == 1 && i == 2, legend({'medido', C.label}, 'Orientation','horizontal', 'Location','northoutside'); end
            if ri == nR, xlabel('t [s]'); end
        end
    end
    title(tl, sprintf('Constante × \\propto V — janelas de 20 s, %s   (linhas: p, q, r em rad/s)', nomes{gI}));
    fn = fullfile(paths.images, sprintf('compare_damping_windows_series_%d.png', gI));
    exportgraphics(f2, fn, 'BackgroundColor','white', 'Resolution', 130);
    fprintf('  Figura: %s\n', fn);
end
