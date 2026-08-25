% show_windows.m — Mostra as janelas de TREINO e VALIDAÇÃO usadas no FEM
% =========================================================================
% Reaproveita as janelas do 3_identification/identify_plant.m. Imprime stats
% de excitação por janela (desvio-padrão de p,q,r — proxy de "quanto excita")
% e plota p,q,r ao longo do log inteiro com as janelas destacadas:
%   verde  = treino   |   vermelho = validação
%
% Uso:  >> show_windows

clear; clc; close all;
here = fileparts(mfilename('fullpath'));                 % .../9_FEM_identification
addpath(fileparts(here));                                % raiz da identification
paths = setup_paths();
img_dir = fullfile(here, 'outputs', 'images');          % saída LOCAL do FEM
if ~exist(img_dir,'dir'), mkdir(img_dir); end

%% ---- CONFIG (mesmas janelas do identify_plant) ----
LOG_FILE = 'logs_concat.mat';
t_trains = {[4,24]; [25,41]; [42,62]; [63,99]; [100,125]};
t_val    = [605, 625];
dt       = 0.1;  SG_ORDER = 2;  SG_FRAME = 7;

%% ---- carregar + filtrar ----
L = load_log_data(fullfile(paths.data, LOG_FILE));
t0 = max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]);
t1 = min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]);
tg = (t0:dt:t1)';
ip = @(t,x) interp1(t,x,tg,'linear');
sg = @(x) sgolayfilt(x,SG_ORDER,SG_FRAME);
P = sg(ip(L.time_IMU,L.gyrX_raw));  Q = sg(ip(L.time_IMU,L.gyrY_raw));  R = sg(ip(L.time_IMU,L.gyrZ_raw));

%% ---- print das janelas ----
fprintf('========================================================================\n');
fprintf('  JANELAS DE IDENTIFICAÇÃO — log "%s" [%.1f, %.1f]s, dt=%.2fs\n', LOG_FILE, t0, t1, dt);
fprintf('========================================================================\n');
fprintf('  %-10s %-12s %7s %8s | std(p)  std(q)  std(r)  [rad/s]\n','tipo','janela','dur(s)','pontos');
fprintf('  %s\n', repmat('-',1,72));
report = @(nm,tw) report_win(nm,tw,tg,P,Q,R);
ttot = 0;
for s=1:numel(t_trains), ttot = ttot + report('treino', t_trains{s}); end
report('valid', t_val);
fprintf('  %s\n', repmat('-',1,72));
fprintf('  Total de treino: %.0f s em %d janelas | validação: %.0f s\n', ttot, numel(t_trains), diff(t_val));

%% ---- plot ----
fig = figure('Name','show_windows','Position',[50 50 1200 720]);
sig = {P,'p (rad/s)'; Q,'q (rad/s)'; R,'r (rad/s)'};
for i=1:3
    subplot(3,1,i); hold on; grid on;
    % faixas das janelas (atrás do sinal)
    yl = [min(sig{i,1}) max(sig{i,1})];
    for s=1:numel(t_trains)
        patch([t_trains{s}(1) t_trains{s}(2) t_trains{s}(2) t_trains{s}(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [0.6 1 0.6], 'EdgeColor','none','FaceAlpha',0.35);
    end
    patch([t_val(1) t_val(2) t_val(2) t_val(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [1 0.6 0.6], 'EdgeColor','none','FaceAlpha',0.35);
    plot(tg, sig{i,1}, 'b', 'LineWidth', 0.8);
    ylabel(sig{i,2});
    if i==1, title('Janelas: verde = treino | vermelho = validação'); end
    if i==3, xlabel('t (s)'); end
end
saveas(fig, fullfile(img_dir,'fem_show_windows.png'));
fprintf('  Figura: %s\n', fullfile(img_dir,'fem_show_windows.png'));

%% ---- helper ----
function dur = report_win(nm, tw, tg, P, Q, R)
    m = (tg>=tw(1)) & (tg<=tw(2));  dur = tw(2)-tw(1);
    fprintf('  %-10s [%4.0f,%4.0f]s %7.0f %8d | %6.3f  %6.3f  %6.3f\n', ...
        nm, tw(1), tw(2), dur, nnz(m), std(P(m)), std(Q(m)), std(R(m)));
end
