% results_id_figures.m — Figuras da identificação para a dissertação (§5.2)
% =========================================================================
% Gera 6 figuras (fundo branco), salvas em outputs/images/ e copiadas p/ Cap5:
%   (1) id_sim_p0     — simulação do modelo inicial P0 vs medições (p,q,r,a)
%   (2) id_cost_eem   — custo do EEM por iteração (decrescente)
%   (3) id_cost_oem   — custo do OEM por iteração (estágios 1s/2s/3s)
%   (4) id_sim_comp   — sobreposição P0 / EEM / OEM vs medições (grade 2x3)
%   (5) id_comp_rates — idem, estilo empilhado por canal (p,q,r) com R^2 no título
%   (6) id_comp_acc   — idem, estilo empilhado por canal (a_x,a_y,a_z)
%
% Simulação: integra as taxas (p,q,r) com o MODELO DE MOTOR FÍSICO (= identify_plant),
% PWM com atraso+lag idêntico ao estimador, e o acelerômetro pelo modelo de sensor
% (não depende de u,v,w). Atitude medida (modo híbrido — sem deriva).
% Históricos de custo: da rodada do identify_plant (arrays em CFG).
% Uso: >> results_id_figures

clear; clc; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
paths = setup_paths();

%% ---------------- CONFIG ----------------
LOG_FILE = 'logs_concat.mat';
dt = 0.1;  DELAY_PWM = 1;
T_VAL = [610 630];   % janela de validacao oficial (= identify_plant)
CAP5  = '/Users/graest/ita-master/artigo/artigo-dh-quad/dissertacao/Cap5';

proj = parameters();
tau_lag   = proj.motor.tau_m;
constants = struct('m', proj.m, 'g', proj.g);

%% ---------------- VETORES DE PARÂMETROS (rodada oficial_2026_teo) ----------------
SS = load(fullfile(paths.outputs,'runs','oficial_2026_teo','summary.mat'));
P0    = SS.summary.P0(:);
P_eem = SS.summary.P_eem(:);
P_oem = SS.summary.P_final(:);

% velocidade estimada p/ os termos aerodinamicos fixos (mesma serie do estimador)
Lv = load_log_data(fullfile(paths.data, LOG_FILE));
t0v = max([min(Lv.time_IMU), min(Lv.time_ATT), min(Lv.time_RCOU)]);
t1v = min([max(Lv.time_IMU), max(Lv.time_ATT), max(Lv.time_RCOU)]);
tcv = (t0v:dt:t1v)';
VEv = estimate_velocity(Lv, tcv);
Vg=VEv.V; vg=VEv.v; wg=VEv.w; Vg(~isfinite(Vg))=0; vg(~isfinite(vg))=0; wg(~isfinite(wg))=0;
setappdata(0,'aero_vel',struct('t',tcv,'V',Vg,'v',vg,'w',wg));

%% ---------------- HISTÓRICOS DE CUSTO (rodada identify_plant) ----------------
eem_cost = [6162.55 2245.99 918.701 661.953 607.713 599.954 599.5 599.495 ...
            599.493 599.493 599.492*ones(1,16)];
oem_1s = [79265.9 5601.13 2696.44 2577.25 2576.72 2576.71*ones(1,7)];
oem_2s = [2634.95 2632.28 2632.26*ones(1,6)];
oem_3s = [2658.47 2656.85 2656.84*ones(1,6)];

%% ---------------- CARREGAR + JANELA DE VALIDAÇÃO ----------------
L  = load_log_data(fullfile(paths.data, LOG_FILE));
t0 = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t1 = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t0:dt:t1)';  ip = @(tt,xx) interp1(tt,xx,tg,'linear');

W = [ip(L.time_RCOU,L.pwm1_raw) ip(L.time_RCOU,L.pwm2_raw) ...
     ip(L.time_RCOU,L.pwm3_raw) ip(L.time_RCOU,L.pwm4_raw)];
% CADEIA DE ATUAÇÃO ÚNICA (2_model/motor_chain.m): atraso + RPM_ss + lag exato.
% W passa a ser RPM e fT/fQ (abaixo) agem sobre RPM — idêntico ao identify_plant.
W = motor_chain(tg, W);

sel = tg>=T_VAL(1) & tg<=T_VAL(2);
time     = tg(sel) - tg(find(sel,1));
pwm      = W(sel,:);
pqr_meas = [ip(L.time_IMU,L.gyrX_raw) ip(L.time_IMU,L.gyrY_raw) ip(L.time_IMU,L.gyrZ_raw)]; pqr_meas = pqr_meas(sel,:);
acc_meas = [ip(L.time_IMU,L.accX_raw) ip(L.time_IMU,L.accY_raw) ip(L.time_IMU,L.accZ_raw)]; acc_meas = acc_meas(sel,:);

%% ---------------- SIMULAR (motor físico, taxas integradas) ----------------
r0 = sim_rates(P0,    time, pwm, pqr_meas, constants, dt);
re = sim_rates(P_eem, time, pwm, pqr_meas, constants, dt);
ro = sim_rates(P_oem, time, pwm, pqr_meas, constants, dt);

R2 = @(z,y) 1 - sum((z-y).^2)/max(sum((z-mean(z)).^2),1e-12);
MEAS = {pqr_meas(:,1),pqr_meas(:,2),pqr_meas(:,3),acc_meas(:,1),acc_meas(:,2),acc_meas(:,3)};
lab  = {'p [rad/s]','q [rad/s]','r [rad/s]','a_x [m/s^2]','a_y [m/s^2]','a_z [m/s^2]'};
sig  = @(R) {R.p,R.q,R.r,R.accX,R.accY,R.accZ};
prnt = @(tag,R) fprintf('  %-4s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n', tag, ...
    R2(MEAS{1},R.p),R2(MEAS{2},R.q),R2(MEAS{3},R.r),R2(MEAS{4},R.accX),R2(MEAS{5},R.accY),R2(MEAS{6},R.accZ));
fprintf('\n  R^2 (validacao 610-630 s):\n         p       q       r      ax      ay      az\n');
prnt('P0',r0); prnt('EEM',re); prnt('OEM',ro);

% --- metricas de validacao por modelo (R2, TIC, decomposicao de Theil) ---
fprintf('\n  ===== METRICAS DE VALIDACAO (R2 | TIC | Theil Ub/Uv/Uc) =====\n');
chn  = {'p','q','r','a_x','a_y','a_z'};
mres = {'P0',r0; 'EEM',re; 'OEM',ro};
for mm = 1:size(mres,1)
    Sm = sig(mres{mm,2});
    fprintf('\n  --- %s ---\n  canal     R2       TIC      Ub     Uv     Uc\n', mres{mm,1});
    for i = 1:6
        M = metr(MEAS{i}, Sm{i});
        fprintf('  %-6s %8.3f %8.3f %6.2f %6.2f %6.2f\n', chn{i}, M.R2, M.TIC, M.Ub, M.Uv, M.Uc);
    end
end

% limites de eixo de aceleração (cliva o az do EEM, que tem k_T~0.5)
ylac = cell(1,3);
for i=4:6, d=MEAS{i}; pad=0.6*(max(d)-min(d))+0.1; ylac{i-3}=[min(d)-pad max(d)+pad]; end

%% ---------------- FIG 1: SIM P0 ----------------
% Estilo p/ dissertação: medido = preto contínuo; Theta_0 = cinza tracejado.
% Sub-legendas (a)–(f) em cada painel; legenda única na base; fontes maiores.
FS1 = 13;
% Paleta única do Cap.5 (medido / Theta_0 / EEM / OEM) — consistente entre figuras
C_MEAS = [0 0 0];          ST_MEAS = '-';    LW_MEAS = 1.3;
C_P0   = [0.45 0.45 0.45]; ST_P0   = '--';   LW_P0   = 1.6;
C_EEM  = [0.00 0.45 0.70]; ST_EEM  = '-.';   LW_EEM  = 1.4;   % azul (Okabe-Ito)
C_OEM  = [0.85 0.37 0.01]; ST_OEM  = '-';    LW_OEM  = 1.7;   % laranja (Okabe-Ito)
f1 = figure('Position',[40 40 1150 760]); set(f1,'Color','w','DefaultAxesFontSize',FS1);
try, f1.Theme='light'; catch, end
s0  = sig(r0);
tl1 = tiledlayout(f1, 2, 3, 'TileSpacing','compact', 'Padding','compact');
sub = {'(a)','(b)','(c)','(d)','(e)','(f)'};
hL1 = gobjects(1,2);
for i=1:6
    ax = nexttile(tl1); hold(ax,'on'); grid(ax,'on');
    h1 = plot(ax, time, MEAS{i}, ST_MEAS, 'Color',C_MEAS, 'LineWidth',LW_MEAS);
    h2 = plot(ax, time, s0{i},   ST_P0,   'Color',C_P0,   'LineWidth',LW_P0);
    if i==1, hL1 = [h1 h2]; end
    ylabel(ax, lab{i});  if i>3, xlabel(ax,'t [s]'); ylim(ax, ylac{i-3}); end
    xlim(ax, [time(1) time(end)]);
    % sub-legenda (a)–(f) + R^2 — texto ancorado no canto sup. esquerdo do painel
    text(ax, 0.02, 0.97, sprintf('%s  R^2 = %.2f', sub{i}, R2(MEAS{i}, s0{i})), ...
         'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'FontSize',FS1, 'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);
end
lg1 = legend(hL1, {'medido','simulado (\Theta_0)'}, 'Orientation','horizontal', 'FontSize',FS1);
lg1.Layout.Tile = 'south';
exportgraphics(f1, fullfile(paths.images,'id_sim_p0.png'), 'BackgroundColor','white', 'Resolution',200);

%% ---------------- FIG 2: CUSTO EEM ----------------
f2 = figure('Position',[60 60 720 460]); set(f2,'Color','w'); try, f2.Theme='light'; catch, end
semilogy(0:numel(eem_cost)-1, eem_cost, '-o', 'LineWidth',1.4, 'MarkerSize',4, 'Color',[0.1 0.45 0.8]);
grid on; xlabel('iteração'); ylabel('custo (soma dos resíduos^2)');
title('Convergência do custo — Fase A (EEM)');
exportgraphics(f2, fullfile(paths.images,'id_cost_eem.png'), 'BackgroundColor','white');

%% ---------------- FIG 3: CUSTO OEM (estágios) ----------------
f3 = figure('Position',[80 80 760 460]); set(f3,'Color','w'); try, f3.Theme='light'; catch, end
oem_all = [oem_1s oem_2s oem_3s];  n1=numel(oem_1s); n2=numel(oem_2s);
semilogy(0:numel(oem_all)-1, oem_all, '-o', 'LineWidth',1.4, 'MarkerSize',4, 'Color',[0.1 0.55 0.3]);
hold on; grid on;
xline(n1-0.5,'--','1s \rightarrow 2s','Color',[.4 .4 .4]);
xline(n1+n2-0.5,'--','2s \rightarrow 3s','Color',[.4 .4 .4]);
xlabel('iteração (estágios 1s, 2s, 3s)'); ylabel('custo (soma dos resíduos^2)');
title('Convergência do custo — Fase B (OEM progressivo)');
exportgraphics(f3, fullfile(paths.images,'id_cost_oem.png'), 'BackgroundColor','white');

%% ---------------- FIG 4: SOBREPOSIÇÃO P0 / EEM / OEM ----------------
f4 = figure('Position',[40 40 1150 760]); set(f4,'Color','w','DefaultAxesFontSize',FS1);
try, f4.Theme='light'; catch, end
se = sig(re); so = sig(ro);
tl4 = tiledlayout(f4, 2, 3, 'TileSpacing','compact', 'Padding','compact');
hL4 = gobjects(1,4);
for i=1:6
    ax = nexttile(tl4); hold(ax,'on'); grid(ax,'on');
    h1 = plot(ax, time, MEAS{i}, ST_MEAS, 'Color',C_MEAS, 'LineWidth',LW_MEAS);
    h2 = plot(ax, time, s0{i},   ST_P0,   'Color',C_P0,   'LineWidth',LW_P0);
    h3 = plot(ax, time, se{i},   ST_EEM,  'Color',C_EEM,  'LineWidth',LW_EEM);
    h4 = plot(ax, time, so{i},   ST_OEM,  'Color',C_OEM,  'LineWidth',LW_OEM);
    if i==1, hL4 = [h1 h2 h3 h4]; end
    ylabel(ax, lab{i});  if i>3, xlabel(ax,'t [s]'); ylim(ax, ylac{i-3}); end
    xlim(ax, [time(1) time(end)]);
    text(ax, 0.02, 0.97, sub{i}, 'Units','normalized', 'HorizontalAlignment','left', ...
         'VerticalAlignment','top', 'FontSize',FS1, 'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);
end
lg4 = legend(hL4, {'medido','\Theta_0','\Theta_{EEM}','\Theta_{OEM}'}, 'Orientation','horizontal', 'FontSize',FS1);
lg4.Layout.Tile = 'south';
exportgraphics(f4, fullfile(paths.images,'id_sim_comp.png'), 'BackgroundColor','white', 'Resolution',200);

%% ---------------- FIG 5: GRADE 3x2 (p,q,r,ax,ay,az) — estilo val_aero ----------------
%  Um unico painel 3x2 com os seis canais, R^2 dos tres modelos em caixa no topo
%  de cada subplot e legenda unica na base (leitura melhor que os empilhados).
f5 = figure('Position',[40 40 1040 900]); set(f5,'Color','w','DefaultAxesFontSize',FS1); try, f5.Theme='light'; catch, end
tl5 = tiledlayout(f5,3,2,'TileSpacing','compact','Padding','compact');
lab6 = {'p [rad/s]','q [rad/s]','r [rad/s]','a_x [m/s^2]','a_y [m/s^2]','a_z [m/s^2]'};
hL5 = [];
for k=1:6
    ax = nexttile(tl5); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    h1=plot(ax, time, MEAS{k}, ST_MEAS, 'Color',C_MEAS, 'LineWidth',1.1);
    h2=plot(ax, time, s0{k},   ST_P0,   'Color',C_P0,   'LineWidth',1.8);
    h3=plot(ax, time, se{k},   ST_EEM,  'Color',C_EEM,  'LineWidth',1.6);
    h4=plot(ax, time, so{k},   ST_OEM,  'Color',C_OEM,  'LineWidth',2.4);
    if k==1, hL5=[h1 h2 h3 h4]; end
    ylabel(ax, lab6{k});  xlim(ax,[time(1) time(end)]);
    if k>=4, ylim(ax, ylac{k-3}); end
    text(ax, 0.03, 1.10, sprintf('R^2:  \\Theta_0 = %.2f   \\Theta_{EEM} = %.2f   {\\color[rgb]{0.85,0.37,0.01}\\Theta_{OEM} = %.2f}', ...
        R2(MEAS{k},s0{k}), R2(MEAS{k},se{k}), R2(MEAS{k},so{k})), ...
        'Units','normalized','VerticalAlignment','top','FontSize',FS1+1,'FontWeight','bold', ...
        'BackgroundColor','w','EdgeColor',[.6 .6 .6],'Margin',3);
    if k>=5, xlabel(ax,'t [s]'); end
end
lg5 = legend(hL5,{'medido','\Theta_0','\Theta_{EEM}','\Theta_{OEM}'},'Orientation','horizontal','FontSize',FS1+2);  lg5.Layout.Tile='south';
exportgraphics(f5, fullfile(paths.images,'id_comp_all.png'), 'BackgroundColor','white', 'Resolution',200);

%% ---------------- COPIAR PARA Cap5 ----------------
for f = {'id_sim_p0','id_cost_eem','id_cost_oem','id_sim_comp','id_comp_all'}
    copyfile(fullfile(paths.images,[f{1} '.png']), fullfile(CAP5,[f{1} '.png']));
end
fprintf('\n  6 figuras salvas e copiadas para Cap5/.\n');


%% ======================= SIMULADOR (motor físico) =======================
function res = sim_rates(P, time, pwm, pqr_meas, constants, dt)
% Integra p,q,r com o modelo de motor FÍSICO (= identify_plant) e calcula o
% acelerômetro pelo modelo de sensor. Atitude medida → sem deriva.
    [fT, fQ] = motor_chain('handles');       % T = kT_rpm·RPM² (mesma cadeia)
    odef = @(t,y) vtol_dynamics(t, y, P, time, pwm, fT, fQ, constants);
    [ts, ys] = ode45(odef, time, pqr_meas(1,:)', odeset('RelTol',1e-6,'AbsTol',1e-9));
    y = interp1(ts, ys, time, 'linear', 'extrap');
    res.p = y(:,1);  res.q = y(:,2);  res.r = y(:,3);

    proj = parameters();  m = constants.m;  r_imu = proj.imu_offset;
    Tmr = (P(5:8)') .* fT(pwm);              % Nx4 empuxos (físico, escalado por k_T)
    T_total = sum(Tmr, 2);
    pd = gradient(res.p, dt);  qd = gradient(res.q, dt);  rd = gradient(res.r, dt);
    [res.accX, res.accY, res.accZ] = accelerometer_model( ...
        res.p, res.q, res.r, zeros(size(res.p)), zeros(size(res.p)), zeros(size(res.p)), ...
        T_total/m, pd, qd, rd, r_imu);
end

function M = metr(z, y)
% R2, TIC (Theil inequality coeff) e decomposicao de Theil (Ub+Uv+Uc=1).
    z = z(:);  y = y(:);  e = z - y;  mse = mean(e.^2) + 1e-12;
    M.R2  = 1 - sum(e.^2) / max(sum((z-mean(z)).^2), 1e-12);
    M.TIC = sqrt(mean(e.^2)) / (sqrt(mean(z.^2)) + sqrt(mean(y.^2)) + 1e-12);
    sz = std(z,1);  sy = std(y,1);  rho = corr(z,y);
    M.Ub = (mean(z)-mean(y))^2 / mse;     % vies (offset sistematico)
    M.Uv = (sz - sy)^2 / mse;             % variancia (amplitude/amortecimento)
    M.Uc = 2*(1-rho)*sz*sy / mse;         % covariancia (erro aleatorio)
end
