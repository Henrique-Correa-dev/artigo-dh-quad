% aero_influence.m — Velocidade estimada (barômetro/EKF + atitude) e influência das
%                    forças e momentos aerodinâmicos frente ao empuxo e aos momentos do modelo
% =========================================================================
% PROBLEMA
%   Nas janelas de validação (voo 7) não há GPS com fix. Para avaliar quanto as
%   forças/momentos aerodinâmicos desprezados valem lá, monta-se uma velocidade
%   ESTIMADA com o que há a bordo:
%     • vertical  w  : velocidade Down do EKF do ArduPilot (XKF1.VD), que sem GPS
%                      é auxiliada só pelo barômetro → válida;
%     • horizontal   : integração da aceleração horizontal quase-estática obtida
%                      da ATITUDE medida e do empuxo da cadeia de motor,
%                        a_NED = R_n←b·[0;0;−T/m] + [0;0;g],
%                      com deriva removida ancorando v = 0 no início e no fim de
%                      cada janela (a aeronave está em pairado nas bordas).
%   O método é AFERIDO no voo 4 (que tem GPS) antes de ser aplicado ao voo 7.
%
% O QUE SE COMPARA
%   Forças (N, em % do empuxo T da cadeia de motor):
%     placa plana vertical  ½ρ w² S·C_N (C_N=1,5)   | sustentação da asa ½ρ V_h² S·CL0 (0,36)
%     arrasto horizontal    ½ρ V_h² S·C_D0 (0,08)   | força H dos rotores 4·k_h·V_h
%   Momentos (N·m, em % do momento de amortecimento do modelo cp·Jx·p etc. e do
%   momento de controle |Mx|,|My|,|Mz| da cadeia de motor):
%     amortecimento da estrutura (AVL): L = ρ V_h S b² |Cl_p| p/4, M = ρ V_h S c̄² |Cm_q| q/2,
%                                       N = ρ V_h S b² |Cn_r| r/4
%     estático (asa+cauda):   q̄ S c̄ (|Cm0| + |Cm_α|·|α|), q̄ S b |Cl_β|·|β|, q̄ S b |Cn_β|·|β|
%                             com α, β limitados a 15° (fora disso a asa está estolada e
%                             a força já está no termo de placa plana)
%
% Uso:  >> aero_influence
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();  g = proj.g;  m = proj.m;
rho = 1.225;  S = proj.wing.S;  b = proj.wing.b;  cbar = proj.wing.c;
CL0 = 0.36;  CD0 = 0.08;  CN_plate = 1.5;
Clp = 0.406;  Cmq = 8.96;  Cnr = 0.070;  Cm0 = 0.08;  Cma = 0.74;  Clb = 0.062;  Cnb = 0.066;   % AVL (módulos)
try, pr = load(fullfile(paths.outputs,'prior_damping.mat')).prior; k_h = pr.k_h; catch, k_h = 0.0156; end
prog = load(fullfile(paths.outputs,'P_identified.mat'));  P = prog.P_final(:);
Jx = P(1); Jy = P(2); Jz = P(3); cp = P(13); cq = P(14); cr = P(15);

L  = load_log_data(fullfile(paths.data,'logs_concat.mat'));  dt = 0.1;
tg = (max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]):dt:min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]))';
ip = @(tt,xx) interp1(tt,xx,tg,'linear');
p = ip(L.time_IMU,L.gyrX_raw); q = ip(L.time_IMU,L.gyrY_raw); r = ip(L.time_IMU,L.gyrZ_raw);
phi = deg2rad(ip(L.time_ATT,L.roll_deg)); th = deg2rad(ip(L.time_ATT,L.pitch_deg)); psi = deg2rad(ip(L.time_ATT,L.yaw_deg));
W4 = [ip(L.time_RCOU,L.pwm1_raw) ip(L.time_RCOU,L.pwm2_raw) ip(L.time_RCOU,L.pwm3_raw) ip(L.time_RCOU,L.pwm4_raw)];
[rpm, fT, fQ] = motor_chain(tg, W4);
Tmr = fT(rpm).*P(5:8)';  Qmr = fQ(rpm).*P(9:12)';  T = sum(Tmr,2);
dyn = vtol_dynamics('get_handles');
[Mx, My, Mz] = dyn.moments(Tmr, Qmr, proj.arms.Lx_r, proj.arms.Lx_l, proj.arms.Ly_f, proj.arms.Ly_r);

%% ---------- EKF VD (barômetro) dos logs brutos → eixo do concat ----------
VD = nan(size(tg));
for k = 1:numel(L.log_names)
    raw = fullfile(paths.data, L.log_names{k});
    try
        Lr = load_log_data(raw);
        t0_log = min([Lr.time_IMU(1), Lr.time_ATT(1), Lr.time_RCOU(1)]);
        t_shift = L.log_starts(k) - t0_log;
        X = load(raw, 'XKF1_0').XKF1_0;             % cols: 2 TimeUS, 7 VN, 8 VE, 9 VD
        tx = X(:,2)/1e6 + t_shift;
        [tx, iu] = unique(tx);
        vd_k = interp1(tx, X(iu,9), tg, 'linear', NaN);
        sel = tg >= L.log_starts(k) & (k==numel(L.log_names) | tg <= L.boundaries(min(k,end)));
        VD(sel) = vd_k(sel);
    catch ME
        fprintf('  (XKF1 do log %d indisponível: %s)\n', k, ME.message);
    end
end

%% ---------- Velocidade horizontal por atitude (integração ancorada) ----------
% aceleração NED do modelo quase-estático: R_n←b·[0;0;−T/m] + g·e3
aN = zeros(size(tg)); aE = aN;
for k = 1:numel(tg)
    cph=cos(phi(k)); sph=sin(phi(k)); cth=cos(th(k)); sth=sin(th(k)); cps=cos(psi(k)); sps=sin(psi(k));
    % terceira coluna de R_n←b (direção de −z_b em NED)
    e3n = [cph*sth*cps + sph*sps;  cph*sth*sps - sph*cps;  cph*cth];
    a = -(T(k)/m)*e3n + [0;0;g];
    aN(k) = a(1); aE(k) = a(2);
end
% integra com ancoragem v=0 nas bordas de SUB-JANELAS de ANCHOR_S s (a aeronave volta
% ao pairado entre doublets; a ancoragem curta remove a deriva por viés de atitude,
% que em janelas longas cresce de forma quadrática). Aferido no voo 4 contra o GPS.
ANCHOR_S = 10;
vint = @(a, mm) local_anchored_int(a, tg, mm, dt, ANCHOR_S);

%% ---------- Janelas ----------
WIN = struct('nome',{'treino (voo 4)','validação 610–630','validação 550–600'}, ...
             'ints',{[4 24; 25 41; 42 62; 63 99; 100 125], [610 630], [550 600]});
fprintf('\n=== aero_influence: velocidades estimadas e influência aerodinâmica ===\n');
fprintf('  Cadeia: T de motor_chain+P_final | AVL: Cl_p %.3f Cm_q %.2f Cn_r %.3f | cp cq cr = %.2f %.2f %.2f\n', Clp, Cmq, Cnr, cp, cq, cr);

% aferição no voo 4: horizontal estimada × GPS
haveGPS = ~isempty(L.time_GPS);
if haveGPS
    tG = L.time_GPS(:); VNg = interp1(tG, L.gps_vn(:), tg, 'linear', NaN); VEg = interp1(tG, L.gps_ve(:), tg, 'linear', NaN);
end

RES = struct();
for wI = 1:numel(WIN)
    ints = WIN(wI).ints;  mm = false(size(tg));
    VN_e = nan(size(tg)); VE_e = VN_e;
    for i = 1:size(ints,1)
        mi = tg>=ints(i,1) & tg<=ints(i,2);  mm = mm | mi;
        VN_e(mi) = vint(aN, mi);  VE_e(mi) = vint(aE, mi);
    end
    Vh = hypot(VN_e, VE_e);
    % corpo
    u = nan(size(tg)); v = u; w = u;
    for k = find(mm)'
        cph=cos(phi(k)); sph=sin(phi(k)); cth=cos(th(k)); sth=sin(th(k)); cps=cos(psi(k)); sps=sin(psi(k));
        Rbn = [ cth*cps, cth*sps, -sth; sph*sth*cps-cph*sps, sph*sth*sps+cph*cps, sph*cth; cph*sth*cps+sph*sps, cph*sth*sps-sph*cps, cph*cth ];
        vb = Rbn*[VN_e(k); VE_e(k); VD(k)];  u(k)=vb(1); v(k)=vb(2); w(k)=vb(3);
    end
    st = @(x) [prctile(abs(x(mm & isfinite(x))),95), max(abs(x(mm & isfinite(x))))];
    sVh = st(Vh);  sw = st(w);
    fprintf('\n  --- %s (%d s) ---\n', WIN(wI).nome, nnz(mm)*dt);
    fprintf('  |V_h| estimada: p95 %.2f | máx %.2f m/s    |w| EKF/baro: p95 %.2f | máx %.2f m/s\n', sVh, sw);
    if haveGPS && all(isfinite(VNg(mm)))
        Vg = hypot(VNg, VEg);
        e = Vh(mm) - Vg(mm);
        fprintf('  AFERIÇÃO × GPS: |V_h| GPS p95 %.2f máx %.2f | corr(est,GPS) %.2f | RMSE %.2f m/s | viés %+.2f m/s\n', ...
            prctile(Vg(mm),95), max(Vg(mm)), corr(Vh(mm),Vg(mm)), rms(e), mean(e));
    end
    % ---- forças (N) vs empuxo ----
    Fplate = 0.5*rho*w.^2*S*CN_plate;   Flift = 0.5*rho*Vh.^2*S*CL0;   Fdrag = 0.5*rho*Vh.^2*S*CD0;   FH = 4*k_h*Vh;
    Ftot = hypot(Fplate, hypot(Flift,Fdrag)+FH);       % soma conservadora
    rF = 100*Ftot./T;
    fprintf('  Forças aero / empuxo T:  p95 %.2f%% | máx %.2f%%   (placa vertical p95 %.2f%%, asa %.2f%%, arrasto %.2f%%, força H %.2f%%)\n', ...
        prctile(rF(mm),95), max(rF(mm)), prctile(100*Fplate(mm)./T(mm),95), prctile(100*Flift(mm)./T(mm),95), prctile(100*Fdrag(mm)./T(mm),95), prctile(100*FH(mm)./T(mm),95));
    % ---- momentos de amortecimento aero vs modelo ----
    Ld_a = rho*Vh*S*b^2*Clp/4.*abs(p);   Md_a = rho*Vh*S*cbar^2*Cmq/2.*abs(q);   Nd_a = rho*Vh*S*b^2*Cnr/4.*abs(r);
    Ld_m = cp*Jx*abs(p);  Md_m = cq*Jy*abs(q);  Nd_m = cr*Jz*abs(r);
    act = mm & abs(p)>0.2;  actq = mm & abs(q)>0.2;  actr = mm & abs(r)>0.1;
    fprintf('  Amortecimento aero (AVL) / amortecimento do modelo (razão dos p95 em manobra): rol %.0f%% | arf %.0f%% | gui %.0f%%\n', ...
        100*prctile(Ld_a(act),95)/prctile(Ld_m(act),95), 100*prctile(Md_a(actq),95)/prctile(Md_m(actq),95), 100*prctile(Nd_a(actr),95)/prctile(Nd_m(actr),95));
    % ---- momentos estáticos aero vs momento de controle ----
    Vtot = hypot(Vh, w);  qbar = 0.5*rho*Vtot.^2;
    alpha = min(abs(atan2(w, max(abs(u),0.3))), deg2rad(15));  beta = min(abs(atan2(v, max(Vtot,0.3))), deg2rad(15));
    Ms_a = qbar*S*cbar.*(Cm0 + Cma*alpha);   Ls_a = qbar*S*b*Clb.*beta;   Ns_a = qbar*S*b*Cnb.*beta;
    Mc = abs(My); Lc = abs(Mx); Nc = abs(Mz);
    fprintf('  Momento estático aero / momento de controle (razão dos p95): rol %.1f%% | arf %.1f%% | gui %.1f%%\n', ...
        100*prctile(Ls_a(mm),95)/prctile(Lc(mm),95), 100*prctile(Ms_a(mm),95)/prctile(Mc(mm),95), 100*prctile(Ns_a(mm),95)/prctile(Nc(mm),95));
    fprintf('  Momento estático aero absoluto: p95 rol %.4f | arf %.4f | gui %.4f N·m   (controle p95: %.3f | %.3f | %.3f N·m)\n', ...
        prctile(Ls_a(mm),95), prctile(Ms_a(mm),95), prctile(Ns_a(mm),95), prctile(Lc(mm),95), prctile(Mc(mm),95), prctile(Nc(mm),95));
    RES(wI).mm = mm; RES(wI).Vh = Vh; RES(wI).w = w; RES(wI).rF = rF; RES(wI).Faero = Ftot; RES(wI).T = T;
    RES(wI).Ld_a = Ld_a; RES(wI).Ld_m = Ld_m; RES(wI).Md_a = Md_a; RES(wI).Md_m = Md_m; RES(wI).Ms_a = Ms_a; RES(wI).Mc = Mc;
end

%% ---------- FIGURA ----------
f = figure('Position',[40 40 1150 1000],'Color','w'); try, f.Theme='light'; catch, end
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
% (a) aferição no voo 4 (janela 63–99)
nexttile; hold on; grid on; sel = tg>=63 & tg<=99;
if haveGPS, plot(tg(sel), hypot(VNg(sel),VEg(sel)), 'k-', 'LineWidth',1.6); end
plot(tg(sel), RES(1).Vh(sel), '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.4);
plot(tg(sel), abs(RES(1).w(sel)), '-', 'Color',[0 0.45 0.7], 'LineWidth',1.2);
ylabel('[m/s]'); legend({'|V_h| GPS','|V_h| estimada (atitude+empuxo, ancorada 10 s)','|w| EKF/barômetro'},'Location','northwest');
text(0.005,1.06,'(a) Aferição no voo 4 (treino, 63–99 s): velocidade horizontal estimada × GPS','Units','normalized','FontWeight','bold');
% (b) validação 610–630: velocidades
nexttile; hold on; grid on; sel = RES(2).mm;
plot(tg(sel), RES(2).Vh(sel), '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.4);
plot(tg(sel), abs(RES(2).w(sel)), '-', 'Color',[0 0.45 0.7], 'LineWidth',1.2);
ylabel('[m/s]'); legend({'|V_h| estimada','|w| EKF/barômetro'},'Location','northwest');
text(0.005,1.06,'(b) Validação 610–630 s (voo 7, sem GPS): velocidades estimadas','Units','normalized','FontWeight','bold');
% (c) forças aero (N) vs empuxo
nexttile; hold on; grid on;
plot(tg(sel), RES(2).Faero(sel), 'k-', 'LineWidth',1.4);
plot(tg(sel), 0.05*RES(2).T(sel), '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
plot(tg(sel), 0.01*RES(2).T(sel), ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylabel('[N]'); legend({'forças aerodinâmicas (soma conservadora)','5% do empuxo','1% do empuxo'},'Location','northwest');
text(0.005,1.06,'(c) Forças aerodinâmicas não modeladas × empuxo, na validação','Units','normalized','FontWeight','bold');
% (d) momentos de arfagem (N·m): controle, amortecimento do modelo, amortecimento aero, estático aero
nexttile; hold on; grid on;
plot(tg(sel), RES(2).Mc(sel), 'k-', 'LineWidth',1.4);
plot(tg(sel), RES(2).Md_m(sel), '-', 'Color',[0.3 0.6 0.3], 'LineWidth',1.3);
plot(tg(sel), RES(2).Md_a(sel), '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.3);
plot(tg(sel), RES(2).Ms_a(sel), '-', 'Color',[0 0.45 0.7], 'LineWidth',1.3);
ylabel('[N·m]'); xlabel('tempo [s]');
legend({'|M_y| de controle (cadeia de motor)','amortecimento do modelo c_q·J_y·|q|','amortecimento da estrutura (AVL C_{mq})','estático da estrutura (C_{m0}, C_{m\alpha})'},'Location','northwest');
text(0.005,1.06,'(d) Momentos de arfagem na validação: o que o modelo usa × o que a aerodinâmica da estrutura daria','Units','normalized','FontWeight','bold');
exportgraphics(f, fullfile(paths.images,'aero_influence.png'), 'BackgroundColor','white','Resolution',150);
fprintf('\n  Figura: %s\n', fullfile(paths.images,'aero_influence.png'));

%% ---------------- helper ----------------
function v = local_anchored_int(a, t, mm, dt, anchor_s)
%LOCAL_ANCHORED_INT integra a em mm, com v=0 nas bordas de sub-janelas de anchor_s s
    idx = find(mm);  a_w = a(idx);  a_w(~isfinite(a_w)) = 0;  t_w = t(idx);
    v = zeros(size(a_w));  n_sub = max(1, round(anchor_s/dt));
    for s0 = 1:n_sub:numel(a_w)
        s1 = min(s0+n_sub-1, numel(a_w));  ii = s0:s1;
        vv = cumtrapz(t_w(ii), a_w(ii));  n = numel(vv);
        if n > 1, vv = vv - vv(1) - (vv(end)-vv(1))*((0:n-1)'/(n-1)); end
        v(ii) = vv;
    end
end
