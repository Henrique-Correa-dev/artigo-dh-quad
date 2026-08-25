% drag_probe.m — O que sobra em a_x, a_y é arrasto induzido dos rotores? (teste barato)
% =========================================================================
% PERGUNTA
%   Beard & McLain (cap. 14.4.2): a única força além de gravidade e empuxo no
%   multirrotor é o arrasto induzido dos rotores, no plano do corpo,
%       f_d^b ≈ −T·C_d·diag(1,1,0)·v^b   ≈  −m·g·C_d·[u; v; 0]
%   Em força específica: a_x ganha −g·C_d·u e a_y ganha −g·C_d·v.
%   O modelo atual não tem esse termo (parameters().k_drag = 0) e é justamente
%   em a_x, a_y que a validação é fraca. Aqui NÃO se identifica nada: só se
%   verifica se o RESÍDUO do acelerômetro horizontal é proporcional à
%   velocidade do corpo, e qual C_d sairia.
%
% DADOS
%   Só o 1º trecho do logs_concat (0–268 s, log "4 25-05-2026") tem GPS com fix,
%   logo velocidade sobre o solo medida (Spd, GCrs a ~5 Hz). Nos outros dois
%   voos o GPS não fixou e o EKF não teve auxílio (PN = PE = 0 no XKF1), então
%   a velocidade do EKF é integração pura e não serve de referência.
%   Hipótese: vento fraco/constante → velocidade sobre o solo ≈ velocidade
%   aerodinâmica; vento constante vira intercepto na regressão (não o coef.).
%
% MÉTODO
%   1. u, v = R_b←n(φ,θ,ψ)·[VN; VE; VD] com VN = Spd·cos(GCrs), VE = Spd·sin(GCrs)
%   2. resíduo r_x = a_x,med − a_x,modelo(sem arrasto)  (modelo = bias + braço da IMU
%      com p,q,r medidos; força específica de CG em x é ZERO no modelo atual)
%   3. regressão r_x = c0 + c1·u  → C_d = −c1/g ; idem em y. Reporta R² da
%      regressão, C_d, e quanto do erro de a_x, a_y desaparece.
%
% Uso:  >> drag_probe
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();  g = proj.g;

L  = load_log_data(fullfile(paths.data,'logs_concat.mat'));
dt = 0.1;
% janelas de treino (têm GPS) + o restante do 1º trecho com Spd > 0
WINS = [4 24; 25 41; 42 62; 63 99; 100 125; 130 180; 200 260];

t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min(268, min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]));
tg = (t_lo:dt:t_hi)';
ip = @(tt,xx) interp1(tt, xx, tg, 'linear');

p = ip(L.time_IMU,L.gyrX_raw); q = ip(L.time_IMU,L.gyrY_raw); r = ip(L.time_IMU,L.gyrZ_raw);
ax = ip(L.time_IMU,L.accX_raw); ay = ip(L.time_IMU,L.accY_raw); az = ip(L.time_IMU,L.accZ_raw);
phi = deg2rad(ip(L.time_ATT,L.roll_deg)); th = deg2rad(ip(L.time_ATT,L.pitch_deg)); psi = deg2rad(ip(L.time_ATT,L.yaw_deg));

% Fonte de velocidade: 'gps' (voo 4; ATENÇÃO: voo em ginásio, GPS com interferência)
% ou 'est' (estimate_velocity: atitude + barômetro, sem GPS). Como TODOS os voos
% foram em ambiente fechado (sem vento), 'est' é a referência mais defensável.
VEL_SRC = 'est';
if strcmp(VEL_SRC,'gps')
    tG = L.time_GPS(:);  spd = L.gps_spd(:);  gcrs = deg2rad(L.gps_gcrs(:));
    VN = interp1(tG, spd.*cos(gcrs), tg, 'linear', NaN);
    VE = interp1(tG, spd.*sin(gcrs), tg, 'linear', NaN);
    VD = interp1(tG, L.gps_vd(:),    tg, 'linear', NaN);
else
    VEst = estimate_velocity(L, tg);  VN = VEst.VN;  VE = VEst.VE;  VD = VEst.VD;
end
u = zeros(size(tg)); v = u; w = u;
for k = 1:numel(tg)
    cph=cos(phi(k)); sph=sin(phi(k)); cth=cos(th(k)); sth=sin(th(k)); cps=cos(psi(k)); sps=sin(psi(k));
    Rbn = [ cth*cps,             cth*sps,            -sth; ...
            sph*sth*cps-cph*sps, sph*sth*sps+cph*cps, sph*cth; ...
            cph*sth*cps+sph*sps, cph*sth*sps-sph*cps, cph*cth ];
    vb = Rbn*[VN(k); VE(k); VD(k)];
    u(k) = vb(1); v(k) = vb(2); w(k) = vb(3);
end

% modelo atual (sem arrasto): força específica de CG em x,y = 0; sobra bias + braço da IMU
pd = gradient(p, dt); qd = gradient(q, dt); rd = gradient(r, dt);
T_m = -az;                                  % irrelevante para x,y (só entra em fz)
[fx_mod, fy_mod, ~] = accelerometer_model(p, q, r, 0*u, 0*v, 0*w, T_m, pd, qd, rd, proj.imu_offset);
rx = ax - fx_mod;   ry = ay - fy_mod;       % resíduo horizontal do modelo atual

% máscara: dentro das janelas, GPS válido, aeronave em voo (Spd>0 em algum momento)
m = false(size(tg));
for i = 1:size(WINS,1), m = m | (tg>=WINS(i,1) & tg<=WINS(i,2)); end
m = m & isfinite(u) & isfinite(v) & isfinite(rx) & isfinite(ry);
fprintf('\n=== drag_probe [%s]: %d amostras válidas em %d janelas (0–268 s) ===\n', VEL_SRC, nnz(m), size(WINS,1));
fprintf('  |u| p50 %.2f p95 %.2f | |v| p50 %.2f p95 %.2f m/s\n', median(abs(u(m))), prctile(abs(u(m)),95), median(abs(v(m))), prctile(abs(v(m)),95));

% suavização leve (o acelerômetro tem vibração ~ >2 Hz que nenhum modelo de corpo rígido segue)
sg = @(x) sgolayfilt(x, 2, 7);
rxs = sg(rx); rys = sg(ry); us = sg(u); vs = sg(v);

% regressões  r = c0 + c1·vel  (c1 = −g·C_d)
reg = @(y,x) [ones(nnz(m),1) x(m)] \ y(m);
cx = reg(rxs, us);  cy = reg(rys, vs);
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
R2x = R2(rxs(m), cx(1)+cx(2)*us(m));   R2y = R2(rys(m), cy(1)+cy(2)*vs(m));
Cd_x = -cx(2)/g;  Cd_y = -cy(2)/g;
% erro-padrão do coeficiente
sex = @(y,x,c) sqrt(sum((y(m)-(c(1)+c(2)*x(m))).^2)/(nnz(m)-2) / sum((x(m)-mean(x(m))).^2));
fprintf('\n  a_x: resíduo = %+.3f %+.3f·u   (R² %.2f)  →  C_d = %.4f ± %.4f   [k_drag = m·g·C_d = %.3f N/(m/s)]\n', ...
    cx, R2x, Cd_x, sex(rxs,us,cx)/g, proj.m*g*Cd_x);
fprintf('  a_y: resíduo = %+.3f %+.3f·v   (R² %.2f)  →  C_d = %.4f ± %.4f   [k_drag = %.3f N/(m/s)]\n', ...
    cy, R2y, Cd_y, sex(rys,vs,cy)/g, proj.m*g*Cd_y);
fprintf('  correlação resíduo×velocidade: x %.2f | y %.2f\n', corr(rxs(m),us(m)), corr(rys(m),vs(m)));

% quanto do erro do acelerômetro some com o termo (mesmo C_d nos dois eixos)
Cd = mean([Cd_x, Cd_y]);
e0x = rms(rx(m)-mean(rx(m)));  e1x = rms(rx(m) + g*Cd*u(m) - mean(rx(m) + g*Cd*u(m)));
e0y = rms(ry(m)-mean(ry(m)));  e1y = rms(ry(m) + g*Cd*v(m) - mean(ry(m) + g*Cd*v(m)));
fprintf('\n  Com C_d médio = %.4f: RMS do resíduo (sem média) cai  a_x %.3f → %.3f m/s² (%.0f%%) | a_y %.3f → %.3f (%.0f%%)\n', ...
    Cd, e0x, e1x, 100*(1-e1x/e0x), e0y, e1y, 100*(1-e1y/e0y));

% referência física: força H do rotor (prior_damping): 4·k_h/m ≈ g·C_d
try
    pr = load(fullfile(paths.outputs,'prior_damping.mat')).prior;
    fprintf('  Referência teórica (força H, prior_damping): k_h = %.4f N/(m/s) por rotor → C_d ≈ 4·k_h/(m·g) = %.4f\n', ...
        pr.k_h, 4*pr.k_h/(proj.m*g));
catch, end
fprintf('  Referência Beard & McLain: valores típicos de C_d para pequenos multirrotores ~0,05–0,3 (adimensional em g·C_d [1/s]).\n');

%% ---------------- FIGURA ----------------
f = figure('Position',[60 60 1100 620],'Color','w'); try, f.Theme='light'; catch, end
subplot(2,2,1); scatter(us(m), rxs(m), 6, 'k', 'filled', 'MarkerFaceAlpha',0.25); hold on; grid on;
xx = linspace(min(us(m)),max(us(m)),50); plot(xx, cx(1)+cx(2)*xx, 'r-', 'LineWidth',2);
xlabel('u [m/s] (GPS → corpo)'); ylabel('resíduo a_x [m/s^2]'); title(sprintf('a_x: C_d = %.3f, R^2 = %.2f', Cd_x, R2x));
subplot(2,2,2); scatter(vs(m), rys(m), 6, 'k', 'filled', 'MarkerFaceAlpha',0.25); hold on; grid on;
xx = linspace(min(vs(m)),max(vs(m)),50); plot(xx, cy(1)+cy(2)*xx, 'r-', 'LineWidth',2);
xlabel('v [m/s]'); ylabel('resíduo a_y [m/s^2]'); title(sprintf('a_y: C_d = %.3f, R^2 = %.2f', Cd_y, R2y));
% série temporal numa janela de treino
sel = tg>=63 & tg<=99;
subplot(2,2,3); plot(tg(sel), rx(sel), 'k', tg(sel), -g*Cd*u(sel)+cx(1), 'r', 'LineWidth',1.2); grid on;
xlabel('t [s]'); ylabel('a_x [m/s^2]'); legend('resíduo medido','−g·C_d·u + c_0','Location','best'); title('janela 63–99 s');
subplot(2,2,4); plot(tg(sel), ry(sel), 'k', tg(sel), -g*Cd*v(sel)+cy(1), 'r', 'LineWidth',1.2); grid on;
xlabel('t [s]'); ylabel('a_y [m/s^2]'); legend('resíduo medido','−g·C_d·v + c_0','Location','best');
exportgraphics(f, fullfile(paths.images,'drag_probe.png'), 'BackgroundColor','white','Resolution',150);
fprintf('\n  Figura: %s\n', fullfile(paths.images,'drag_probe.png'));

%% ---------------- SEGUNDA HIPÓTESE: braço da IMU (termo de Euler α × r) ----------------
% Os picos do resíduo coincidem com os doublets (ṗ, q̇ grandes), não com u, v.
% eul_y = ṙ·rx − ṗ·rz ; eul_x = q̇·rz − ṙ·ry. Se rz (IMU acima/abaixo do CG) estiver
% errado, sobra −Δrz·ṗ em a_y e +Δrz·q̇ em a_x. Regressão conjunta [ṗ ou q̇, u ou v].
X_y = [ones(nnz(m),1), sg(pd(m)), vs(m)];   cyy = X_y \ rys(m);
X_x = [ones(nnz(m),1), sg(qd(m)), us(m)];   cxx = X_x \ rxs(m);
R2yy = R2(rys(m), X_y*cyy);  R2xx = R2(rxs(m), X_x*cxx);
fprintf('\n=== Hipótese 2: braço da IMU ===\n');
fprintf('  a_y: resíduo = %+.3f %+.4f·ṗ %+.3f·v   (R² %.2f)  →  Δrz = %+.3f m (rz atual %.3f)\n', cyy, R2yy, -cyy(2), proj.imu_offset(3));
fprintf('  a_x: resíduo = %+.3f %+.4f·q̇ %+.3f·u   (R² %.2f)  →  Δrz = %+.3f m\n', cxx, R2xx, cxx(2));
fprintf('  correlação resíduo×ṗ (y): %.2f | resíduo×q̇ (x): %.2f\n', corr(rys(m),sg(pd(m))), corr(rxs(m),sg(qd(m))));

%% ---------------- TERCEIRA HIPÓTESE: força horizontal proporcional ao COMANDO ----------------
% Picos de a_x, a_y sincronizados com os doublets → força ∝ momento comandado?
% Candidatos físicos: eixo do empuxo levemente inclinado (rotor i com inclinação
% ε_i produz f_h ≈ T_i·ε_i, e o DIFERENCIAL de empuxo que gera M_x, M_y vira força
% lateral ∝ M), batimento das pás, thrust tilt. Regressão do resíduo contra
% M_x, M_y, M_z calculados pela MESMA cadeia da identificação (motor_chain + P_final).
prog = load(fullfile(paths.outputs,'P_identified.mat'));  P = prog.P_final(:);
pwm4 = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
[rpm, fT, fQ] = motor_chain(tg, pwm4);
Tmr = fT(rpm) .* P(5:8)';   Qmr = fQ(rpm) .* P(9:12)';
dyn = vtol_dynamics('get_handles');
[Mx, My, Mz] = dyn.moments(Tmr, Qmr, proj.arms.Lx_r, proj.arms.Lx_l, proj.arms.Ly_f, proj.arms.Ly_r);
Ttot = sum(Tmr,2);
Mxs = sg(Mx); Mys = sg(My); Mzs = sg(Mz);
% regressões: cada eixo contra os três momentos (+ intercepto)
Xm = [ones(nnz(m),1) Mxs(m) Mys(m) Mzs(m)];
cmx = Xm \ rxs(m);   cmy = Xm \ rys(m);
R2mx = R2(rxs(m), Xm*cmx);   R2my = R2(rys(m), Xm*cmy);
fprintf('\n=== Hipótese 3: força horizontal ∝ momento comandado ===\n');
fprintf('  a_x: resíduo = %+.3f %+.3f·Mx %+.3f·My %+.3f·Mz   (R² %.2f)\n', cmx, R2mx);
fprintf('  a_y: resíduo = %+.3f %+.3f·Mx %+.3f·My %+.3f·Mz   (R² %.2f)\n', cmy, R2my);
fprintf('  correlações: a_x×Mx %.2f, a_x×My %.2f | a_y×Mx %.2f, a_y×My %.2f | a_x×ṗ %.2f, a_y×q̇ %.2f\n', ...
    corr(rxs(m),Mxs(m)), corr(rxs(m),Mys(m)), corr(rys(m),Mxs(m)), corr(rys(m),Mys(m)), corr(rxs(m),sg(pd(m))), corr(rys(m),sg(qd(m))));
% momento → força lateral equivale a um braço: f = M/ℓ_eq  →  ℓ_eq = 1/(m·coef)
fprintf('  Leitura física: coeficiente [m/s² por N·m] × massa = 1/braço equivalente. a_x×My: 1/(m·%.3f) = %.2f m | a_y×Mx: %.2f m\n', ...
    abs(cmx(3)), 1/(proj.m*max(abs(cmx(3)),1e-9)), 1/(proj.m*max(abs(cmy(2)),1e-9)));
% modelo combinado (tudo junto): quanto do resíduo é explicável linearmente?
Xall = [ones(nnz(m),1) us(m) vs(m) sg(pd(m)) sg(qd(m)) sg(rd(m)) Mxs(m) Mys(m) Mzs(m) sg(Ttot(m))];
cax = Xall \ rxs(m);  cay = Xall \ rys(m);
fprintf('  Regressão com TUDO (u,v,ṗ,q̇,ṙ,Mx,My,Mz,T): R² a_x %.2f | a_y %.2f  →  o resto é não-linear/ruído/vibração\n', ...
    R2(rxs(m),Xall*cax), R2(rys(m),Xall*cay));
% espectro do resíduo: onde está a energia?
[Pxx,fx] = pwelch(rx(m)-mean(rx(m)), 128, 64, 256, 1/dt);  [Pyy,fy] = pwelch(ry(m)-mean(ry(m)), 128, 64, 256, 1/dt);
cum = @(Pw,fw,f0) sum(Pw(fw<=f0))/sum(Pw);
fprintf('  Energia do resíduo abaixo de 1 Hz: a_x %.0f%% | a_y %.0f%%; abaixo de 2 Hz: a_x %.0f%% | a_y %.0f%% (Nyquist 5 Hz)\n', ...
    100*cum(Pxx,fx,1), 100*cum(Pyy,fy,1), 100*cum(Pxx,fx,2), 100*cum(Pyy,fy,2));

%% ---------------- Qual regressor carrega o R² de a_x? (empuxo total / a_z → desalinhamento da IMU) ----------------
Tts = sg(Ttot);  azs = sg(az);
fprintf('\n=== Hipótese 4: desalinhamento do eixo do sensor (a_x ∝ T ou a_z) ===\n');
fprintf('  correlação a_x×T %.2f | a_x×a_z %.2f | a_y×T %.2f | a_y×a_z %.2f\n', corr(rxs(m),Tts(m)), corr(rxs(m),azs(m)), corr(rys(m),Tts(m)), corr(rys(m),azs(m)));
cT_x = [ones(nnz(m),1) azs(m)] \ rxs(m);   cT_y = [ones(nnz(m),1) azs(m)] \ rys(m);
fprintf('  a_x = %+.3f %+.4f·a_z  (R² %.2f)  → ângulo equivalente de desalinhamento em torno de y: %.2f°\n', cT_x, R2(rxs(m),[ones(nnz(m),1) azs(m)]*cT_x), rad2deg(cT_x(2)));
fprintf('  a_y = %+.3f %+.4f·a_z  (R² %.2f)  → em torno de x: %.2f°\n', cT_y, R2(rys(m),[ones(nnz(m),1) azs(m)]*cT_y), rad2deg(-cT_y(2)));
% e com T do modelo (independente do sensor)
cTm_x = [ones(nnz(m),1) Tts(m)/proj.m] \ rxs(m);
fprintf('  a_x = %+.3f %+.4f·(T/m)  (R² %.2f)  → %.2f°\n', cTm_x, R2(rxs(m),[ones(nnz(m),1) Tts(m)/proj.m]*cTm_x), rad2deg(-cTm_x(2)));
