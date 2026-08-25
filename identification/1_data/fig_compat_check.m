% fig_compat_check.m — Verificação de compatibilidade cinemática (Cap3 da tese)
% =========================================================================
% Compara TRÊS estimativas de atitude (phi, theta, psi) numa janela de voo:
%   (1) GIRO puro      — integra a cinemática de Euler a partir do pqr corrigido
%                        (sem bias) e filtrado (SG). Tende a derivar no tempo.
%   (2) EKF do Pixhawk — referência (fusão a bordo; yaw via GSF, compass off).
%   (3) Complementar   — funde acelerômetro (direção da gravidade -> roll/pitch)
%                        com o giro. O accel NÃO observa o yaw, então em psi o
%                        complementar coincide com a integração pura (só o EKF/GSF
%                        segura o yaw) — isso explica o yaw ser o canal fraco.
%
% Critério (Jategaonkar cap.10): se a atitude reconstruída da cinemática bate com
% a referência, os dados (pqr, accel, atitude) são cinematicamente CONSISTENTES.
%
% Uso:  >> fig_compat_check

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();
set(groot, 'defaultAxesFontSize', 11);

%% ---------------- CONFIG ----------------
LOG_FILE = 'logs_concat.mat';
dt = 0.1;  g = 9.81;  SG_ORDER = 2;  SG_FRAME = 7;
W     = [605 665];                 % janela de voo (contígua, com manobras)
W_CAL = [136 150];                 % calibração (drone parado) p/ o bias
TAU   = 1.0;                       % constante do filtro complementar [s]
DEST = {paths.images, '/Users/graest/ita-master/artigo/artigo-dh-quad/dissertacao/Cap3'};

%% ---------------- CARREGAR + BIAS ----------------
L = load_log_data(fullfile(paths.data, LOG_FILE));
tIMU = L.time_IMU;
mcal = (tIMU>=W_CAL(1)) & (tIMU<=W_CAL(2));
biasG = [mean(L.gyrX_raw(mcal)); mean(L.gyrY_raw(mcal)); mean(L.gyrZ_raw(mcal))];
macA = (L.time_ATT>=W_CAL(1)) & (L.time_ATT<=W_CAL(2));
phic = deg2rad(mean(L.roll_deg(macA)));  thc = deg2rad(mean(L.pitch_deg(macA)));
exp_acc = [ g*sin(thc); -g*sin(phic)*cos(thc); -g*cos(phic)*cos(thc) ];
biasA = [mean(L.accX_raw(mcal)); mean(L.accY_raw(mcal)); mean(L.accZ_raw(mcal))] - exp_acc;

%% ---------------- GRADE + SINAIS NA JANELA ----------------
tg = (W(1):dt:W(2))';  N = numel(tg);
ipw = @(t,x) interp1(t, x, tg, 'linear');
if exist('sgolayfilt','file')==2, sgf=@(x) sgolayfilt(x,SG_ORDER,SG_FRAME);
else, warning('sgolayfilt ausente — SG local'); sgf=@(x) sg_local(x,SG_ORDER,SG_FRAME); end
mwin = (tIMU>=W(1)) & (tIMU<=W(2));
if max(diff(tIMU(mwin))) > 0.2, warning('a janela [%g %g] cruza um gap entre logs.', W); end

% giro corrigido (sem bias) + filtrado
p = sgf(ipw(tIMU, L.gyrX_raw) - biasG(1));
q = sgf(ipw(tIMU, L.gyrY_raw) - biasG(2));
r = sgf(ipw(tIMU, L.gyrZ_raw) - biasG(3));
% giro RAW (sem bias, sem filtro) — mostra o efeito da correção
praw = ipw(tIMU, L.gyrX_raw);  qraw = ipw(tIMU, L.gyrY_raw);  rraw = ipw(tIMU, L.gyrZ_raw);
% acelerômetro (sem bias do solo) + filtrado — para a direção da gravidade
ax = sgf(ipw(tIMU, L.accX_raw) - biasA(1));
ay = sgf(ipw(tIMU, L.accY_raw) - biasA(2));
az = sgf(ipw(tIMU, L.accZ_raw) - biasA(3));
% atitude do EKF (graus) — referência
roll_e = ipw(L.time_ATT, L.roll_deg);  pitch_e = ipw(L.time_ATT, L.pitch_deg);  yaw_e = ipw(L.time_ATT, L.yaw_deg);

%% ---------------- (1) INTEGRAÇÃO DO GIRO (RK4): RAW e corrigido ----------------
att0 = [deg2rad(roll_e(1)); deg2rad(pitch_e(1)); deg2rad(yaw_e(1))];   % IC = EKF em t0
[PHI,  TH,  PSI ] = integ_att(att0, p,    q,    r,    dt);   % corrigido (sem bias) + filtrado
[PHIr, THr, PSIr] = integ_att(att0, praw, qraw, rraw, dt);   % RAW (sem bias, sem filtro)

%% ---------------- (3) FILTRO COMPLEMENTAR (accel + giro) ----------------
alpha = TAU/(TAU+dt);
ac = [deg2rad(roll_e(1)); deg2rad(pitch_e(1)); deg2rad(yaw_e(1))];
PHc = zeros(N,1); THc = zeros(N,1); PSc = zeros(N,1);
PHc(1)=ac(1); THc(1)=ac(2); PSc(1)=ac(3);
for k = 2:N
    pred = ac + dt*euler_kin(ac, [p(k);q(k);r(k)]);     % predição pelo giro
    phi_a = atan2(-ay(k), -az(k));                       % roll/pitch pela gravidade
    th_a  = atan2( ax(k), hypot(ay(k), az(k)) );
    ac = [ alpha*pred(1) + (1-alpha)*phi_a;              % funde roll
           alpha*pred(2) + (1-alpha)*th_a;               % funde pitch
           pred(3) ];                                    % yaw: só giro (accel não observa)
    PHc(k)=ac(1); THc(k)=ac(2); PSc(k)=ac(3);
end

%% ---------------- PLOT ----------------
phi_g = rad2deg(PHI);   th_g = rad2deg(TH);   psi_g = rad2deg(unwrap(PSI));
phi_r = rad2deg(PHIr);  th_r = rad2deg(THr);  psi_r = rad2deg(unwrap(PSIr));
phi_c = rad2deg(PHc);   th_c = rad2deg(THc);  psi_c = rad2deg(unwrap(PSc));
yaw_e_u = rad2deg(unwrap(deg2rad(yaw_e)));
S = {phi_g, roll_e,  phi_c,  phi_r, 'φ (rolagem)'; ...
     th_g,  pitch_e, th_c,   th_r,  'θ (arfagem)'; ...
     psi_g, yaw_e_u, psi_c,  psi_r, 'ψ (guinada)'};
rms = @(a,b) sqrt(mean((a-b).^2));
fig = figure('Color','w','Position',[60 60 1250 820]);
for i = 1:3
    subplot(3,1,i); hold on; grid on;
    plot(tg, S{i,2}, 'k-', 'LineWidth', 1.8);                          % EKF (referência)
    plot(tg, S{i,4}, '-', 'Color', [0.92 0.6 0],  'LineWidth', 1.1);   % giro RAW
    plot(tg, S{i,1}, '-', 'Color', [0 0.45 0.9],  'LineWidth', 1.3);   % giro corrigido+filtrado
    plot(tg, S{i,3}, '-', 'Color', [0.85 0.2 0.2],'LineWidth', 1.3);   % complementar
    ylabel([S{i,5} ' [deg]']);
    text(0.01, 0.90, sprintf('RMS vs EKF:  raw = %.1f°   corrig = %.1f°   compl = %.1f°', ...
        rms(S{i,4},S{i,2}), rms(S{i,1},S{i,2}), rms(S{i,3},S{i,2})), ...
        'Units','normalized', 'FontSize', 9, 'BackgroundColor', [1 1 1 .6]);
    if i==1  % (sem título — a legenda do LaTeX descreve a figura)
        legend({'EKF (Pixhawk)','giro RAW','giro corrigido+filtrado','complementar (accel+giro)'}, 'Location','best'); end
    if i==3, xlabel('t [s]'); end
end
save_to(fig, 'compat_check.png', DEST);

fprintf('Janela [%g %g]s | bias giro=[%+.4f %+.4f %+.4f] | RMS vs EKF (raw|corrig|compl) deg:\n', W, biasG);
fprintf('  phi:   %.1f | %.1f | %.1f\n  theta: %.1f | %.1f | %.1f\n  psi:   %.1f | %.1f | %.1f\n', ...
    rms(phi_r,roll_e),  rms(phi_g,roll_e),  rms(phi_c,roll_e), ...
    rms(th_r,pitch_e),  rms(th_g,pitch_e),  rms(th_c,pitch_e), ...
    rms(psi_r,yaw_e_u), rms(psi_g,yaw_e_u), rms(psi_c,yaw_e_u));
set(groot, 'defaultAxesFontSize', 'remove');

%% ---------------- HELPERS ----------------
function d = euler_kin(att, pqr)
    phi = att(1); th = att(2);  p = pqr(1); q = pqr(2); r = pqr(3);
    sp = sin(phi); cp = cos(phi); tt = tan(th); ct = cos(th);
    d = [ p + (q*sp + r*cp)*tt;  q*cp - r*sp;  (q*sp + r*cp)/ct ];
end

function [PHI, TH, PSI] = integ_att(att0, p, q, r, dt)
% Integra a cinemática de Euler (RK4) a partir de att0, com taxas p,q,r na grade.
    N = numel(p);  PHI = zeros(N,1); TH = zeros(N,1); PSI = zeros(N,1);
    att = att0(:);  PHI(1)=att(1); TH(1)=att(2); PSI(1)=att(3);
    for k = 1:N-1
        rk=[p(k);q(k);r(k)]; rk1=[p(k+1);q(k+1);r(k+1)]; rm=0.5*(rk+rk1);
        k1=euler_kin(att,rk); k2=euler_kin(att+0.5*dt*k1,rm);
        k3=euler_kin(att+0.5*dt*k2,rm); k4=euler_kin(att+dt*k3,rk1);
        att = att + dt/6*(k1+2*k2+2*k3+k4);
        PHI(k+1)=att(1); TH(k+1)=att(2); PSI(k+1)=att(3);
    end
end

function save_to(fig, name, dests)
    lighten(fig);
    for i = 1:numel(dests)
        if ~exist(dests{i}, 'dir'), mkdir(dests{i}); end
        try, exportgraphics(fig, fullfile(dests{i},name), 'Resolution', 200);
        catch, saveas(fig, fullfile(dests{i},name)); end
    end
    fprintf('  salvo: %s\n', name);
end

function lighten(fig)
    set(fig, 'Color', 'w');
    set(findall(fig,'Type','axes'), 'Color','w','XColor','k','YColor','k', 'GridColor',[.15 .15 .15]);
    set(findall(fig,'Type','text'), 'Color', 'k');
    lg = findall(fig,'Type','legend'); for i=1:numel(lg), set(lg(i),'TextColor','k','Color','w'); end
end

function y = sg_local(x, order, frame)
    x = x(:); m = (frame-1)/2; t = (-m:m)'; A = t.^(0:order);
    Hc = A * pinv(A); c = Hc(m+1, :)'; y = conv(x, flipud(c), 'same');
end
