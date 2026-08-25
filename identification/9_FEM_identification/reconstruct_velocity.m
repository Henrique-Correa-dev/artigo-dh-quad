% reconstruct_velocity.m — Flight-Path Reconstruction (FPR) por fusão IMU+GPS
% =========================================================================
% Jategaonkar (2015), Cap. 10 (Data Compatibility / Flight-Path Reconstruction).
% Estima a velocidade NED (e u,v,w body) fundindo:
%   • ACELERÔMETRO (integra, via atitude medida + gravidade) → suave, alta taxa
%   • GPS velocidade (corrige a deriva do accel)              → âncora absoluta
%   • bias do accel estimado como estado (random walk)
% Saída = velocidade MODEL-INDEPENDENT (não usa k_T/arrasto), sem drift onde há GPS.
%
% ⚠️ Só confiável onde há GPS (< ~268 s). Além disso, vira predição pura → deriva.
%
% EKF (loosely-coupled INS/GPS):
%   estado x = [Vn,Ve,Vd, b_ax,b_ay,b_az]  (6)
%   pred:  V̇ = R_body→ned(att)·(a_medido − b) + g_ned ;  ḃ = 0
%   upd:   z = V_gps ,  H = [I_3  0_3]
%
% Saída salva: outputs/uvw_fpr.mat  (tg, u,v,w, V_ned, valid_gps)
% Uso:  >> reconstruct_velocity
% =========================================================================

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));  paths = setup_paths();
out_dir = fullfile(here,'outputs'); if ~exist(out_dir,'dir'), mkdir(out_dir); end
img_dir = fullfile(out_dir,'images'); if ~exist(img_dir,'dir'), mkdir(img_dir); end

%% ---- CONFIG ----
LOG_FILE = 'logs_concat.mat';
dt   = 0.05;                 % grade fina (FPR ganha com taxa alta do IMU)
g    = 9.81;
T_MAX = 130;                 % corta em t≤T_MAX (cobre as janelas de treino 4–125s)
% Ruídos (tuning)
SIG_A   = 0.5;               % ruído do accel → processo na velocidade (m/s²/√Hz)
SIG_B   = 0.02;              % random walk do bias do accel (m/s²/√s·√Hz)
SIG_GPS = 0.4;               % ruído da velocidade GPS (m/s)
B0_P    = 1.0;               % incerteza inicial do bias (m/s²)

%% ---- carregar ----
d = load(fullfile(paths.data, LOG_FILE));  L = d.L;
t0 = max([min(L.time_IMU),min(L.time_ATT)]);
t1 = min([max(L.time_IMU),max(L.time_ATT)]);
% CLIPA: sem GPS a FPR deriva (polui). E corta em T_MAX (janelas de treino).
if isfield(L,'time_GPS') && ~isempty(L.time_GPS)
    t1 = min(t1, max(L.time_GPS) + 1);     % +1s só p/ ver a linha de "fim GPS"
end
t1 = min(t1, T_MAX);
tg = (t0:dt:t1)';  N = numel(tg);
ip = @(t,x) interp1(t,x,tg,'linear','extrap');
ax = ip(L.time_IMU,L.accX_raw); ay = ip(L.time_IMU,L.accY_raw); az = ip(L.time_IMU,L.accZ_raw);
phi = deg2rad(ip(L.time_ATT,L.roll_deg)); th = deg2rad(ip(L.time_ATT,L.pitch_deg)); ps = deg2rad(ip(L.time_ATT,L.yaw_deg));

% GPS: interpola na grade SÓ dentro da cobertura (fora → sem update)
hasG = isfield(L,'gps_vn') && ~isempty(L.gps_vn) && ~isempty(L.time_GPS);
if hasG
    tGmin = min(L.time_GPS); tGmax = max(L.time_GPS);
    vn = interp1(L.time_GPS,L.gps_vn,tg,'linear'); ve = interp1(L.time_GPS,L.gps_ve,tg,'linear'); vd = interp1(L.time_GPS,L.gps_vd,tg,'linear');
    valid_gps = (tg>=tGmin) & (tg<=tGmax) & isfinite(vn);
else
    error('reconstruct_velocity:noGPS','logs_concat sem velocidade GPS.');
end
fprintf('=== FPR IMU+GPS | [%.0f,%.0f]s | dt=%.2f | GPS em [%.0f,%.0f]s (%d/%d amostras) ===\n', ...
    t0,t1,dt, tGmin,tGmax, nnz(valid_gps), N);

%% ---- EKF ----
g_ned = [0;0;g];
x = zeros(6,1);                                   % [V(3); bias(3)]
P = diag([1 1 1, B0_P B0_P B0_P].^2);
Qc = diag([SIG_A SIG_A SIG_A, SIG_B SIG_B SIG_B].^2);   % densidade espectral
Rg = (SIG_GPS^2)*eye(3);
Hk = [eye(3) zeros(3)];
Vned = zeros(N,3); bias = zeros(N,3); innov = nan(N,3);
for k=1:N
    Rbn = dcm_body2ned(phi(k),th(k),ps(k));
    am  = [ax(k);ay(k);az(k)];
    % --- predição ---
    Vdot = Rbn*(am - x(4:6)) + g_ned;
    x(1:3) = x(1:3) + dt*Vdot;                    % bias: ḃ=0
    F = [zeros(3) -Rbn; zeros(3) zeros(3)];
    Phi = eye(6) + F*dt;
    P = Phi*P*Phi' + Qc*dt;
    % --- update (se GPS) ---
    if valid_gps(k)
        z = [vn(k);ve(k);vd(k)];
        y = z - Hk*x;  S = Hk*P*Hk' + Rg;  K = P*Hk'/S;
        x = x + K*y;  P = (eye(6)-K*Hk)*P;
        innov(k,:) = y';
    end
    Vned(k,:) = x(1:3)';  bias(k,:) = x(4:6)';
end

%% ---- velocidade body (u,v,w) = R_ned→body · V_ned ----
u=zeros(N,1); v=zeros(N,1); w=zeros(N,1);
for k=1:N
    Rnb = dcm_body2ned(phi(k),th(k),ps(k))';      % ned→body
    uvw = Rnb*Vned(k,:)'; u(k)=uvw(1); v(k)=uvw(2); w(k)=uvw(3);
end

%% ---- relatório ----
m=valid_gps;
fprintf('  bias accel estimado (fim, janela GPS): [%.3f %.3f %.3f] m/s²\n', bias(find(m,1,'last'),:));
fprintf('  |V| (corpo) na janela GPS: média %.2f  máx %.2f m/s\n', mean(vecnorm([u(m) v(m) w(m)],2,2)), max(vecnorm([u(m) v(m) w(m)],2,2)));
rms_nan = @(x) sqrt(mean(x(~isnan(x)).^2));    % RMS ignorando NaN (sem toolbox)
fprintf('  u (frente): média %.2f  máx %.2f m/s | RMS inovação V: [%.2f %.2f %.2f] m/s\n', ...
    mean(u(m)), max(abs(u(m))), rms_nan(innov(:,1)),rms_nan(innov(:,2)),rms_nan(innov(:,3)));

%% ---- salvar (pro OEM/FEM usar) ----
save(fullfile(out_dir,'uvw_fpr.mat'),'tg','u','v','w','Vned','valid_gps','bias','dt');
fprintf('  Salvo: %s  (use no OEM como u,v,w medidos, só onde valid_gps)\n', fullfile(out_dir,'uvw_fpr.mat'));

%% ---- plots ----
fig=figure('Name','FPR IMU+GPS','Position',[40 40 1300 820]);
lab={'V_N','V_E','V_D'}; gpsv={vn,ve,vd};
for c=1:3
    subplot(3,3,c); hold on; grid on;
    plot(tg(m),gpsv{c}(m),'.','Color',[.6 .6 .6],'DisplayName','GPS cru');
    plot(tg,Vned(:,c),'b','LineWidth',1.2,'DisplayName','FPR (fundido)');
    ylabel(lab{c}); if c==1, legend('Location','best'); title('NED: GPS cru vs FPR'); end
    if tGmax<=tg(end), xline(tGmax,'r--','HandleVisibility','off'); end
    xlim([tg(1) tg(end)]);
end
uvw={u,v,w}; luvw={'u','v','w'};
for c=1:3
    subplot(3,3,3+c); hold on; grid on;
    plot(tg,uvw{c},'b','LineWidth',1.2);
    if tGmax<=tg(end), xline(tGmax,'r--','fim GPS'); end
    xlim([tg(1) tg(end)]);
    ylabel([luvw{c} ' (m/s)']); if c==1, title('Body u,v,w (FPR)'); end
end
for c=1:3
    subplot(3,3,6+c); hold on; grid on;
    plot(tg,bias(:,c),'m','LineWidth',1.2);
    if tGmax<=tg(end), xline(tGmax,'r--'); end
    xlim([tg(1) tg(end)]);
    ylabel(['bias a' luvw{c} ' (m/s²)']); xlabel('t (s)');
    if c==1, title('Bias do accel estimado'); end
end
saveas(fig,fullfile(img_dir,'fpr_velocity.png'));
fprintf('  Figura: %s\n', fullfile(img_dir,'fpr_velocity.png'));

%% ---- helper: DCM body→ned (3-2-1) ----
function Rbn = dcm_body2ned(phi,th,ps)
    cphi=cos(phi);sphi=sin(phi);cth=cos(th);sth=sin(th);cps=cos(ps);sps=sin(ps);
    Rnb = [ cth*cps,                cth*sps,               -sth;
            sphi*sth*cps-cphi*sps,  sphi*sth*sps+cphi*cps,  sphi*cth;
            cphi*sth*cps+sphi*sps,  cphi*sth*sps-sphi*cps,  cphi*cth];
    Rbn = Rnb';
end
