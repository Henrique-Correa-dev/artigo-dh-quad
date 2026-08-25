% identify_kdrag.m — Identifica o arrasto translacional k_drag por regressão
% =========================================================================
% Usa a velocidade da FPR (reconstruct_velocity → uvw_fpr.mat, model-independent)
% nas janelas de treino (<130 s, com GPS). Modelo do livro (Eq. 6.31):
%
%   accX_medido = braço_x(p,q,r,ṗ,q̇,ṙ) − (k_drag/m)·u + bias_x
%   accY_medido = braço_y(...)          − (k_drag/m)·v + bias_y
%
% → resíduo_x = accX − braço_x  =  −(k_drag/m)·u + bias_x   (idem y)
% Regressão linear conjunta [k_drag/m, bias_x, bias_y] → CRB exato (mín. quadrados).
%
% Uso:  >> identify_kdrag   (rode reconstruct_velocity antes)
% =========================================================================

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));  paths = setup_paths();
img_dir = fullfile(here,'outputs','images'); if ~exist(img_dir,'dir'), mkdir(img_dir); end

%% ---- CONFIG ----
LOG_FILE = 'logs_concat.mat';
t_trains = {[4,24]; [25,41]; [42,62]; [63,99]; [100,125]};   % <130 s (FPR válida)
dt = 0.1;  SG_ORDER=2; SG_FRAME=7;

%% ---- carregar log + FPR ----
proj = parameters();  m = proj.m;  r_imu = proj.imu_offset;
rx=r_imu(1); ry=r_imu(2); rz=r_imu(3);
L = load_log_data(fullfile(paths.data, LOG_FILE));
F = load(fullfile(here,'outputs','uvw_fpr.mat'));   % tg,u,v,w,valid_gps

t0=max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]);
t1=min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]);
tg=(t0:dt:t1)';
ip=@(t,x)interp1(t,x,tg,'linear'); sg=@(x)sgolayfilt(x,SG_ORDER,SG_FRAME);
p=sg(ip(L.time_IMU,L.gyrX_raw)); q=sg(ip(L.time_IMU,L.gyrY_raw)); r=sg(ip(L.time_IMU,L.gyrZ_raw));
axm=sg(ip(L.time_IMU,L.accX_raw)); aym=sg(ip(L.time_IMU,L.accY_raw));
pd=gradient(p,dt); qd=gradient(q,dt); rd=gradient(r,dt);
% velocidade FPR interpolada na grade (só onde válida)
uF=interp1(F.tg,F.u,tg,'linear'); vF=interp1(F.tg,F.v,tg,'linear');
validF=interp1(F.tg,double(F.valid_gps),tg,'linear')>0.5 & isfinite(uF);

% braço do IMU (eul + cen) — mesma fórmula do accelerometer_model
lev_x = qd*rz - rd*ry + p.*q*ry + p.*r*rz - (q.^2+r.^2)*rx;
lev_y = rd*rx - pd*rz + p.*q*rx + q.*r*rz - (p.^2+r.^2)*ry;
res_x = axm - lev_x;          % deve = -(k_drag/m)·u + bias_x
res_y = aym - lev_y;          % deve = -(k_drag/m)·v + bias_y

%% ---- monta a regressão nas janelas de treino ----
A=[]; b=[]; Uall=[]; Vall=[]; RXall=[]; RYall=[];
for s=1:numel(t_trains)
    sel=(tg>=t_trains{s}(1))&(tg<=t_trains{s}(2))&validF;
    uu=uF(sel); vv=vF(sel); rxw=res_x(sel); ryw=res_y(sel); n=nnz(sel);
    % linhas x: [-u, 1, 0]; linhas y: [-v, 0, 1]
    A=[A; [-uu, ones(n,1), zeros(n,1)]; [-vv, zeros(n,1), ones(n,1)]]; %#ok<AGROW>
    b=[b; rxw; ryw]; %#ok<AGROW>
    Uall=[Uall;uu]; Vall=[Vall;vv]; RXall=[RXall;rxw]; RYall=[RYall;ryw]; %#ok<AGROW>
end
Ntot=numel(b);
fprintf('=== identify_kdrag | %d janelas | %d pontos | FPR u,v,w ===\n', numel(t_trains), Ntot/2);
fprintf('  faixa de u: [%.2f, %.2f] m/s | v: [%.2f, %.2f] m/s\n', min(Uall),max(Uall),min(Vall),max(Vall));

%% ---- mínimos quadrados + CRB ----
theta = A\b;                         % [k_drag/m; bias_x; bias_y]
resid = b - A*theta;
sig2  = (resid'*resid)/(Ntot-3);
Cov   = sig2 * inv(A'*A);            %#ok<MINV>
se    = sqrt(diag(Cov));
kdm   = theta(1);  k_drag = kdm*m;
se_kd = se(1)*m;
fprintf('\n--- Resultado ---\n');
fprintf('  k_drag   = %.4f  N/(m/s)   (std %.4f → %.1f%%)\n', k_drag, se_kd, 100*se_kd/max(abs(k_drag),1e-9));
fprintf('  bias_x   = %.4f m/s²       (FPR ~-0.18 | hardcoded -0.30)\n', theta(2));
fprintf('  bias_y   = %.4f m/s²       (FPR ~-0.20 | hardcoded -0.20)\n', theta(3));
fprintf('  σ resíduo= %.3f m/s² | R² da regressão = %.3f\n', sqrt(sig2), 1 - (resid'*resid)/sum((b-mean(b)).^2));

crbpc = 100*se_kd/max(abs(k_drag),1e-9);
if crbpc < 20
    fprintf('  >>> CRB < 20%% → k_drag IDENTIFICÁVEL (Sato). Arrasto capturado!\n');
else
    fprintf('  >>> CRB %.0f%% > 20%% → fraco (esperado: hover, u pequeno). Documenta a tentativa.\n', crbpc);
end

%% ---- plot: resíduo vs velocidade + reta ----
fig=figure('Name','identify_kdrag','Position',[60 60 1150 480]);
uu=linspace(min(Uall),max(Uall),50); vv=linspace(min(Vall),max(Vall),50);
subplot(1,2,1); hold on; grid on;
plot(Uall,RXall,'.','Color',[.5 .6 1],'MarkerSize',4);
plot(uu, -kdm*uu+theta(2),'r-','LineWidth',2);
xlabel('u (m/s)'); ylabel('accX − braço  (m/s²)'); title(sprintf('ROLL: inclinação = −k_{drag}/m  (k_{drag}=%.3f, %.0f%%)',k_drag,crbpc));
subplot(1,2,2); hold on; grid on;
plot(Vall,RYall,'.','Color',[1 .6 .6],'MarkerSize',4);
plot(vv, -kdm*vv+theta(3),'r-','LineWidth',2);
xlabel('v (m/s)'); ylabel('accY − braço  (m/s²)'); title('PITCH/lateral: inclinação = −k_{drag}/m');
saveas(fig,fullfile(img_dir,'kdrag_regression.png'));
fprintf('  Figura: %s\n', fullfile(img_dir,'kdrag_regression.png'));
