%% integra_giro_euler.m
%  TESTE SEPARADO — Atitude por INTEGRAÇÃO PURA do giroscópio.
%
%  Integra a cinemática de Euler (sequência 3-2-1) usando as taxas p,q,r
%  medidas pela IMU e REALIMENTA a própria atitude calculada: a derivada
%  depende de phi/theta atuais, então o estado integrado entra de volta na
%  cinemática a cada passo (dead-reckoning de malha fechada no próprio estado).
%
%  NÃO há correção de sensor nenhum (sem accel, sem mag, sem GPS). É o oposto
%  do EKF: mostra o que acontece quando você confia 100%% no giroscópio.
%
%  Cinemática (mesma usada na dinâmica do VTOL):
%     phi_dot   = p + (q*sin(phi) + r*cos(phi)) * tan(theta)
%     theta_dot =      q*cos(phi) - r*sin(phi)
%     psi_dot   =     (q*sin(phi) + r*cos(phi)) / cos(theta)
%
%  Esperado: deriva crescente = (viés do giro)*t + ruído integrado. Como não há
%  referência absoluta, os 3 eixos derivam; o resultado afasta-se do ATK do log.
%
%  Dados:  ../1_data/logs_concat.mat   (struct L: IMU 25 Hz + ATT 10 Hz)
%  Janela: 1 a 110 s
%  Rodar:  >> integra_giro_euler

clear; clc; close all;
HERE = fileparts(mfilename('fullpath'));
LOGF = fullfile(HERE, '..', '1_data', 'logs_concat.mat');
W    = [1 110];      % janela [s]
DT   = 0.005;        % passo de integração (fino p/ RK4)

%% ---- carrega e recorta ----
S = load(LOGF);  L = S.L;
im = L.time_IMU >= W(1) & L.time_IMU <= W(2);
at = L.time_ATT >= W(1) & L.time_ATT <= W(2);

t_imu = L.time_IMU(im);
p = L.gyrX_raw(im);  q = L.gyrY_raw(im);  r = L.gyrZ_raw(im);   % rad/s

t_att     = L.time_ATT(at);
roll_log  = L.roll_deg(at);  pitch_log = L.pitch_deg(at);  yaw_log = L.yaw_deg(at);

fprintf('Janela [%g, %g] s | IMU: %d amostras | ATT: %d amostras\n', ...
    W(1), W(2), numel(t_imu), numel(t_att));

%% ---- grade de integração + interpolação do giro ----
tg = (W(1):DT:W(2))';  Ng = numel(tg);
pg = interp1(t_imu, p, tg, 'linear', 'extrap');
qg = interp1(t_imu, q, tg, 'linear', 'extrap');
rg = interp1(t_imu, r, tg, 'linear', 'extrap');

% Condição inicial = atitude do log (ATT) no instante inicial
phi0   = deg2rad(interp1(t_att, roll_log,  tg(1), 'linear', 'extrap'));
theta0 = deg2rad(interp1(t_att, pitch_log, tg(1), 'linear', 'extrap'));
psi0   = deg2rad(interp1(t_att, yaw_log,   tg(1), 'linear', 'extrap'));
fprintf('CI (do ATT em t=%.1fs): phi=%.2f  theta=%.2f  psi=%.2f deg\n', ...
    tg(1), rad2deg(phi0), rad2deg(theta0), rad2deg(psi0));

%% ---- integração RK4 (realimentando a própria atitude) ----
PHI = zeros(Ng,1);  TH = zeros(Ng,1);  PS = zeros(Ng,1);
PHI(1) = phi0;  TH(1) = theta0;  PS(1) = psi0;

pqr = @(tt)[interp1(tg,pg,tt,'linear','extrap'); ...
            interp1(tg,qg,tt,'linear','extrap'); ...
            interp1(tg,rg,tt,'linear','extrap')];
f = @(tt, y) euler_kin(y, pqr(tt));     % y = [phi;theta;psi] realimentado

for k = 1:Ng-1
    y = [PHI(k); TH(k); PS(k)];  h = DT;  tt = tg(k);
    k1 = f(tt,       y);
    k2 = f(tt+h/2,   y + h/2*k1);
    k3 = f(tt+h/2,   y + h/2*k2);
    k4 = f(tt+h,     y + h*k3);
    yn = y + h/6*(k1 + 2*k2 + 2*k3 + k4);
    PHI(k+1) = yn(1);  TH(k+1) = yn(2);  PS(k+1) = yn(3);
end

%% ---- métricas (vs ATT interpolado na grade) ----
rl = interp1(t_att, roll_log,  tg, 'linear', 'extrap');
pl = interp1(t_att, pitch_log, tg, 'linear', 'extrap');
yl = interp1(t_att, yaw_log,   tg, 'linear', 'extrap');
e_phi = rad2deg(PHI) - rl;
e_th  = rad2deg(TH)  - pl;
e_ps  = wrap180(rad2deg(PS) - yl);

fprintf('\n=== Integração pura do giro vs ATT (deriva esperada) ===\n');
fprintf('  roll  : RMS=%6.2f deg | deriva final=%+7.2f deg\n', rms(e_phi), e_phi(end));
fprintf('  pitch : RMS=%6.2f deg | deriva final=%+7.2f deg\n', rms(e_th),  e_th(end));
fprintf('  yaw   : RMS=%6.2f deg | deriva final=%+7.2f deg\n', rms(e_ps),  e_ps(end));

%% ---- gráfico ----
fig = figure('Color','w','Position',[60 60 1250 760]);
nm  = {'\phi (roll)','\theta (pitch)','\psi (yaw)'};
sim = {rad2deg(PHI), rad2deg(TH), rad2deg(PS)};
log = {rl, pl, yl};
err = {e_phi, e_th, e_ps};
for i = 1:3
    subplot(3,2,2*i-1); hold on; grid on;
    plot(tg, log{i}, 'b-',  'LineWidth',1.3, 'DisplayName','ATT (log / EKF Pixhawk)');
    plot(tg, sim{i}, 'r--', 'LineWidth',1.3, 'DisplayName','Integração pura do giro');
    ylabel([nm{i} ' [°]']);
    if i==1, legend('Location','best'); title('Atitude: integração do giro vs log'); end
    if i==3, xlabel('t [s]'); end

    subplot(3,2,2*i); hold on; grid on;
    plot(tg, err{i}, 'm-', 'LineWidth',1); yline(0,'k--');
    ylabel(['\Delta' nm{i} ' [°]']);
    title(sprintf('Erro (RMS=%.2f°, final=%+.2f°)', rms(err{i}), err{i}(end)));
    if i==3, xlabel('t [s]'); end
end
sgtitle('Integração pura da cinemática de Euler (p,q,r) — sem correção de sensor');

out = fullfile(HERE, 'integra_giro_euler_cmp.png');
saveas(fig, out);
fprintf('\nFigura salva: %s\n', out);

%% ===================== subfunções =====================
function dy = euler_kin(y, pqr)
    phi = y(1);  theta = y(2);
    w = pqr;  p = w(1);  q = w(2);  r = w(3);
    sp = sin(phi);  cp = cos(phi);  ct = cos(theta);
    if abs(ct) < 1e-7, ct = 1e-7*sign(ct + 1e-12); end
    tt = sin(theta)/ct;
    dy = [ p + (q*sp + r*cp)*tt;
               q*cp - r*sp;
              (q*sp + r*cp)/ct ];
end

function a = wrap180(a)
    a = mod(a + 180, 360) - 180;
end
