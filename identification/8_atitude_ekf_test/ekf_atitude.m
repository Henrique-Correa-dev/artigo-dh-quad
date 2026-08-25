%% ekf_atitude.m
%  Recriação do EKF de atitude do Pixhawk (ArduPilot) — giro + accel + MAG.
%
%  ─── O QUE O PIXHAWK FAZ (pesquisado) ───────────────────────────────────────
%  ArduPilot EKF3: EKF de 24 estados, formulação error-state com vetor de rotação
%  aplicado ao quatérnio (ref.: Pittelkau, "Rotation Vector in Attitude Estimation",
%  J. Guidance 2003), fundindo IMU (giro+accel) como ENTRADA da predição e
%  GPS + magnetômetro + barômetro como MEDIDAS de correção. ATT = saída do filtro.
%
%  ─── O QUE ESTE SCRIPT RECRIA ────────────────────────────────────────────────
%  EKF quaterniônico de 7 estados  x = [q0 q1 q2 q3 | bgx bgy bgz]
%  (quatérnio NED->corpo + viés do giro), com:
%     • PREDIÇÃO     : integra o giro, descontando o viés estimado.
%     • CORREÇÃO 1   : acelerômetro como referência de gravidade  -> ROLL/PITCH.
%     • CORREÇÃO 2   : heading do magnetômetro (tilt-compensado)  -> YAW.
%  É o mesmo princípio do Pixhawk para atitude; falta só o GPS (ver nota abaixo).
%
%  MAG: usa o compass ATIVO (MAG_1 = compass 2 externo; o test_flight.m já coloca
%  esse no logs_concat.mat). O mag deste drone é ruidoso (~20° de espalhamento de
%  heading, |mag| varia ~20%) -> fundido com ruído alto (sig_yaw), igual ao EKF
%  real (EK3_YAW_M_NSE=0.5 rad). Ele LIMITA a deriva do yaw, não a zera.
%
%  GPS: logs_concat.mat agora tem velocidade do GPS, mas fundi-la exige estados de
%  velocidade (um INS-EKF completo de ~15 estados). Neste voo (hover, GPS<2.5 m/s)
%  ela ajudaria pouco o tilt; então este EKF de ATITUDE não a usa (próximo nível).
%
%  Dados:  ../1_data/logs_concat.mat   | Janela: 1 a 110 s
%  Rodar:  >> ekf_atitude

clear; clc; close all;
HERE = fileparts(mfilename('fullpath'));
LOGF = fullfile(HERE, '..', '1_data', 'logs_concat.mat');
W    = [1 110];     % janela [s]
g    = 9.81;

% ---- chaves ----
USE_MAG = true;     % false = só giro+accel (yaw deriva) p/ comparar

% ---- pesos (do PARM real; ajustáveis) ----
sig_gyro  = 0.015;  % EK3_GYRO_P_NSE  [rad/s]
sig_gbias = 0.001;  % EK3_GBIAS_P_NSE [rad/s]
sig_acc   = 0.5;    % ruído-base do accel-como-gravidade [m/s^2]
k_acc_dyn = 3.0;    % infla R do accel ~ | |a|-g |
acc_gate  = 4.0;    % rejeita correção de accel se | |a|-g | > acc_gate
sig_yaw   = 0.35;   % ruído do heading do mag [rad] (~20°, medido) ~ EK3_YAW_M_NSE
mag_dec   = -0.374; % COMPASS_DEC [rad] (~-21.4°, declinação do Brasil)
mag_amp_tol = 0.40; % rejeita mag se | |m|-mediana | > tol*mediana
yaw_gate    = deg2rad(60);   % rejeita inovação de yaw acima disto

%% ---- carrega e recorta ----
S = load(LOGF);  L = S.L;
im = L.time_IMU >= W(1) & L.time_IMU <= W(2);
at = L.time_ATT >= W(1) & L.time_ATT <= W(2);

t   = L.time_IMU(im);  N = numel(t);
gyr = [L.gyrX_raw(im), L.gyrY_raw(im), L.gyrZ_raw(im)];
acc = [L.accX_raw(im), L.accY_raw(im), L.accZ_raw(im)];

t_att    = L.time_ATT(at);
roll_log = L.roll_deg(at);  pitch_log = L.pitch_deg(at);  yaw_log = L.yaw_deg(at);

% MAG interpolado na grade da IMU (NaN fora do range -> gated)
has_mag = isfield(L,'magX') && ~isempty(L.time_MAG);
if USE_MAG && has_mag
    mxi = interp1(L.time_MAG, L.magX, t);
    myi = interp1(L.time_MAG, L.magY, t);
    mzi = interp1(L.time_MAG, L.magZ, t);
    mag_med = median(sqrt(mxi.^2+myi.^2+mzi.^2), 'omitnan');
else
    mxi = nan(N,1); myi = nan(N,1); mzi = nan(N,1); mag_med = NaN;
end
fprintf('Janela [%g, %g] s | IMU: %d | ATT: %d | MAG: %s | |mag|~%.0f mGauss\n', ...
    W(1), W(2), N, numel(t_att), string(USE_MAG && has_mag), mag_med);

%% ---- inicialização (do ATT no início) ----
phi0 = deg2rad(interp1(t_att, roll_log,  t(1), 'linear', 'extrap'));
th0  = deg2rad(interp1(t_att, pitch_log, t(1), 'linear', 'extrap'));
ps0  = deg2rad(interp1(t_att, yaw_log,   t(1), 'linear', 'extrap'));
x = [eul2quat_(phi0, th0, ps0); 0; 0; 0];
P = diag([1e-3*ones(1,4), (0.02)^2*ones(1,3)]);

%% ---- loop do EKF (rate da IMU) ----
EUL = zeros(N,3); BIAS = zeros(N,3);
USED_A = false(N,1); USED_M = false(N,1);
AMG = zeros(N,1); PSIM = nan(N,1);
dt_med = median(diff(t));

for k = 1:N
    if k == 1, dt = dt_med; else, dt = t(k) - t(k-1); end
    if dt <= 0 || dt > 0.5, dt = dt_med; end

    qk = x(1:4);  w = gyr(k,:)' - x(5:7);

    % ===== PREDIÇÃO =====
    dth = w*dt; ang = norm(dth);
    if ang > 1e-8, dq = [cos(ang/2); (dth/ang)*sin(ang/2)]; else, dq = [1; 0.5*dth]; end
    qp = quatmul_(qk, dq);  qp = qp/norm(qp);
    Om = Omega_(w); Xi = Xi_(qk);
    F  = eye(7) + [0.5*Om, -0.5*Xi; zeros(3,4), zeros(3)]*dt;
    Gg = [0.5*Xi; zeros(3,3)];
    Q  = Gg*(sig_gyro^2)*Gg'*dt + blkdiag(zeros(4), (sig_gbias^2)*eye(3))*dt;
    x(1:4) = qp;  P = F*P*F' + Q;  x(1:4) = x(1:4)/norm(x(1:4));

    % ===== CORREÇÃO 1: acelerômetro (gravidade -> roll/pitch) =====
    am = norm(acc(k,:));  AMG(k) = am - g;
    if abs(am - g) <= acc_gate
        q0=x(1);q1=x(2);q2=x(3);q3=x(4);
        h = -g*[2*(q1*q3-q0*q2); 2*(q2*q3+q0*q1); 1-2*(q1^2+q2^2)];
        Hq = [ 2*g*q2, -2*g*q3,  2*g*q0, -2*g*q1;
              -2*g*q1, -2*g*q0, -2*g*q3, -2*g*q2;
                  0,    4*g*q1,  4*g*q2,     0   ];
        H = [Hq, zeros(3,3)];
        s_a = sig_acc + k_acc_dyn*abs(am-g);
        K = (P*H')/(H*P*H' + (s_a^2)*eye(3));
        x = x + K*(acc(k,:)' - h);  x(1:4)=x(1:4)/norm(x(1:4));
        P = (eye(7)-K*H)*P;  P=(P+P')/2;  USED_A(k)=true;
    end

    % ===== CORREÇÃO 2: magnetômetro (heading -> yaw) =====
    if USE_MAG && ~isnan(mxi(k))
        m = [mxi(k); myi(k); mzi(k)];  amp = norm(m);
        q0=x(1);q1=x(2);q2=x(3);q3=x(4);
        % roll/pitch atuais do FILTRO para tilt-compensar
        phi_e = atan2(2*(q0*q1+q2*q3), 1-2*(q1^2+q2^2));
        th_e  = asin(max(min(2*(q0*q2-q1*q3),1),-1));
        cph=cos(phi_e); sph=sin(phi_e); cth=cos(th_e); sth=sin(th_e);
        Xh = m(1)*cth + m(2)*sph*sth + m(3)*cph*sth;
        Yh = m(2)*cph - m(3)*sph;
        psi_meas = atan2(-Yh, Xh) + mag_dec;   PSIM(k) = psi_meas;
        % yaw previsto + inovação (wrap)
        a = 2*(q0*q3+q1*q2);  b = 1-2*(q2^2+q3^2);
        psi_pred = atan2(a, b);
        yk = atan2(sin(psi_meas-psi_pred), cos(psi_meas-psi_pred));
        if abs(amp-mag_med) <= mag_amp_tol*mag_med && abs(yk) <= yaw_gate
            den = a^2 + b^2;
            Hy = [2*b*q3, 2*b*q2, 2*b*q1+4*a*q2, 2*b*q0+4*a*q3, 0,0,0]/den;
            Ky = (P*Hy')/(Hy*P*Hy' + sig_yaw^2);
            x = x + Ky*yk;  x(1:4)=x(1:4)/norm(x(1:4));
            P = (eye(7)-Ky*Hy)*P;  P=(P+P')/2;  USED_M(k)=true;
        end
    end

    BIAS(k,:) = x(5:7)';
    [phi, th, ps] = quat2eul_(x(1:4));  EUL(k,:) = [phi, th, ps];
end

%% ---- métricas (vs ATT) ----
rl = interp1(t_att, roll_log,  t, 'linear', 'extrap');
pl = interp1(t_att, pitch_log, t, 'linear', 'extrap');
yl = interp1(t_att, yaw_log,   t, 'linear', 'extrap');
e_phi = rad2deg(EUL(:,1)) - rl;
e_th  = rad2deg(EUL(:,2)) - pl;
e_ps  = wrap180(rad2deg(EUL(:,3)) - yl);

if USE_MAG, magtag = '+mag'; else, magtag = ''; end
fprintf('\n=== EKF de atitude (giro+accel%s) vs ATT do log ===\n', magtag);
fprintf('  correções: accel %d, mag %d (de %d amostras)\n', sum(USED_A), sum(USED_M), N);
fprintf('  roll  : RMS=%6.2f deg | final=%+7.2f deg\n', rms(e_phi), e_phi(end));
fprintf('  pitch : RMS=%6.2f deg | final=%+7.2f deg\n', rms(e_th),  e_th(end));
if USE_MAG
    fprintf('  yaw   : RMS=%6.2f deg | final=%+7.2f deg   <- agora LIMITADO pelo mag\n', rms(e_ps), e_ps(end));
else
    fprintf('  yaw   : RMS=%6.2f deg | final=%+7.2f deg   <- DERIVA (mag desligado)\n', rms(e_ps), e_ps(end));
end
fprintf('  viés do giro estimado (final): [%+.4f %+.4f %+.4f] rad/s\n', BIAS(end,:));

%% ---- gráfico principal ----
fig = figure('Color','w','Position',[50 50 1280 780]);
nm  = {'\phi (roll)','\theta (pitch)','\psi (yaw)'};
est = {rad2deg(EUL(:,1)), rad2deg(EUL(:,2)), rad2deg(EUL(:,3))};
lg  = {rl, pl, yl};  er = {e_phi, e_th, e_ps};
if USE_MAG, ynote='(corrigido pelo mag)'; else, ynote='(SÓ giro: deriva)'; end
note= {'(corrigido pelo accel)','(corrigido pelo accel)', ynote};
for i = 1:3
    subplot(3,2,2*i-1); hold on; grid on;
    plot(t, lg{i},  'b-',  'LineWidth',1.3, 'DisplayName','ATT (EKF do Pixhawk)');
    plot(t, est{i}, 'r--', 'LineWidth',1.3, 'DisplayName','EKF recriado');
    ylabel([nm{i} ' [°]']);
    if i==1, legend('Location','best'); title('Atitude: EKF recriado vs log'); end
    text(0.02,0.92, note{i}, 'Units','normalized','Color',[.3 .3 .3],'FontSize',8);
    if i==3, xlabel('t [s]'); end
    subplot(3,2,2*i); hold on; grid on;
    plot(t, er{i}, 'm-', 'LineWidth',1); yline(0,'k--');
    ylabel(['\Delta' nm{i} ' [°]']);
    title(sprintf('Erro (RMS=%.2f°, final=%+.2f°)', rms(er{i}), er{i}(end)));
    if i==3, xlabel('t [s]'); end
end
sgtitle(sprintf('EKF de atitude recriado (giro+accel%s) vs ATT do Pixhawk', magtag));
saveas(fig, fullfile(HERE,'ekf_atitude_cmp.png'));

%% ---- gráfico diagnóstico ----
fig2 = figure('Color','w','Position',[80 80 1200 760]);
subplot(2,2,1); hold on; grid on;
plot(t, rad2deg(BIAS(:,1)),'r', t, rad2deg(BIAS(:,2)),'g', t, rad2deg(BIAS(:,3)),'b');
ylabel('viés giro [°/s]'); legend('b_x','b_y','b_z','Location','best');
title('Viés do giroscópio estimado');
subplot(2,2,2); hold on; grid on;
plot(t, AMG, 'k'); yline([acc_gate -acc_gate],'r--'); plot(t(~USED_A), AMG(~USED_A),'r.');
ylabel('|a|-g [m/s^2]'); title('Gating do accel (vermelho = rejeitado)');
subplot(2,2,3); hold on; grid on;
plot(t, mod(yl,360), 'b-', 'LineWidth',1.3, 'DisplayName','ATT yaw');
plot(t, mod(rad2deg(EUL(:,3)),360), 'r--','LineWidth',1.3,'DisplayName','EKF yaw');
plot(t, mod(rad2deg(PSIM),360), '.', 'Color',[.6 .6 .6], 'MarkerSize',4, 'DisplayName','heading mag (cru)');
ylabel('yaw [°]'); xlabel('t [s]'); legend('Location','best');
title('Yaw: mag cru (ruidoso) -> EKF filtra -> casa com ATT');
subplot(2,2,4); hold on; grid on;
plot(t, e_ps, 'm-'); yline(0,'k--');
ylabel('\Delta yaw [°]'); xlabel('t [s]');
title(sprintf('Erro de yaw vs ATT (RMS=%.2f°)', rms(e_ps)));
saveas(fig2, fullfile(HERE,'ekf_atitude_diag.png'));

fprintf('Figuras salvas em %s\n', HERE);

%% ===================== subfunções =====================
function q = eul2quat_(phi, th, ps)
    cr=cos(phi/2); sr=sin(phi/2); cp=cos(th/2); sp=sin(th/2); cy=cos(ps/2); sy=sin(ps/2);
    q = [cr*cp*cy + sr*sp*sy; sr*cp*cy - cr*sp*sy; cr*sp*cy + sr*cp*sy; cr*cp*sy - sr*sp*cy];
    q = q/norm(q);
end
function [phi, th, ps] = quat2eul_(q)
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    phi = atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1^2 + q2^2));
    th  = asin(max(min(2*(q0*q2 - q1*q3), 1), -1));
    ps  = atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2^2 + q3^2));
end
function r = quatmul_(a, b)
    r = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
         a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
         a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
         a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)];
end
function Om = Omega_(w)
    Om = [ 0 -w(1) -w(2) -w(3); w(1) 0 w(3) -w(2); w(2) -w(3) 0 w(1); w(3) w(2) -w(1) 0];
end
function Xi = Xi_(q)
    Xi = [-q(2) -q(3) -q(4); q(1) -q(4) q(3); q(4) q(1) -q(2); -q(3) q(2) q(1)];
end
function a = wrap180(a), a = mod(a + 180, 360) - 180; end
