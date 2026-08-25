%FILTRO_COMPLEMENTAR  Estimação de atitude (roll/pitch): accel puro vs giro vs
%                     filtro complementar, comparado com o ATT do log.
%
%  AUTOCONTIDO — não depende do pipeline de identificação. Só usa o loader de
%  log (load_log_data) e setup_paths. Nada aqui altera o modelo/identificação.
%
%  Teoria (Quan 2017, "Introduction to Multicopter Design and Control", §9.1.1):
%    - Acelerômetro mede força específica. Em baixa dinâmica (Eq. 9.1 com termos
%      de velocidade/arrasto desprezados):  a_xb ≈  g·sinθ ,  a_yb ≈ -g·cosθ·sinφ.
%    - Logo (Eq. 9.2):  θ_acc = asin(a_xb/g) ,  φ_acc = -asin(a_yb/(g·cosθ_acc)).
%      → estável (sem deriva), mas RUIDOSO e ruim em manobra (a ≠ g).
%    - Giroscópio: integra a cinemática de Euler → suave e bom em manobra, mas
%      DERIVA (acumula bias).
%    - Filtro complementar: passa-alta no giro + passa-baixa no accel.
%        x_cf = α·(x_cf + ẋ_giro·dt) + (1-α)·x_acc ,  α = τ/(τ+dt).
%
%  Saída: figura ..._cmp.png + RMSE de cada método vs ATT.

clear; clc; close all;

% ---------------------------------------------------------------- config
HERE = fileparts(mfilename('fullpath'));
addpath(fullfile(HERE, '..'));        % raiz → setup_paths
setup_paths();                        % adiciona 2_model (load_log_data) etc.

LOG_FILE = 'logs_concat.mat';
t_window = [1, 110];                  % [s] janela de análise
g        = 9.80665;                   % [m/s^2]
TAU      = 1.0;                        % [s] constante de tempo do compl. filter
                                      %  (maior = confia mais no giro; menor = no accel)

% ---------------------------------------------------------------- dados
paths = setup_paths();
L = load_log_data(fullfile(paths.data, LOG_FILE));

im = (L.time_IMU >= t_window(1)) & (L.time_IMU <= t_window(2));
t   = L.time_IMU(im);
acc = [L.accX_raw(im), L.accY_raw(im), L.accZ_raw(im)];   % m/s^2 (força específica, FRD)
gyr = [L.gyrX_raw(im), L.gyrY_raw(im), L.gyrZ_raw(im)];   % rad/s  (p,q,r)
N   = numel(t);
dt  = [diff(t); median(diff(t))];     % passo por amostra

% Referência: ATT do log (graus → rad), reamostrado na grade do IMU
phi_ref = deg2rad(interp1(L.time_ATT, L.roll_deg,  t, 'linear', 'extrap'));
th_ref  = deg2rad(interp1(L.time_ATT, L.pitch_deg, t, 'linear', 'extrap'));

% ---------------------------------------------------------------- (1) ACCEL PURO  (Quan Eq. 9.2)
clamp = @(x) max(min(x, 1), -1);
th_acc  = asin( clamp(acc(:,1) ./ g) );
phi_acc = -asin( clamp(acc(:,2) ./ (g .* cos(th_acc))) );

% ---------------------------------------------------------------- (2) GIRO PURO (integra Euler → DERIVA)
phi_gyro = zeros(N,1);  th_gyro = zeros(N,1);
phi_gyro(1) = phi_ref(1);  th_gyro(1) = th_ref(1);     % IC = ATT inicial
for k = 1:N-1
    p = gyr(k,1); q = gyr(k,2); r = gyr(k,3);
    ph = phi_gyro(k); te = th_gyro(k);
    phi_dot = p + sin(ph)*tan(te)*q + cos(ph)*tan(te)*r;
    th_dot  =     cos(ph)*q        - sin(ph)*r;
    phi_gyro(k+1) = ph + phi_dot*dt(k);
    th_gyro(k+1)  = te + th_dot *dt(k);
end

% ---------------------------------------------------------------- (3) FILTRO COMPLEMENTAR
phi_cf = zeros(N,1);  th_cf = zeros(N,1);
phi_cf(1) = phi_ref(1);  th_cf(1) = th_ref(1);
for k = 1:N-1
    p = gyr(k,1); q = gyr(k,2); r = gyr(k,3);
    ph = phi_cf(k); te = th_cf(k);
    % predição (giro, cinemática de Euler)
    phi_dot = p + sin(ph)*tan(te)*q + cos(ph)*tan(te)*r;
    th_dot  =     cos(ph)*q        - sin(ph)*r;
    phi_pred = ph + phi_dot*dt(k);
    th_pred  = te + th_dot *dt(k);
    % correção (accel, passa-baixa) — α depende do dt local
    a = TAU / (TAU + dt(k));
    phi_cf(k+1) = a*phi_pred + (1-a)*phi_acc(k+1);
    th_cf(k+1)  = a*th_pred  + (1-a)*th_acc(k+1);
end

% ---------------------------------------------------------------- métricas (RMSE em graus vs ATT)
rmse = @(e) sqrt(mean(e.^2, 'omitnan'));
r2d  = @(x) x*180/pi;
E = struct();
E.acc_phi  = rmse(r2d(phi_acc  - phi_ref));   E.acc_th  = rmse(r2d(th_acc  - th_ref));
E.gyro_phi = rmse(r2d(phi_gyro - phi_ref));   E.gyro_th = rmse(r2d(th_gyro - th_ref));
E.cf_phi   = rmse(r2d(phi_cf   - phi_ref));   E.cf_th   = rmse(r2d(th_cf   - th_ref));

fprintf('\n=== RMSE vs ATT do log  (janela %g–%g s, %d amostras, TAU=%.2fs) ===\n', ...
        t_window(1), t_window(2), N, TAU);
fprintf('              ROLL (phi)    PITCH (theta)\n');
fprintf('  Accel puro   %6.2f deg     %6.2f deg\n', E.acc_phi,  E.acc_th);
fprintf('  Giro puro    %6.2f deg     %6.2f deg   <- deriva\n', E.gyro_phi, E.gyro_th);
fprintf('  Compl.filter %6.2f deg     %6.2f deg   <- melhor\n', E.cf_phi,  E.cf_th);

% ---------------------------------------------------------------- plot
figure('Position', [80 80 1180 720], 'Color', 'w');
tt = t - t(1);

subplot(2,1,1); hold on; grid on;
plot(tt, r2d(th_acc),  'Color',[.85 .55 .55], 'LineWidth',.6);
plot(tt, r2d(th_gyro), 'Color',[.50 .50 .50], 'LineWidth',1.0);
plot(tt, r2d(th_cf),   'r-', 'LineWidth',1.4);
plot(tt, r2d(th_ref),  'Color',[0 0.35 0.95], 'LineWidth',1.2);
ylabel('\theta  pitch [deg]');
title(sprintf('Estimação de atitude: accel puro vs giro vs filtro complementar (\\tau=%.1fs)', TAU));
legend({sprintf('accel (%.2f°)',E.acc_th), sprintf('giro (%.2f°)',E.gyro_th), ...
        sprintf('compl. (%.2f°)',E.cf_th), 'ATT (log)'}, 'Location','best');

subplot(2,1,2); hold on; grid on;
plot(tt, r2d(phi_acc),  'Color',[.85 .55 .55], 'LineWidth',.6);
plot(tt, r2d(phi_gyro), 'Color',[.50 .50 .50], 'LineWidth',1.0);
plot(tt, r2d(phi_cf),   'r-', 'LineWidth',1.4);
plot(tt, r2d(phi_ref),  'Color',[0 0.35 0.95], 'LineWidth',1.2);
ylabel('\phi  roll [deg]'); xlabel('t [s]');
legend({sprintf('accel (%.2f°)',E.acc_phi), sprintf('giro (%.2f°)',E.gyro_phi), ...
        sprintf('compl. (%.2f°)',E.cf_phi), 'ATT (log)'}, 'Location','best');

out_png = fullfile(HERE, 'filtro_complementar_cmp.png');
exportgraphics(gcf, out_png, 'Resolution', 130);
fprintf('\nFigura salva: %s\n', out_png);
