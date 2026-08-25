% measure_kv.m — medir do voo a derivada de empuxo com a velocidade axial
% =========================================================================
% k_v = −∂T/∂w é o número que domina o amortecimento angular no pairado
% (mecanismo [A] do damping_budget). Hoje ele vem de teoria com hipóteses sobre
% a pá (corda, inclinação de sustentação, arrasto de perfil), e os dois limites
% da teoria são muito distantes:
%     quase-estático (elemento de pá)  k_v = (T0/v_h)·κ ,  κ = aσ/(16λ_i+aσ)
%     influxo congelado                k_v = T0/v_h                (κ = 1)
% Aqui k_v é MEDIDO: para uma dada rotação comandada, o empuxo real cai quando o
% drone sobe. O acelerômetro dá o empuxo (a_z = −T/m·cosφcosθ, força específica)
% e o EKF com barômetro dá a velocidade vertical w. A regressão
%     T_medido − T_modelo(RPM) = −4·k_v·w
% devolve k_v sem nenhuma hipótese sobre a pá.
%
% Uso:  >> measure_kv
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
WINS = [4 125; 400 480; 525 650];        % trechos de voo com subidas e descidas

L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';
ip = @(tt,xx) interp1(tt, xx, tg, 'linear');
W4 = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
[rpm, fT] = motor_chain(tg, W4);
T_mod = sum(fT(rpm), 2);                              % empuxo da bancada, sem escala
az  = ip(L.time_IMU, L.accZ_raw);
phi = deg2rad(ip(L.time_ATT, L.roll_deg));
th  = deg2rad(ip(L.time_ATT, L.pitch_deg));
VE  = estimate_velocity(L, tg);
w   = VE.VD;                                          % positivo para BAIXO (NED)

% empuxo medido pelo acelerômetro: a_z = −(T/m)·cos(eps) + bias
bias_z = p.imu_bias(3);
T_meas = -(az - bias_z) * p.m;                        % [N] ao longo de −z do corpo

sel = isfinite(w) & isfinite(T_meas) & abs(phi) < deg2rad(25) & abs(th) < deg2rad(25) & T_mod > 5;
inw = false(size(tg));
for k = 1:size(WINS,1), inw = inw | (tg>=WINS(k,1) & tg<=WINS(k,2)); end
sel = sel & inw;

% Regressão: T_meas = kT·T_mod − 4·k_v·(−w)   [w do NED é para baixo; subir é w<0]
%   subir (w_up = −w > 0) reduz o empuxo →  T_meas = kT·T_mod − 4·k_v·w_up
w_up = -w;
X = [T_mod(sel), w_up(sel)];  y = T_meas(sel);
b = X\y;
kT_hat = b(1);  kv_hat = -b(2)/4;
res = y - X*b;  R2 = 1 - sum(res.^2)/sum((y-mean(y)).^2);
n = sum(sel);  se = sqrt(diag(inv(X'*X)) * sum(res.^2)/(n-2));

% teoria, para comparar
R = 0.127; A = pi*R^2; Nb = 2; cb = 0.025; a_l = 5.0;
sig = Nb*cb/(pi*R);  T0 = p.m*p.g/4;
rpm_h = interp1(p.bench.T_grams*1e-3*p.g, p.bench.RPM, T0, 'linear');
Om = rpm_h*2*pi/60;  vh = sqrt(T0/(2*p.rho*A));  lam = vh/(Om*R);
kap = a_l*sig/(16*lam + a_l*sig);
kv_qs = (T0/vh)*kap;   kv_fr = T0/vh;

fprintf('\n  Amostras: %d (%.0f s de voo, |φ|,|θ| < 25°)\n', n, n*0.1);
fprintf('  Faixa de w para cima: %.2f a %.2f m/s (p95 |w| = %.2f)\n', min(w_up(sel)), max(w_up(sel)), prctile(abs(w_up(sel)),95));
fprintf('\n  Regressão  T_medido = k_T·T_bancada − 4·k_v·w_subida     (R² = %.3f)\n', R2);
fprintf('    k_T = %.3f ± %.3f\n', kT_hat, se(1));
fprintf('    k_v = %.3f ± %.3f N/(m/s) por rotor\n', kv_hat, se(2)/4);
fprintf('\n  Teoria:  quase-estático %.3f (κ = %.3f)  |  influxo congelado %.3f  |  T0/v_h = %.3f\n', ...
    kv_qs, kap, kv_fr, T0/vh);
fprintf('  κ implícito na medida: %.3f\n', kv_hat/(T0/vh));
ly = mean([p.arms.Lx_r p.arms.Lx_l]);  lx = mean([p.arms.Ly_f p.arms.Ly_r]);
fprintf('\n  → c_p do influxo com k_v medido: %.2f 1/s   (teoria quase-estática %.2f)\n', ...
    4*kv_hat*ly^2/p.J.Jx, 4*kv_qs*ly^2/p.J.Jx);
fprintf('  → c_q do influxo com k_v medido: %.2f 1/s   (teoria quase-estática %.2f)\n', ...
    4*kv_hat*lx^2/p.J.Jy, 4*kv_qs*lx^2/p.J.Jy);
save(fullfile(paths.outputs,'measure_kv.mat'), 'kv_hat','kT_hat','se','R2','kv_qs','kv_fr','n');

%% figura
f = figure('Position',[40 40 1150 450],'Color','w'); try, f.Theme='light'; catch, end
subplot(1,2,1); hold on; grid on;
resid_T = y - kT_hat*X(:,1);                 % o que sobra depois de tirar o empuxo comandado
scatter(w_up(sel), resid_T, 6, 'filled', 'MarkerFaceAlpha',0.15, 'MarkerFaceColor',[0 0.45 0.7]);
ww = linspace(min(w_up(sel)), max(w_up(sel)), 50)';
plot(ww, -4*kv_hat*ww, 'r-', 'LineWidth',2.2);
plot(ww, -4*kv_qs*ww, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.6);
plot(ww, -4*kv_fr*ww, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.6);
xlabel('velocidade de subida w [m/s]'); ylabel('empuxo medido − comandado [N]');
legend({'dado', sprintf('ajuste: k_v = %.2f', kv_hat), sprintf('quase-estático %.2f', kv_qs), ...
        sprintf('influxo congelado %.2f', kv_fr)}, 'Location','southwest','FontSize',8);
title('derivada de empuxo com a velocidade axial');

subplot(1,2,2); hold on; grid on;
[bn, ed] = discretize(w_up(sel), -3:0.5:3);
mu = accumarray(bn(~isnan(bn)), resid_T(~isnan(bn)), [], @mean, NaN);
sd = accumarray(bn(~isnan(bn)), resid_T(~isnan(bn)), [], @std, NaN);
ctr = (ed(1:end-1)+ed(2:end))/2;
errorbar(ctr(1:numel(mu)), mu, sd, 'ko-', 'LineWidth',1.4, 'MarkerFaceColor','k');
plot(ww, -4*kv_hat*ww, 'r-', 'LineWidth',2);
xlabel('velocidade de subida w [m/s]'); ylabel('empuxo medido − comandado [N]');
title('médias por faixa de w (barras: desvio)');
fn = fullfile(paths.images,'measure_kv.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',140);
fprintf('  Figura: %s\n', fn);
