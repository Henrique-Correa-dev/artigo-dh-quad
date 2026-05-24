%% test_dynamics.m — Teste rapido do modelo vtol_dynamics.m
clear; close all; clc;

%% Carregar dados
load(fullfile('..','identification','log_data.mat'));
IMU.TimeS = double(IMU.TimeUS)/1e6;
RCOU.TimeS = double(RCOU.TimeUS)/1e6;
ATT.TimeS = double(ATT.TimeUS)/1e6;

idx = IMU.I==0; dt=0.1;
ts = max(min(IMU.TimeS(idx)),min(RCOU.TimeS));
te = min(max(IMU.TimeS(idx)),max(RCOU.TimeS));
t = (ts:dt:te)';

pwm = [interp1(RCOU.TimeS,double(RCOU.C1),t,'linear'), ...
       interp1(RCOU.TimeS,double(RCOU.C2),t,'linear'), ...
       interp1(RCOU.TimeS,double(RCOU.C3),t,'linear'), ...
       interp1(RCOU.TimeS,double(RCOU.C4),t,'linear')];

p_m = interp1(IMU.TimeS(idx),IMU.GyrX(idx),t,'linear');
q_m = interp1(IMU.TimeS(idx),IMU.GyrY(idx),t,'linear');
r_m = interp1(IMU.TimeS(idx),IMU.GyrZ(idx),t,'linear');
roll_m = interp1(ATT.TimeS,ATT.Roll,t,'linear');
pitch_m = interp1(ATT.TimeS,ATT.Pitch,t,'linear');
yaw_m = interp1(ATT.TimeS,ATT.Yaw,t,'linear');

%% Setup parametros (CAD + bancada)
params = setup_params();

%% Condicao inicial em t=150s
t0 = 150; tf = 155;
i0 = find(t>=t0,1);
x0 = [p_m(i0); q_m(i0); r_m(i0); ...
      deg2rad(roll_m(i0)); deg2rad(pitch_m(i0)); deg2rad(yaw_m(i0)); ...
      0; 0; 0];

fprintf('=== TESTE DO MODELO vtol_dynamics.m ===\n');
fprintf('  Integracao: %.1f a %.1f s\n', t0, tf);
fprintf('  CI: p=%.3f q=%.3f r=%.3f phi=%.1f° theta=%.1f° psi=%.1f°\n', ...
    x0(1:3), rad2deg(x0(4:6)));

%% Integrar
idx_sim = (t>=t0) & (t<=tf);
[t_sim, x_sim] = ode45(@(tt,xx) vtol_dynamics(tt,xx,t,pwm,params), t(idx_sim), x0);

fprintf('  %d passos do ode45\n', length(t_sim));
fprintf('  Estados finais:\n');
fprintf('    p=%.4f q=%.4f r=%.4f [rad/s]\n', x_sim(end,1:3));
fprintf('    phi=%.2f° theta=%.2f° psi=%.2f° [deg]\n', rad2deg(x_sim(end,4:6)));
fprintf('    u=%.2f v=%.2f w=%.2f [m/s]\n', x_sim(end,7:9));

%% Comparar com medido
t_comp = t(idx_sim);
p_sim = interp1(t_sim, x_sim(:,1), t_comp, 'linear');
q_sim = interp1(t_sim, x_sim(:,2), t_comp, 'linear');
r_sim = interp1(t_sim, x_sim(:,3), t_comp, 'linear');

R2 = @(meas,pred) 1 - sum((meas-pred).^2)/sum((meas-mean(meas)).^2);
fprintf('\n  R2 (janela de 5s, parametros do CAD, k=1, sem amortecimento):\n');
fprintf('    p: %.3f\n', R2(p_m(idx_sim), p_sim));
fprintf('    q: %.3f\n', R2(q_m(idx_sim), q_sim));
fprintf('    r: %.3f\n', R2(r_m(idx_sim), r_sim));

%% Plot
figure('Position',[50 50 1000 600],'Color','w');

subplot(3,1,1);
plot(t_comp, p_m(idx_sim), 'b-', t_comp, p_sim, 'r--', 'LineWidth', 1.2);
ylabel('p [rad/s]'); grid on; legend('Medido','Simulado');
title(sprintf('Teste vtol\\_dynamics.m  (%.0f-%.0fs, params CAD)', t0, tf));

subplot(3,1,2);
plot(t_comp, q_m(idx_sim), 'b-', t_comp, q_sim, 'r--', 'LineWidth', 1.2);
ylabel('q [rad/s]'); grid on; legend('Medido','Simulado');

subplot(3,1,3);
plot(t_comp, r_m(idx_sim), 'b-', t_comp, r_sim, 'r--', 'LineWidth', 1.2);
ylabel('r [rad/s]'); xlabel('Tempo [s]'); grid on; legend('Medido','Simulado');

saveas(gcf, 'fig11_test_dynamics.png');
fprintf('\n  Modelo OK! Figura salva.\n');
