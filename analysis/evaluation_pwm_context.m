%% evaluation_pwm_context.m
% =========================================================================
% Avaliacao dos 4 PWMs no contexto da dinamica do drone:
%   - PWMs (entrada do sistema)
%   - Empuxo/torque estimado por motor (com k_T identificado)
%   - Gyro (p, q, r)
%   - Atitude e setpoint do autopilot (DesRoll vs Roll, etc)
%   - Tracking error
%   - GPS speed (proxy de vento)
%
% Objetivo: entender PORQUE os PWMs sao tao assimetricos no hover
% (motor weakness vs wind vs aerodinamica vs CG)
% =========================================================================

clear; close all; clc;
cd(fileparts(mfilename('fullpath')));

%% 1. Carregar dados
load(fullfile('..','identification','log_data.mat'));

g_acc = 9.80665;
pwm_bench = [1000; 1200; 1400; 1600; 1800; 2000];
T_bench   = [0; 1.4022; 3.2154; 5.2163; 7.6884; 8.2664];
Q_bench   = [0; 0.034; 0.070; 0.115; 0.171; 0.176];

% Conversoes de tempo
ATT.TimeS  = double(ATT.TimeUS) / 1e6;
IMU.TimeS  = double(IMU.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS)/ 1e6;
GPS.TimeS  = double(GPS.TimeUS) / 1e6;

idx_imu = IMU.I == 0;
time_IMU  = IMU.TimeS(idx_imu);
gyrX = IMU.GyrX(idx_imu);
gyrY = IMU.GyrY(idx_imu);
gyrZ = IMU.GyrZ(idx_imu);

idx_gps = GPS.I == 0;
time_GPS = GPS.TimeS(idx_gps);
gps_spd  = GPS.Spd(idx_gps);
gps_crs  = GPS.GCrs(idx_gps);
gps_vz   = GPS.VZ(idx_gps);

%% 2. Janela de analise
t_ini = 150; t_fim = 200;
dt = 0.1;
t = (t_ini:dt:t_fim)';

% Reamostrar
p = interp1(time_IMU, gyrX, t, 'linear');
q = interp1(time_IMU, gyrY, t, 'linear');
r = interp1(time_IMU, gyrZ, t, 'linear');

Roll      = interp1(ATT.TimeS, ATT.Roll,    t, 'linear');
Pitch     = interp1(ATT.TimeS, ATT.Pitch,   t, 'linear');
Yaw       = interp1(ATT.TimeS, ATT.Yaw,     t, 'linear');
DesRoll   = interp1(ATT.TimeS, ATT.DesRoll, t, 'linear');
DesPitch  = interp1(ATT.TimeS, ATT.DesPitch,t, 'linear');
DesYaw    = interp1(ATT.TimeS, ATT.DesYaw,  t, 'linear');
ErrRP     = interp1(ATT.TimeS, ATT.ErrRP,   t, 'linear');
ErrYaw    = interp1(ATT.TimeS, ATT.ErrYaw,  t, 'linear');

pwm1 = interp1(RCOU.TimeS, double(RCOU.C1), t, 'linear');
pwm2 = interp1(RCOU.TimeS, double(RCOU.C2), t, 'linear');
pwm3 = interp1(RCOU.TimeS, double(RCOU.C3), t, 'linear');
pwm4 = interp1(RCOU.TimeS, double(RCOU.C4), t, 'linear');

gps_spd_i = interp1(time_GPS, gps_spd, t, 'linear', 'extrap');
gps_vz_i  = interp1(time_GPS, gps_vz,  t, 'linear', 'extrap');

%% 3. Computar empuxo/torque por motor (com modelo de bancada, k=1)
T_ref_fn = @(pwm) max(0, interp1(pwm_bench, T_bench, ...
    min(max(pwm,1000),2000), 'makima', 'extrap'));
Q_ref_fn = @(pwm) max(0, interp1(pwm_bench, Q_bench, ...
    min(max(pwm,1000),2000), 'makima', 'extrap'));

T_per_motor = [T_ref_fn(pwm1), T_ref_fn(pwm2), T_ref_fn(pwm3), T_ref_fn(pwm4)];
Q_per_motor = [Q_ref_fn(pwm1), Q_ref_fn(pwm2), Q_ref_fn(pwm3), Q_ref_fn(pwm4)];

T_total = sum(T_per_motor, 2);

%% 4. Estatisticas
fprintf('\n========================================================\n');
fprintf(' ASSIMETRIA DE PWM NO HOVER (janela %d-%ds)\n', t_ini, t_fim);
fprintf('========================================================\n');
fprintf('  Motor   Posicao             PWM medio   PWM std    %% saturado(>1900)\n');
labels = {'C1', 'C2', 'C3', 'C4'};
posic  = {'Front-Dir', 'Rear-Esq', 'Front-Esq', 'Rear-Dir'};
pwm_all = [pwm1, pwm2, pwm3, pwm4];
for k = 1:4
    pwm_k = pwm_all(:,k);
    fprintf('  %-6s  %-18s %8.1f    %7.1f    %5.1f%%\n', ...
        labels{k}, posic{k}, mean(pwm_k), std(pwm_k), ...
        100*sum(pwm_k>1900)/numel(pwm_k));
end

fprintf('\n  T_total medio = %.2f N  (peso ~ %.2f N) -> %.1f%% do peso\n', ...
    mean(T_total), 1.6011*9.81, 100*mean(T_total)/(1.6011*9.81));
fprintf('  GPS speed media = %.2f m/s (>>0 indica vento ou translacao)\n', mean(gps_spd_i));
fprintf('  GPS Vz media    = %.2f m/s (>0 = subindo)\n', mean(gps_vz_i));

% Tracking error do autopilot
fprintf('\n  Erro de tracking (Des - Real):\n');
fprintf('    Roll:  media=%+.2f  std=%.2f  deg\n', mean(DesRoll-Roll), std(DesRoll-Roll));
fprintf('    Pitch: media=%+.2f  std=%.2f  deg\n', mean(DesPitch-Pitch), std(DesPitch-Pitch));
fprintf('    Yaw:   media=%+.2f  std=%.2f  deg\n', mean(DesYaw-Yaw), std(DesYaw-Yaw));
fprintf('  ErrRP medio (estimativa do autopilot) = %.4f\n', mean(ErrRP));

%% 5. Diagnostico de assimetria
fprintf('\n========================================================\n');
fprintf(' DIAGNOSTICO DE ASSIMETRIA\n');
fprintf('========================================================\n');

% Mxyz medios (com geometria conhecida)
lx_f = 0.31118; lx_r = 0.34287; ly = 0.232;
% Sinais por motor (C1 F-Dir, C2 R-Esq, C3 F-Esq, C4 R-Dir)
y_sign = [-1, +1, +1, -1];  % Mx contribuicao = -y_i * T_i
x_sign = [+1, -1, +1, -1];  % My contribuicao = +x_i * T_i (com sinal de frente/tras)
d_sign = [+1, +1, -1, -1];  % Mz contribuicao = +d_i * Q_i (CW/CCW)
x_arm  = [lx_f, lx_r, lx_f, lx_r];

Mx = ly * sum(T_per_motor .* y_sign, 2);
My = sum(T_per_motor .* (x_sign .* x_arm), 2);
Mz = sum(Q_per_motor .* d_sign, 2);

fprintf('  Momentos medios (com bancada k=1):\n');
fprintf('    Mx (roll)  = %+.4f Nm  (>0 = drone tende a rolar pra direita)\n', mean(Mx));
fprintf('    My (pitch) = %+.4f Nm  (>0 = drone tende a pitchar pra cima)\n', mean(My));
fprintf('    Mz (yaw)   = %+.4f Nm  (>0 = drone tende a girar CW visto de cima)\n', mean(Mz));

% Se Mx≠0 em media, autopilot compensa algo (vento ou motor weak)
% Se k_T fosse o problema, multiplicaria T_per_motor por k_T -> Mx mudaria

fprintf('\n  Interpretacao:\n');
if abs(mean(My)) > 0.2
    fprintf('    [!!] My medio = %+.3f Nm e GRANDE -- frame esta gerando pitch.\n', mean(My));
    fprintf('         Origem possivel: aerodinamica da asa (este eh VTOL?),\n');
    fprintf('         vento de cauda, CG deslocado em x.\n');
end
if abs(mean(Mx)) > 0.1
    fprintf('    [!]  Mx medio = %+.3f Nm -- assimetria lateral.\n', mean(Mx));
    fprintf('         Origem possivel: motor mais forte de um lado,\n');
    fprintf('         vento lateral, CG deslocado em y.\n');
end
if abs(mean(Mz)) > 0.05
    fprintf('    [!]  Mz medio = %+.3f Nm -- assimetria de yaw.\n', mean(Mz));
    fprintf('         Origem possivel: motor com Q_q diferente,\n');
    fprintf('         vento gerando yaw.\n');
end

%% 6. Figuras
fig = figure('Position',[80 80 1500 900], 'Color','w');
sgtitle(sprintf('Evaluation: PWMs e contexto -- janela %d-%ds', t_ini, t_fim), ...
    'FontWeight','bold');

% --- 1) PWMs dos 4 motores ---
ax1 = subplot(5,1,1);
plot(t, pwm1, 'LineWidth',1.4); hold on;
plot(t, pwm2, 'LineWidth',1.4);
plot(t, pwm3, 'LineWidth',1.4);
plot(t, pwm4, 'LineWidth',1.4);
ylabel('PWM [us]'); grid on; ylim([1000, 2000]);
legend('C1 (F-Dir)','C2 (R-Esq)','C3 (F-Esq)','C4 (R-Dir)','Location','best');
title('Entradas dos 4 motores');

% --- 2) Empuxo por motor ---
ax2 = subplot(5,1,2);
plot(t, T_per_motor(:,1), 'LineWidth',1.2); hold on;
plot(t, T_per_motor(:,2), 'LineWidth',1.2);
plot(t, T_per_motor(:,3), 'LineWidth',1.2);
plot(t, T_per_motor(:,4), 'LineWidth',1.2);
yline(1.6011*9.81/4, 'k--', 'm·g/4 (ideal)', 'LineWidth',1.2);
ylabel('T_i [N]'); grid on;
title('Empuxo por motor (bancada, k=1)');

% --- 3) Setpoint vs medido (Roll e Pitch) ---
ax3 = subplot(5,1,3);
yyaxis left
plot(t, Roll,    'b-',  'LineWidth',1.4); hold on;
plot(t, DesRoll, 'b--', 'LineWidth',1.0);
ylabel('Roll [deg]', 'Color','b'); set(gca,'YColor','b');
yyaxis right
plot(t, Pitch,    'r-',  'LineWidth',1.4);
plot(t, DesPitch, 'r--', 'LineWidth',1.0);
ylabel('Pitch [deg]', 'Color','r'); set(gca,'YColor','r');
grid on;
legend('Roll real','Roll des','Pitch real','Pitch des','Location','best');
title('Atitude: setpoint vs medido (autopilot ja esta corrigindo erros)');

% --- 4) Gyro p, q, r ---
ax4 = subplot(5,1,4);
plot(t, p, 'LineWidth',1.2); hold on;
plot(t, q, 'LineWidth',1.2);
plot(t, r, 'LineWidth',1.2);
ylabel('[rad/s]'); grid on;
legend('p (roll)','q (pitch)','r (yaw)','Location','best');
title('Velocidades angulares (corpo)');

% --- 5) GPS speed (proxy de vento) e Vz ---
ax5 = subplot(5,1,5);
yyaxis left
plot(t, gps_spd_i, 'g-', 'LineWidth',1.4);
ylabel('Vel. horizontal GPS [m/s]', 'Color',[0 0.6 0]);
set(gca,'YColor',[0 0.6 0]);
yyaxis right
plot(t, gps_vz_i, 'm-', 'LineWidth',1.4); hold on;
yline(0, 'k--', 'LineWidth',0.8);
ylabel('Vz GPS [m/s]', 'Color','m');
xlabel('Tempo [s]'); grid on;
title('GPS: velocidade horizontal (vento/translacao) e vertical');

linkaxes([ax1 ax2 ax3 ax4 ax5], 'x');

out_dir = fullfile(pwd, 'images');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
saveas(fig, fullfile(out_dir, 'eval_pwm_context.png'));
fprintf('\nFigura salva: %s\n', fullfile(out_dir, 'eval_pwm_context.png'));

%% 7. Imprimir resumo final
fprintf('\n========================================================\n');
fprintf(' SINTOMAS NO HOVER:\n');
fprintf('========================================================\n');
fprintf('  PWMs:  C1=%.0f  C2=%.0f  C3=%.0f  C4=%.0f us (medio)\n', ...
    mean(pwm1), mean(pwm2), mean(pwm3), mean(pwm4));
fprintf('  Razao PWM max/min = %.2f  (ideal seria 1.0)\n', ...
    max([mean(pwm1) mean(pwm2) mean(pwm3) mean(pwm4)]) / ...
    min([mean(pwm1) mean(pwm2) mean(pwm3) mean(pwm4)]));
fprintf('  Empuxo total / peso = %.2f  (>1 = subindo; <1 = descendo)\n', ...
    mean(T_total) / (1.6011*9.81));
fprintf('  GPS speed medio = %.2f m/s (se for ~hover, esperado <1 m/s)\n', mean(gps_spd_i));
