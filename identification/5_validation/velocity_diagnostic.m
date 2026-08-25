% velocity_diagnostic.m — Tem AIRSPEED suficiente pra a ASA importar?
% =========================================================================
% O resíduo de PITCH (q) tem oscilação 0.6 Hz que p não tem → suspeita de
% momento de arfagem da ASA, M_y = ½ρV²·S·c̄·(Cm0+Cmα·α). MAS isso exige
% velocidade longitudinal e α — que não medimos direto. Este script:
%   1. lê a velocidade NED do GPS (logs_concat.mat: gps_vn/ve/vd) — só <268 s;
%   2. rotaciona NED→corpo pela atitude medida → u, v, w (assume vento ≈ 0);
%   3. calcula V=|vel|, α=atan2(w,u), β=atan2(v,V), q̄=½ρV²;
%   4. estima o momento de asa disponível (q̄·S·c̄) vs o momento dos motores (~Ly·ΔT);
%   5. marca as janelas de treino.
%
% DECISÃO:
%   • V alto (>~5 m/s) e M_y_wing ~ comparável aos momentos de motor → é ASA:
%     vale modelar M_y(α), mas identificar/validar SÓ em janelas com GPS.
%   • V baixo (quase-hover) → asa ∝V² desprezível → o resíduo de q é
%     ROTOR/INFLOW (não precisa de airspeed); modela-se outra coisa.
%
% ⚠️ GPS dá velocidade de SOLO; airspeed = solo − vento. Sem dado de vento,
%    assumimos vento ≈ 0 (limitação). Sinal de VD a conferir (afeta só α).
%
% Uso:  >> velocity_diagnostic

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---- CONFIG ----
LOG_FILE = 'logs_concat.mat';
t_trains = {[4, 24]; [42, 62]; [100, 125]};   % janelas de treino (têm GPS)
rho      = 1.225;                              % densidade do ar (kg/m³, nível do mar)
proj     = parameters();
S  = proj.wing.S;   cbar = proj.wing.c;        % asa (m², m)
Ly = 0.5*(proj.arms.Ly_f + proj.arms.Ly_r);    % braço médio de pitch (m)

%% ---- carregar log (direto: load_log_data não repassa a velocidade GPS) ----
d = load(fullfile(paths.data, LOG_FILE));  L = d.L;
if ~isfield(L,'gps_vn') || isempty(L.gps_vn)
    error('velocity_diagnostic:noGPS', 'logs_concat.mat não tem velocidade GPS (gps_vn).');
end
tG = L.time_GPS(:);
vn = L.gps_vn(:);  ve = L.gps_ve(:);  vd = L.gps_vd(:);

% atitude (graus) interpolada nos instantes do GPS
roll  = deg2rad(interp1(L.time_ATT, L.roll_deg,  tG, 'linear', 'extrap'));
pitch = deg2rad(interp1(L.time_ATT, L.pitch_deg, tG, 'linear', 'extrap'));
yaw   = deg2rad(interp1(L.time_ATT, L.yaw_deg,   tG, 'linear', 'extrap'));

%% ---- rotaciona velocidade NED → corpo (DCM 3-2-1) ----
N = numel(tG);  u=zeros(N,1); v=zeros(N,1); w=zeros(N,1);
for k = 1:N
    cphi=cos(roll(k)); sphi=sin(roll(k));
    cth =cos(pitch(k)); sth=sin(pitch(k));
    cps =cos(yaw(k));  sps=sin(yaw(k));
    % R_body<-NED (aeronáutico, yaw-pitch-roll)
    R = [ cth*cps,                cth*sps,               -sth;
          sphi*sth*cps-cphi*sps,  sphi*sth*sps+cphi*cps,  sphi*cth;
          cphi*sth*cps+sphi*sps,  cphi*sth*sps-sphi*cps,  cphi*cth];
    uvw = R * [vn(k); ve(k); vd(k)];
    u(k)=uvw(1); v(k)=uvw(2); w(k)=uvw(3);
end

V     = sqrt(u.^2 + v.^2 + w.^2);
Vh    = sqrt(vn.^2 + ve.^2);            % velocidade horizontal de solo
alpha = atan2(w, u) * 180/pi;           % ângulo de ataque (graus)
beta  = atan2(v, max(V,1e-3)) * 180/pi; % derrapagem (graus)
qbar  = 0.5 * rho * V.^2;               % pressão dinâmica (Pa)
My_wing_perCm = qbar * S * cbar;        % momento de asa por unidade de Cm (N·m)

%% ---- resumo ----
fprintf('========================================================================\n');
fprintf('  VELOCIDADE / AIRSPEED — trecho com GPS [%.0f-%.0f s], %d fixes\n', tG(1), tG(end), N);
fprintf('========================================================================\n');
fprintf('  V (corpo):   média %.2f  máx %.2f m/s\n', mean(V), max(V));
fprintf('  Vh (solo):   média %.2f  máx %.2f m/s\n', mean(Vh), max(Vh));
fprintf('  α:           [%.1f, %.1f]°   |  β: [%.1f, %.1f]°\n', min(alpha),max(alpha), min(beta),max(beta));
fprintf('  q̄=½ρV²:      média %.2f  máx %.2f Pa\n', mean(qbar), max(qbar));
fprintf('  M_y asa/Cm:  média %.3f  máx %.3f N·m  (×Cm~0.1-0.3 = momento real)\n', ...
    mean(My_wing_perCm), max(My_wing_perCm));
fprintf('  (compare: momento de pitch dos motores ~ Ly·ΔT ≈ %.2f-1 N·m p/ ΔT~3 N, Ly=%.2f)\n', Ly*3, Ly);

fprintf('\n  --- Por janela de treino ---\n');
fprintf('  %-14s %8s %8s %10s\n','janela','V méd','V máx','q̄ máx');
for s = 1:numel(t_trains)
    m = (tG>=t_trains{s}(1)) & (tG<=t_trains{s}(2));
    if any(m)
        fprintf('  [%4.0f,%4.0f]s    %7.2f  %7.2f  %9.1f\n', ...
            t_trains{s}(1), t_trains{s}(2), mean(V(m)), max(V(m)), max(qbar(m)));
    else
        fprintf('  [%4.0f,%4.0f]s    --- sem GPS ---\n', t_trains{s}(1), t_trains{s}(2));
    end
end

if max(V) < 3
    fprintf('\n  >>> VEREDITO: V baixo (quase-hover). Asa ∝V² DESPREZÍVEL.\n');
    fprintf('      O resíduo de q NÃO é asa → é ROTOR/INFLOW (não precisa de airspeed).\n');
else
    fprintf('\n  >>> VEREDITO: há avanço (V máx %.1f m/s). Asa PODE importar.\n', max(V));
    fprintf('      Modele M_y(α) e identifique/valide SÓ em janelas com GPS (<268 s).\n');
end

%% ---- plots ----
fig = figure('Name','velocity_diagnostic','Position',[50 50 1200 850]);
mark = @() arrayfun(@(s) xline(t_trains{s}(1),'g--'), 1:numel(t_trains));
mark2= @() arrayfun(@(s) xline(t_trains{s}(2),'g--'), 1:numel(t_trains));

subplot(4,1,1); hold on; grid on;
plot(tG, V, 'b', 'LineWidth',1.3, 'DisplayName','|V| (corpo)');
plot(tG, Vh,'c--','DisplayName','V horiz (solo)');
mark(); mark2(); ylabel('V (m/s)'); legend('Location','best');
title('Velocidade — verde = janelas de treino');

subplot(4,1,2); hold on; grid on;
plot(tG, u, 'r','DisplayName','u (frente)');
plot(tG, v, 'g','DisplayName','v (lateral)');
plot(tG, w, 'b','DisplayName','w (baixo)');
mark(); mark2(); ylabel('u,v,w (m/s)'); legend('Location','best');

subplot(4,1,3); hold on; grid on;
plot(tG, alpha,'m','DisplayName','\alpha');
plot(tG, beta, 'k','DisplayName','\beta');
mark(); mark2(); ylabel('\alpha, \beta (°)'); legend('Location','best');

subplot(4,1,4); hold on; grid on;
plot(tG, My_wing_perCm*0.2, 'b','LineWidth',1.3,'DisplayName','M_y asa (Cm=0.2)');
yline(Ly*3,'r--','momento motor ~1 N·m');
mark(); mark2(); ylabel('M_y (N·m)'); xlabel('t (s)'); legend('Location','best');

saveas(fig, fullfile(paths.images,'velocity_diagnostic.png'));
fprintf('\n  Figura: %s\n', fullfile(paths.images,'velocity_diagnostic.png'));
