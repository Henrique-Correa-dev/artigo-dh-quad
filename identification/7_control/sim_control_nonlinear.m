% sim_control_nonlinear.m — Malha fechada do MESMO controlador na PLANTA NÃO-LINEAR
%
% Caminho realista (Diagrama B):
%   control_law → forces_to_pwm (ponte) → PWM → vtol_dynamics (NL) → estados
%
% A cada passo: controlador comanda u=[T,Mx,My,Mz]; a PONTE aloca em PWM;
% a planta NÃO-LINEAR reconverte PWM→forças e integra os 9 estados. A dupla
% conversão expõe a imperfeição da alocação (que a realimentação corrige).
%
% Altitude e trajetória: integra a posição INERCIAL (NED) via rotação
% R_body→NED de [u v w]. ḣ = -ż_ned.  (no hover nivelado → ḣ = -w, igual linear)
%
% Compara com o resultado do MODELO LINEAR (sim_control_linear_result.mat).
%
% Uso:  >> design_control; sim_control_linear; sim_control_nonlinear

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ===================== CARREGAR =====================
lm = load(fullfile(paths.outputs,'linear_model.mat'));
G  = load(fullfile(paths.control,'control_gains.mat'));  K = G.K;
proj = parameters();  constants = struct('m',proj.m,'g',proj.g);
[fT,fQ] = motor_models();
dyn = vtol_dynamics('get_handles');
P   = lm.P;

BRIDGE_MODE = 'nonlinear';                 % 'linear' | 'nonlinear'
bridge = forces_to_pwm(lm, BRIDGE_MODE);

%% ===================== CENÁRIO (igual ao linear) =====================
dt = 0.01;  T_end = 25;
t  = (0:dt:T_end)';  N = numel(t);
sp_h   = 2.0 * (t >= 1);
sp_u   = 3.0 * (t >= 8);
sp_psi = deg2rad(30) * (t >= 16);

%% ===================== INTEGRAÇÃO (RK4, ZOH no controlador) =====================
x = zeros(9,1);                 % [p q r φ θ ψ u v w]
pos = zeros(3,1);               % posição inercial NED [N;E;D]
cs = struct('int_h',0,'int_u',0);

X = zeros(N,9);  H = zeros(N,1);  U = zeros(N,4);  PWM = zeros(N,4);  POS = zeros(N,3);
for k = 1:N
    h = -pos(3);                            % altitude = -D
    sp = struct('h',sp_h(k),'u',sp_u(k),'psi',sp_psi(k));
    [u_cmd, cs] = control_law(x, h, sp, K, cs, dt);     % controlador (ZOH)
    pwm = bridge(u_cmd);                                % PONTE → PWM
    [Tt,Mx,My,Mz] = dyn.forces(pwm', P, fT, fQ);        % PLANTA: PWM → forças reais

    X(k,:)=x';  H(k)=h;  U(k,:)=[Tt Mx My Mz];  PWM(k,:)=pwm';  POS(k,:)=pos';

    % derivadas (PWM/forças fixos no passo = ZOH)
    fx = @(xx) dyn.rigid_body(xx, Tt, Mx, My, Mz, P, constants);
    fp = @(xx) body2ned_vel(xx);            % ṗos_ned = R·[u;v;w]
    k1=fx(x);            p1=fp(x);
    k2=fx(x+dt/2*k1);    p2=fp(x+dt/2*k1);
    k3=fx(x+dt/2*k2);    p3=fp(x+dt/2*k2);
    k4=fx(x+dt*k3);      p4=fp(x+dt*k3);
    x   = x   + dt/6*(k1+2*k2+2*k3+k4);
    pos = pos + dt/6*(p1+2*p2+2*p3+p4);
end

%% ===================== MÉTRICAS =====================
fprintf('==========================================================\n');
fprintf('  MALHA FECHADA — PLANTA NÃO-LINEAR (ponte: %s)\n', BRIDGE_MODE);
fprintf('==========================================================\n');
report_step('Altitude h', t, H,             sp_h,   1,  'm');
report_step('Vel fwd  u', t, X(:,7),        sp_u,   8,  'm/s');
report_step('Heading psi',t, rad2deg(X(:,6)), rad2deg(sp_psi), 16, '°');
fprintf('  Pitch máx: %.1f° | Roll máx: %.2f° | PWM faixa: [%.0f, %.0f] µs\n', ...
    rad2deg(max(abs(X(:,5)))), rad2deg(max(abs(X(:,4)))), min(PWM(:)), max(PWM(:)));

% Comparação com o linear
L = [];
f_lin = fullfile(paths.control,'sim_control_linear_result.mat');
if exist(f_lin,'file'), tmp = load(f_lin); L = tmp.res_lin; end

%% ===================== PLOT =====================
fig = figure('Name','sim_control_nonlinear','Position',[60 40 1300 850]);
sb = @(i) subplot(3,2,i);
spc = [1 0.85 0.1];   % cor do setpoint (âmbar — bem visível no fundo escuro)

sb(1); hold on; grid on;
plot(t, sp_h,'--','Color',spc,'LineWidth',1.8,'DisplayName','setpoint');
if ~isempty(L), plot(L.t,L.H,'Color',[.6 .6 .6],'LineWidth',1.6,'DisplayName','linear'); end
plot(t, H,'b-','LineWidth',1.4,'DisplayName','não-linear');
ylabel('Altitude (m)'); legend('Location','best'); title('Altitude');

sb(2); hold on; grid on;
plot(t, sp_u,'--','Color',spc,'LineWidth',1.8,'DisplayName','setpoint');
if ~isempty(L), plot(L.t,L.X(:,7),'Color',[.6 .6 .6],'LineWidth',1.6,'DisplayName','linear'); end
plot(t, X(:,7),'b-','LineWidth',1.4,'DisplayName','não-linear');
ylabel('Vel forward (m/s)'); legend('Location','best'); title('Velocidade forward');

sb(3); hold on; grid on;
plot(t, rad2deg(sp_psi),'--','Color',spc,'LineWidth',1.8,'DisplayName','setpoint');
if ~isempty(L), plot(L.t,rad2deg(L.X(:,6)),'Color',[.6 .6 .6],'LineWidth',1.6,'DisplayName','linear'); end
plot(t, rad2deg(X(:,6)),'b-','LineWidth',1.4,'DisplayName','não-linear');
ylabel('Heading (°)'); legend('Location','best'); title('Heading');

sb(4); hold on; grid on;
if ~isempty(L), plot(L.t,rad2deg(L.X(:,5)),'Color',[.6 .6 .6],'LineWidth',1.6,'DisplayName','\theta linear'); end
plot(t, rad2deg(X(:,5)),'b-','LineWidth',1.4,'DisplayName','\theta NL');
plot(t, rad2deg(X(:,4)),'g-','LineWidth',1.0,'DisplayName','\phi NL');
ylabel('Atitude (°)'); legend('Location','best'); title('Pitch/Roll');

sb(5); hold on; grid on;
plot(t, PWM(:,1),'DisplayName','M1'); plot(t, PWM(:,2),'DisplayName','M2');
plot(t, PWM(:,3),'DisplayName','M3'); plot(t, PWM(:,4),'DisplayName','M4');
yline(1000,'k:'); yline(2000,'k:');
ylabel('PWM (µs)'); xlabel('t (s)'); legend('Location','best'); title('PWM por motor (saída da ponte)');

sb(6); hold on; grid on; axis equal;
plot3(POS(:,2), POS(:,1), -POS(:,3),'b-','LineWidth',1.4);   % E, N, h
plot3(POS(1,2),POS(1,1),-POS(1,3),'go','MarkerFaceColor','g');
plot3(POS(end,2),POS(end,1),-POS(end,3),'rs','MarkerFaceColor','r');
xlabel('Leste (m)'); ylabel('Norte (m)'); zlabel('Alt (m)');
title('Trajetória 3D'); view(45,25); grid on;

sgtitle(sprintf('Controle cascata — PLANTA NÃO-LINEAR vs LINEAR (ponte %s)', BRIDGE_MODE));
out = fullfile(paths.images,'sim_control_nonlinear.png');
saveas(fig, out);
fprintf('\n  Figura: %s\n', out);
fprintf('==========================================================\n');


%% ===================== HELPERS =====================
function vned = body2ned_vel(x)
% Velocidade inercial NED = R_body→NED · [u;v;w]  (Euler ZYX φ,θ,ψ)
    phi=x(4); th=x(5); ps=x(6); u=x(7); v=x(8); w=x(9);
    cphi=cos(phi); sphi=sin(phi); cth=cos(th); sth=sin(th); cps=cos(ps); sps=sin(ps);
    R = [ cth*cps,  sphi*sth*cps-cphi*sps,  cphi*sth*cps+sphi*sps; ...
          cth*sps,  sphi*sth*sps+cphi*cps,  cphi*sth*sps-sphi*cps; ...
         -sth,      sphi*cth,               cphi*cth ];
    vned = R*[u;v;w];
end

function report_step(name, t, y, sp, t_step, unit)
    final = sp(end);
    y_ss = mean(y(t >= (t(end)-2)));
    after = t >= t_step;
    if final ~= 0, os = (max(abs(y(after)))-abs(final))/abs(final)*100; else, os = NaN; end
    fprintf('  %-12s alvo=%6.2f %-4s  final=%6.2f  erro=%+6.3f  overshoot=%5.1f%%\n', ...
        name, final, unit, y_ss, y_ss-final, os);
end
