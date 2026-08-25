% sim_control_linear.m — Malha fechada do controle cascata no MODELO LINEAR
%
% Integra o modelo linear AUMENTADO (9 estados + altitude h) com a lei
% control_law.m. Aplica degraus de setpoint em altitude, vel. forward e
% heading — valida que a ESTRUTURA de controle fecha e segue referência.
%
% Modelo linear: u = [T,Mx,My,Mz] entra DIRETO (sem alocação) — Diagrama A.
%   ẋ = A·x + B·(u_cmd - u0),   ḣ = -w
%
% Uso:  >> design_control   (gera ganhos)   depois   >> sim_control_linear

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ===================== CARREGAR MODELO + GANHOS =====================
lm = load(fullfile(paths.outputs, 'linear_model.mat'));
A = lm.A;  B = lm.B;  u0 = lm.u0;          % u0 = [m·g;0;0;0]
G  = load(fullfile(paths.control, 'control_gains.mat'));  K = G.K;

%% ===================== CENÁRIO DE SETPOINTS =====================
dt = 0.01;  T_end = 25;
t  = (0:dt:T_end)';  N = numel(t);

% Degraus: altitude (1s) → vel forward (8s) → heading (16s)
sp_h   = 2.0  * (t >= 1);     % subir 2 m
sp_u   = 3.0  * (t >= 8);     % 3 m/s pra frente
sp_psi = deg2rad(30) * (t >= 16);   % girar 30°

%% ===================== ALOCAÇÃO DE MOTORES (p/ saturação física) =====================
%  A planta recebe o [T,τx,τy,τz] EFETIVO que os motores conseguem entregar:
%  comando -> empuxo/motor (B⁻¹) -> PWM -> satura [1000,2000] µs -> empuxo
%  efetivo -> [T,τ] efetivo (B). Reproduz o limite real dos atuadores na malha.
pp = parameters();  Pf = lm.P(:);                       % P do linear_model (OEM final)
Lxr=pp.arms.Lx_r; Lxl=pp.arms.Lx_l; Lyf=pp.arms.Ly_f; Lyr=pp.arms.Ly_r;
cMcT = mean(Pf(9:12))*pp.motor.kQ_rpm / (mean(Pf(5:8))*pp.motor.kT_rpm);
Balloc = [1 1 1 1; -Lxr Lxl Lxl -Lxr; Lyf -Lyr Lyf -Lyr; cMcT cMcT -cMcT -cMcT];
kTm  = Pf(5:8) .* pp.motor.kT_rpm;                      % N/RPM² por motor (coluna 4×1)
CR=pp.motor.CR;  Om_b=pp.motor.Omega_b;  PWM_LIM=[1200 2000];   % idle 1200 µs (rotor mín. armado)
Tmax1 = max(kTm) * (CR*PWM_LIM(2)+Om_b)^2;              % empuxo máx por motor [N]

%% ===================== INTEGRAÇÃO (RK4) =====================
x = zeros(9,1);  h = 0;                 % parte do hover nivelado
cs = struct('int_h',0,'int_u',0);

X = zeros(N,9);  H = zeros(N,1);  U = zeros(N,4);  TH_c = zeros(N,1);
PWM = zeros(N,4);  PWMdem = zeros(N,4);  Usat = false(N,1);
for k = 1:N
    sp = struct('h',sp_h(k),'u',sp_u(k),'psi',sp_psi(k));
    [u_cmd, cs, dbg] = control_law(x, h, sp, K, cs, dt);

    % --- alocação -> PWM -> SATURAÇÃO -> [T,τ] efetivo (o que a planta recebe) ---
    Ti    = Balloc \ u_cmd;                             % 4×1 empuxo demandado/motor
    pwm   = (sqrt(max(Ti,0)./kTm) - Om_b) / CR;         % empuxo -> PWM demandado
    pwmS  = min(max(pwm, PWM_LIM(1)), PWM_LIM(2));      % satura [1000,2000] µs
    TiS   = kTm .* (CR*pwmS + Om_b).^2;                 % PWM saturado -> empuxo efetivo
    u_eff = Balloc * TiS;                               % [T,τ] efetivo na planta

    X(k,:) = x';  H(k) = h;  U(k,:) = u_eff';  TH_c(k) = dbg.theta_c;
    PWM(k,:) = pwmS';  PWMdem(k,:) = pwm';
    Usat(k)  = any(pwm < PWM_LIM(1) | pwm > PWM_LIM(2));

    % derivadas (modelo linear aumentado) — usa o comando EFETIVO (saturado)
    f  = @(xx,hh) deal(A*xx + B*(u_eff - u0), -xx(9));   % ẋ , ḣ=-w
    [k1x,k1h] = f(x,         h);
    [k2x,k2h] = f(x+dt/2*k1x, h+dt/2*k1h);
    [k3x,k3h] = f(x+dt/2*k2x, h+dt/2*k2h);
    [k4x,k4h] = f(x+dt*k3x,   h+dt*k3h);
    x = x + dt/6*(k1x+2*k2x+2*k3x+k4x);
    h = h + dt/6*(k1h+2*k2h+2*k3h+k4h);
end

%% ===================== MÉTRICAS =====================
fprintf('==========================================================\n');
fprintf('  MALHA FECHADA — MODELO LINEAR\n');
fprintf('==========================================================\n');
report_step('Altitude h', t, H,           sp_h,   1,  'm');
report_step('Vel fwd  u', t, X(:,7),      sp_u,   8,  'm/s');
report_step('Heading psi',t, rad2deg(X(:,6)), rad2deg(sp_psi), 16, '°');
fprintf('  Pitch máx: %.1f° | Roll máx: %.2f° | T faixa: [%.1f, %.1f] N\n', ...
    rad2deg(max(abs(X(:,5)))), rad2deg(max(abs(X(:,4)))), min(U(:,1)), max(U(:,1)));

%% ===================== SATURAÇÃO DOS MOTORES (resumo) =====================
nsat = sum(Usat);
fprintf('  PWM demandado: [%.0f, %.0f] µs | aplicado (saturado): [%.0f, %.0f] µs\n', ...
    min(PWMdem(:)), max(PWMdem(:)), min(PWM(:)), max(PWM(:)));
fprintf('  Empuxo/motor lim=%.2f N | %d/%d instantes com saturação (%.1f%%)\n', ...
    Tmax1, nsat, N, 100*nsat/N);

%% ===================== PLOT =====================
fig = figure('Name','sim_control_linear','Position',[60 15 1250 1050]);
set(fig,'Color','w'); try, fig.Theme = 'light'; catch, end
sb = @(i) subplot(4,2,i);
spc = [1 0.85 0.1];   % cor do setpoint (âmbar — bem visível no fundo escuro)

sb(1); hold on; grid on;
plot(t, sp_h,'--','Color',spc,'LineWidth',1.8,'DisplayName','setpoint');
plot(t, H,'b-','LineWidth',1.5,'DisplayName','h');
ylabel('Altitude (m)'); legend('Location','best'); title('Altitude');

sb(2); hold on; grid on;
plot(t, sp_u,'--','Color',spc,'LineWidth',1.8,'DisplayName','setpoint');
plot(t, X(:,7),'b-','LineWidth',1.5,'DisplayName','u');
ylabel('Vel forward (m/s)'); legend('Location','best'); title('Velocidade forward');

sb(3); hold on; grid on;
plot(t, rad2deg(sp_psi),'--','Color',spc,'LineWidth',1.8,'DisplayName','setpoint');
plot(t, rad2deg(X(:,6)),'b-','LineWidth',1.5,'DisplayName','\psi');
ylabel('Heading (°)'); legend('Location','best'); title('Heading');

sb(4); hold on; grid on;
plot(t, rad2deg(TH_c),'r--','LineWidth',1.0,'DisplayName','\theta_c (comando)');
plot(t, rad2deg(X(:,5)),'b-','LineWidth',1.5,'DisplayName','\theta');
plot(t, rad2deg(X(:,4)),'g-','LineWidth',1.0,'DisplayName','\phi');
ylabel('Atitude (°)'); legend('Location','best'); title('Pitch/Roll (interno)');

sb(5); hold on; grid on;
plot(t, U(:,1),'b-','LineWidth',1.3,'DisplayName','T (N)');
yline(u0(1),':','peso m·g','Color',[0.8 0.8 0.8],'HandleVisibility','off');
ylabel('Empuxo T (N)'); xlabel('t (s)');
legend('Location','best'); title('Empuxo total');

sb(6); hold on; grid on;
plot(t, U(:,2),'r-','DisplayName','Mx'); plot(t, U(:,3),'b-','DisplayName','My');
plot(t, U(:,4),'g-','DisplayName','Mz');
ylabel('Momentos (N·m)'); legend('Location','best'); title('Momentos');

subplot(4,2,[7 8]); hold on; grid on;
co = lines(4);
for mi = 1:4
    plot(t, PWMdem(:,mi),':','Color',[co(mi,:) 0.45],'LineWidth',0.9,'HandleVisibility','off');
    plot(t, PWM(:,mi),'-','Color',co(mi,:),'LineWidth',1.2,'DisplayName',sprintf('M%d',mi));
end
yline(2000,'r--','HandleVisibility','off'); yline(1200,'r--','HandleVisibility','off');
ylim([850 2150]);
ylabel('PWM (\mus)'); xlabel('t (s)'); legend('Location','best','Orientation','horizontal');
title('PWM dos motores — saturado em [1200, 2000] \mus (tracejado = demandado)');

sgtitle('Controle cascata — malha fechada no MODELO LINEAR (com saturação dos atuadores)');
out = fullfile(paths.images,'sim_control_linear.png');
exportgraphics(fig, out, 'BackgroundColor','white');

% Salvar resultados p/ sobreposição no sim NÃO-LINEAR
res_lin = struct('t',t,'H',H,'X',X,'U',U,'sp_h',sp_h,'sp_u',sp_u,'sp_psi',sp_psi);
save(fullfile(paths.control,'sim_control_linear_result.mat'),'res_lin');
fprintf('\n  Figura: %s\n', out);
fprintf('==========================================================\n');


%% ===================== HELPER =====================
function report_step(name, t, y, sp, t_step, unit)
    final = sp(end);
    idx_ss = t >= (t(end) - 2);          % média dos últimos 2 s
    y_ss = mean(y(idx_ss));
    err  = y_ss - final;
    % overshoot relativo ao degrau
    after = t >= t_step;
    if final ~= 0
        os = (max(abs(y(after))) - abs(final)) / abs(final) * 100;
    else
        os = NaN;
    end
    fprintf('  %-12s alvo=%6.2f %-4s  final=%6.2f  erro=%+6.3f  overshoot=%5.1f%%\n', ...
        name, final, unit, y_ss, err, os);
end
