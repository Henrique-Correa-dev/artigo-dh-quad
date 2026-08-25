% sim_control_scenarios.m — Comparação do controlador LINEAR × NÃO-LINEAR em
% TRÊS níveis de agressividade da manobra (mesma lei, mesma atuação saturada):
%   S1  "estável e igual"  — manobra suave  → linear ≈ não-linear
%   S2  "varia um pouco"   — manobra média  → divergência pequena
%   S3  "varia muito"      — manobra forte  → divergência grande (acoplamento)
%
% Manobra combinada (avanço + guinada), escalada por um fator α que cresce de
% S1→S3: pulso de velocidade de avanço (gera arfagem q) + spin de guinada
% (sustenta r). O produto q·r excita a ROLAGEM via -Γ2·q·r — termo que o modelo
% linear (hover, desacoplado) ignora — e a inclinação reduz o empuxo vertical
% (T·cosθ), fazendo a altitude afundar no NL. Ambos crescem com α.
%
% Cada figura mostra ENTRADAS (T, momentos, PWM), SAÍDAS (h,u,ψ,φ — linear vs
% NL) e SETPOINTS. Uso:  >> design_control;  sim_control_scenarios

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths=setup_paths();
lm=load(fullfile(paths.outputs,'linear_model.mat'));
G =load(fullfile(paths.control,'control_gains.mat'));
proj=parameters(); [fT,fQ]=motor_models(); dyn=vtol_dynamics('get_handles');
M=struct('A',lm.A,'B',lm.B,'u0',lm.u0,'P',lm.P,'K',G.K,'dyn',dyn, ...
         'bridge',forces_to_pwm(lm,'nonlinear'),'fT',fT,'fQ',fQ, ...
         'constants',struct('m',proj.m,'g',proj.g),'dt',0.01,'PWM_LIM',[1200 2000]);

T_end=10; t=(0:M.dt:T_end)'; h0=3;
% três níveis: u_amp = 5α (m/s),  yaw_rate = 50α (°/s)
alpha=[0.25 0.55 0.90];
tag ={'s1_igual','s2_pouco','s3_muito'};
ttl ={'S1 — manobra suave (estável e igual)', ...
      'S2 — manobra média (varia um pouco)', ...
      'S3 — manobra forte (varia muito)'};

fprintf('======================================================================\n');
fprintf(' cenário      | u_amp  yaw   | máx|Δh| máx|Δu| máx|Δψ| máx|Δφ| | estável?\n');
fprintf('-------------+--------------+--------------------------------+---------\n');
for i=1:3
    a=alpha(i);
    s=struct('name',tag{i},'h0',h0,'t',t, ...
        'sp_h',h0*ones(size(t)), ...
        'sp_u',5*a*((t>=3)&(t<7)), ...
        'sp_psi',deg2rad(50*a)*t);
    RL=cl_loop(false,s,M);  RN=cl_loop(true,s,M);

    dh =max(abs(RL.H-RN.H));
    du =max(abs(RL.X(:,7)-RN.X(:,7)));
    dps=rad2deg(max(abs(RL.X(:,6)-RN.X(:,6))));
    dph=rad2deg(max(abs(RL.X(:,4)-RN.X(:,4))));
    stab=(rad2deg(max(abs(RN.X(:,4))))<30) && (dh<2);
    fprintf(' %-11s | %4.1f  %4.0f° | %6.3f %6.3f %6.2f %6.2f | %s\n', ...
        tag{i}, 5*a, 50*a, dh,du,dps,dph, string(stab));

    plot_scenario(tag{i}, ttl{i}, s, RL, RN, M, paths);
end
fprintf('======================================================================\n');


%% ===================== PLOT (2×3: saídas + entradas) =====================
function plot_scenario(tag, ttl, s, RL, RN, M, paths)
    t=s.t;
    amber=[0.93 0.69 0.13]; graymed=[.5 .5 .5]; blue=[0 0.30 0.85];
    fig=figure('Name',['sim_ctrl_' tag],'Position',[55 45 1280 820]);
    set(fig,'Color','w'); try, fig.Theme='light'; catch, end

    % ---- SAÍDAS (com setpoint): h, u, ψ ----  (eixos FIXOS p/ comparar S1..S3)
    subplot(2,3,1); hold on; grid on;
    plot(t,s.sp_h,'--','Color',amber,'LineWidth',1.7,'DisplayName','setpoint');
    plot(t,RL.H,'-','Color',graymed,'LineWidth',1.9,'DisplayName','linear');
    plot(t,RN.H,'-','Color',blue,'LineWidth',1.3,'DisplayName','não-linear');
    ylabel('Altitude h (m)'); ylim([1.3 3.3]); legend('Location','best'); title('Altitude (saída)');

    subplot(2,3,2); hold on; grid on;
    plot(t,s.sp_u,'--','Color',amber,'LineWidth',1.7,'DisplayName','setpoint');
    plot(t,RL.X(:,7),'-','Color',graymed,'LineWidth',1.9,'DisplayName','linear');
    plot(t,RN.X(:,7),'-','Color',blue,'LineWidth',1.3,'DisplayName','não-linear');
    ylabel('Vel. avanço u (m/s)'); ylim([-12 6]); legend('Location','best'); title('Velocidade (saída)');

    subplot(2,3,3); hold on; grid on;
    plot(t,rad2deg(s.sp_psi),'--','Color',amber,'LineWidth',1.7,'DisplayName','setpoint');
    plot(t,rad2deg(RL.X(:,6)),'-','Color',graymed,'LineWidth',1.9,'DisplayName','linear');
    plot(t,rad2deg(RN.X(:,6)),'-','Color',blue,'LineWidth',1.3,'DisplayName','não-linear');
    ylabel('Guinada \psi (°)'); ylim([0 470]); legend('Location','best'); title('Guinada (saída)');

    % ---- ROLAGEM (acoplamento, sp=0) + ENTRADAS: momentos, PWM ----
    subplot(2,3,4); hold on; grid on;
    plot(t,rad2deg(RL.X(:,4)),'-','Color',graymed,'LineWidth',1.9,'DisplayName','linear');
    plot(t,rad2deg(RN.X(:,4)),'-','Color',blue,'LineWidth',1.3,'DisplayName','não-linear');
    yline(0,':','Color',[.7 .7 .7],'HandleVisibility','off');
    ylabel('Rolagem \phi (°)'); ylim([-20 6]); xlabel('t (s)'); legend('Location','best');
    title('Rolagem — setpoint \phi=0 (acoplamento)');

    subplot(2,3,5); hold on; grid on;
    plot(t,RN.U(:,2),'-','Color',[0.85 0.1 0.1],'LineWidth',1.3,'DisplayName','M_x');
    plot(t,RN.U(:,3),'-','Color',[0.1 0.5 0.85],'LineWidth',1.3,'DisplayName','M_y');
    plot(t,RN.U(:,4),'-','Color',[0.2 0.6 0.2],'LineWidth',1.3,'DisplayName','M_z');
    ylabel('Momentos (N·m)'); ylim([-1.4 2]); xlabel('t (s)'); legend('Location','best');
    title('Momentos de controle (entrada, NL)');

    subplot(2,3,6); hold on; grid on;
    co=lines(4);
    for mi=1:4, plot(t,RN.PWM(:,mi),'Color',co(mi,:),'LineWidth',1.0,'DisplayName',sprintf('M%d',mi)); end
    yline(M.PWM_LIM(2),'r--','HandleVisibility','off'); yline(M.PWM_LIM(1),'r--','HandleVisibility','off');
    ylim([M.PWM_LIM(1)-160 M.PWM_LIM(2)+160]);
    ylabel('PWM (\mus)'); xlabel('t (s)'); legend('Location','best','Orientation','horizontal');
    title('PWM motores (entrada, NL — saturado)');

    sgtitle(['Controle LINEAR × NÃO-LINEAR — ' ttl]);
    out=fullfile(paths.images,['sim_ctrl_' tag '.png']);
    exportgraphics(fig,out,'BackgroundColor','white');
    fprintf('   figura: %s\n', out);
end
