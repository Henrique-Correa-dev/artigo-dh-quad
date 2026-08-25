% sim_control_coupling.m — Limiar de divergência LINEAR × NÃO-LINEAR por
% ACOPLAMENTO GIROSCÓPICO.
%
% O modelo linear de hover despreza os termos giroscópicos (p·q, q·r, ...).
% Manobra que os excita: SPIN de guinada contínuo (sustenta r) + um pulso de
% velocidade de avanço (gera q). O produto q·r gera, via -Γ2·q·r, um torque de
% ROLAGEM que o modelo linear não enxerga (prevê φ≡0). Mostra-se:
%   (a) um caso representativo estável (spin 45°/s): φ_linear≈0 vs φ_NL≈6°;
%   (b) a curva de LIMIAR: o desvio cresce com o spin e, além de ~55°/s, o
%       controlador (projetado no linear) PERDE ESTABILIDADE na planta NL —
%       limite que o próprio modelo linear não prevê.
%
% Uso:  >> design_control;  sim_control_coupling

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths=setup_paths();
lm=load(fullfile(paths.outputs,'linear_model.mat'));
G =load(fullfile(paths.control,'control_gains.mat'));
proj=parameters(); [fT,fQ]=motor_models(); dyn=vtol_dynamics('get_handles');
M=struct('A',lm.A,'B',lm.B,'u0',lm.u0,'P',lm.P,'K',G.K,'dyn',dyn, ...
         'bridge',forces_to_pwm(lm,'nonlinear'),'fT',fT,'fQ',fQ, ...
         'constants',struct('m',proj.m,'g',proj.g),'dt',0.01,'PWM_LIM',[1200 2000]);

% Γ2 do modelo (corpo rígido) — o ganho do torque giroscópico de rolagem
Jx=M.P(1); Jy=M.P(2); Jz=M.P(3); Jxz=M.P(4);  g0=Jx*Jz-Jxz^2;
G2=(Jz*(Jz-Jy)+Jxz^2)/g0;

%% ===================== (a) CASO REPRESENTATIVO (spin 45°/s) =====================
yr0=45;  T_end=10;  t=(0:M.dt:T_end)';
sc=struct('name','coup','h0',3,'t',t, ...
    'sp_h',3*ones(size(t)),'sp_u',3*((t>=3)&(t<6.5)),'sp_psi',deg2rad(yr0)*t);
RL=cl_loop(false,sc,M);  RN=cl_loop(true,sc,M);

phiL=rad2deg(RL.X(:,4)); phiN=rad2deg(RN.X(:,4));
qN=rad2deg(RN.X(:,2));   rN=rad2deg(RN.X(:,3));
dist=-G2*(RN.X(:,2).*RN.X(:,3));            % torque giroscópico de rolagem (rad/s²)

%% ===================== (b) VARREDURA DE LIMIAR =====================
rates=[0 10 20 30 40 45 50 55 60 70];
dphi=zeros(size(rates)); stab=false(size(rates));
for i=1:numel(rates)
    s=struct('name','sw','h0',3,'t',t, ...
        'sp_h',3*ones(size(t)),'sp_u',3*((t>=3)&(t<6.5)),'sp_psi',deg2rad(rates(i))*t);
    rl=cl_loop(false,s,M); rn=cl_loop(true,s,M);
    dphi(i)=rad2deg(max(abs(rl.X(:,4)-rn.X(:,4))));
    stab(i)=(rad2deg(max(abs(rn.X(:,4))))<25) && (rad2deg(max(abs(rn.X(:,5))))<40) ...
            && (max(abs(rl.H-rn.H))<1.5);     % "controlado": rolagem regulada <25°
end
i_unst=find(~stab,1);                        % 1ª taxa instável
if isempty(i_unst), thr=NaN; else, thr=mean(rates([i_unst-1 i_unst])); end

fprintf('==========================================================\n');
fprintf('  ACOPLAMENTO GIROSCÓPICO — limiar linear×NL\n');
fprintf('  Γ2 = %.3f | caso (a) spin %d°/s: φ_lin máx=%.2f°, φ_NL máx=%.2f°\n', ...
    G2, yr0, max(abs(phiL)), max(abs(phiN)));
fprintf('  Limiar de estabilidade (controlador linear na planta NL): ~%.0f°/s\n', thr);
fprintf('  rate(°/s):'); fprintf(' %4d',rates); fprintf('\n');
fprintf('  Δφ  (°)  :'); fprintf(' %4.1f',min(dphi,99)); fprintf('\n');
fprintf('  estável  :'); fprintf('  %s',string(stab)); fprintf('\n');

%% ===================== PLOT =====================
fig=figure('Name','sim_control_coupling','Position',[60 50 1250 850]);
set(fig,'Color','w'); try, fig.Theme='light'; catch, end
amber=[0.93 0.69 0.13]; graymed=[.5 .5 .5]; blue=[0 0.30 0.85]; red=[0.85 0.1 0.1];

% (1) Rolagem: linear (≡0) vs NL
subplot(2,2,1); hold on; grid on;
plot(t,phiL,'-','Color',graymed,'LineWidth',1.9,'DisplayName','linear');
plot(t,phiN,'-','Color',blue,'LineWidth',1.4,'DisplayName','não-linear');
yline(0,':','Color',[.7 .7 .7],'HandleVisibility','off');
ylabel('Rolagem \phi (°)'); legend('Location','best');
title(sprintf('Rolagem induzida (spin %d°/s) — linear prevê \\phi\\equiv0',yr0));

% (2) As taxas simultâneas q e r (NL)
subplot(2,2,2); hold on; grid on;
plot(t,qN,'-','Color',[0.1 0.6 0.2],'LineWidth',1.4,'DisplayName','q (arfagem)');
plot(t,rN,'-','Color',[0.6 0.2 0.7],'LineWidth',1.4,'DisplayName','r (guinada)');
ylabel('Taxas (°/s)'); legend('Location','best'); title('Taxas simultâneas (não-linear)');

% (3) O torque giroscópico de rolagem -Γ2·q·r (o que o linear ignora)
subplot(2,2,3); hold on; grid on;
plot(t,dist,'-','Color',red,'LineWidth',1.5,'DisplayName','-\Gamma_2\,q\,r');
yline(0,':','Color',[.7 .7 .7],'HandleVisibility','off');
ylabel('Torque giro. de rolagem (rad/s^2)'); xlabel('t (s)');
legend('Location','best'); title('Termo -\Gamma_2 q r (ausente no modelo linear)');

% (4) Curva de limiar: desvio de rolagem vs spin + região de tombamento
subplot(2,2,4); hold on; grid on;
ymax=1.25*max(dphi(stab));
if ~isnan(thr)
    xb=[thr max(rates) max(rates) thr];
    patch(xb,[0 0 ymax ymax],red,'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
    xline(thr,'--','Color',red,'LineWidth',1.4,'DisplayName',sprintf('limiar ~%.0f°/s',thr));
end
plot(rates(stab),dphi(stab),'-o','Color',blue,'LineWidth',1.6,'MarkerFaceColor',blue,'DisplayName','estável');
plot(rates(~stab),min(dphi(~stab),ymax),'x','Color',red,'MarkerSize',11,'LineWidth',2,'DisplayName','tombamento (NL)');
ylim([0 ymax]); xlim([0 max(rates)]);
ylabel('máx|\phi_{lin}-\phi_{NL}| (°)'); xlabel('taxa de spin (°/s)');
legend('Location','northwest'); title('Limiar: desvio de rolagem vs spin');

sgtitle('Controle LINEAR × NÃO-LINEAR — acoplamento giroscópico (taxas altas)');
out=fullfile(paths.images,'sim_compare_coupling.png');
exportgraphics(fig,out,'BackgroundColor','white');
fprintf('  Figura: %s\n',out);
fprintf('==========================================================\n');
