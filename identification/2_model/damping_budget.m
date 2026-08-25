% damping_budget.m — de onde vem o amortecimento c_p, c_q, c_r do multirrotor
% =========================================================================
% Monta o BALANÇO COMPLETO do amortecimento angular no pairado, mecanismo por
% mecanismo, cada um calculado da geometria e de dado medido, sem simplificar
% nenhum deles em "um coeficiente efetivo". No fim compara a soma com os
% valores identificados por EEM (erro de equação, sem integração) e por OEM
% (erro de saída, com integração).
%
% MECANISMOS (todos em 1/s, já divididos pela inércia do eixo)
%
%  [A] INFLUXO — empuxo diferencial pela rotação do corpo.
%      Girando em rolagem a p, o rotor em +y sobe a p·ly e o em −y desce; o
%      empuxo de cada um responde com −k_v por m/s de velocidade axial:
%          c_p = 4·k_v·ly²/J_x        (Padfield 2007, Leishman 2006)
%      k_v = (T0/v_h)·κ com κ = aσ/(16λ_i + aσ). O κ é a correção QUASE-ESTÁTICA
%      do elemento de pá: sem ela (influxo congelado) k_v seria T0/v_h·1, que é
%      1/κ ≈ 4 vezes maior. Qual dos dois vale depende da frequência: o tempo de
%      relaxação do influxo é τ = (4/3π)·R/v_i, calculado abaixo.
%
%  [B] MOMENTO DE CUBO — a pá vê variação cíclica de incidência pela rotação do
%      corpo. Integrando a sustentação da seção ao longo da pá e em azimute,
%          M/q por rotor = N_b·ρ·a·c_pá·Ω·R⁴/16
%      É o análogo do amortecimento de rotor rígido de helicóptero, sem batimento.
%
%  [C] FORÇA H — o plano dos rotores fica a uma altura h acima do CG, então uma
%      taxa de rolagem dá velocidade NO PLANO do rotor de p·h, e a força H se
%      opõe:  c = 4·k_h·h²/J.
%
%  [D] ESTEIRA (download) — parte da estrutura fica na esteira dos rotores e
%      sofre uma força vertical de placa plana. Girando o corpo, cada elemento
%      sobe ou desce em relação à esteira, e como a força é ½ρC_D S w², a
%      derivada em torno de w já é LINEAR: ∂F/∂w = 2·D/w. O download D é medido
%      pelo excesso de empuxo no pairado (T/mg − 1). É o único termo que depende
%      de v_i, não da velocidade de voo, e é o que justifica um coeficiente
%      CONSTANTE: v_i é fixado pelo peso, não pela velocidade de translação.
%
%  [E] AERODINÂMICA DA ESTRUTURA na velocidade de voo, com C_lp, C_mq e C_nr
%      MEDIDOS nos doublets de asa fixa (11_fixed_wing), não os da AVL.
%
% Uso:  >> damping_budget
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
rho = p.rho;  g = p.g;  m = p.m;

%% ---------------- rotor (hélice 1045, 2 pás) — mesmas hipóteses do prior_damping
R = 0.127;  A = pi*R^2;  Nb = 2;  c_blade = 0.025;  a_lift = 5.0;  Cd0 = 0.03;
pitch_in = 4.5;
sigma = Nb*c_blade/(pi*R);
theta75 = atan(pitch_in/(0.75*2*R/0.0254*pi));
T0 = m*g/4;
rpm_h = interp1(p.bench.T_grams*1e-3*g, p.bench.RPM, T0, 'linear');
Om = rpm_h*2*pi/60;  OmR = Om*R;
vh = sqrt(T0/(2*rho*A));  lam = vh/OmR;
kap = a_lift*sigma/(16*lam + a_lift*sigma);
k_v = (T0/vh)*kap;                                   % N/(m/s) quase-estático
k_v_frozen = T0/vh;                                  % limite de influxo congelado
k_h = rho*A*OmR*(sigma/4)*(Cd0 + a_lift*theta75*lam);
tau_inflow = (4/(3*pi))*R/vh;                        % relaxação do influxo [s]

ly = mean([p.arms.Lx_r, p.arms.Lx_l]);      % braço lateral  (rolagem)
lx = mean([p.arms.Ly_f, p.arms.Ly_r]);      % braço longitudinal (arfagem)
d2 = lx^2 + ly^2;
h_rot = abs(p.wing.Z_CG);                   % plano dos rotores acima do CG [m]
J = [p.J.Jx, p.J.Jy, p.J.Jz];

%% ---------------- download medido pelo excesso de empuxo no pairado
DL_FRAC = 0.078;              % T/mg − 1 na janela de validação (k_T identificado)
D_dl = DL_FRAC*m*g;           % [N]
w_wake = 1.5*vh;              % velocidade da esteira onde a estrutura está (v_i a 2v_i)
kD = 2*D_dl/w_wake;           % [N/(m/s)] derivada da força de placa plana

%% ---------------- aerodinâmica da estrutura, coeficientes MEDIDOS em asa fixa
FWl = load(fullfile(paths.outputs,'fw_longitudinal.mat'));
FWt = load(fullfile(paths.outputs,'fw_lateral.mat'));
gi = @(S_,nm) mean(S_.Pid(strcmp(S_.pn,nm),:));
Clp = gi(FWt,'Clp');  Cmq = gi(FWl,'Cmq');  Cnr = gi(FWt,'Cnr');
S = p.wing.S;  b = p.wing.b;  cbar = p.wing.c;
V_ENV = 1.29;                 % p95 da velocidade estimada no trecho de treino

%% ---------------- balanço
B = zeros(5,3);   % [mecanismo x eixo]
B(1,:) = [4*k_v*ly^2, 4*k_v*lx^2, 0] ./ J;                       % A influxo
B(2,:) = [1 1 0]*Nb*rho*a_lift*c_blade*Om*R^4/16*4 ./ J(1:3);    % B momento de cubo
B(2,3) = 0;
B(3,:) = [4*k_h*h_rot^2, 4*k_h*h_rot^2, 4*k_h*d2] ./ J;          % C força H
B(4,:) = [kD*ly^2, kD*lx^2, 0] ./ J;                             % D esteira
B(5,:) = [0.25*rho*S*b^2*abs(Clp), 0.25*rho*S*cbar^2*abs(Cmq), ...
          0.25*rho*S*b^2*abs(Cnr)] * V_ENV ./ J;                 % E aerodinâmica

mech = {'A influxo (empuxo diferencial)','B momento de cubo do rotor', ...
        'C força H (rotor acima do CG)','D esteira sobre a estrutura', ...
        'E aerodinâmica da estrutura'};
eixo = {'c_p (rolagem)','c_q (arfagem)','c_r (guinada)'};

EEM = [2.412, 2.188, 0.646];      % outputs/runs/aero_off_2026 (erro de equação)
OEM = [5.500, 3.861, 0.777];      % outputs/P_identified.mat  (erro de saída, oficial)

fprintf('\n  ======================================================================\n');
fprintf('   ROTOR: %.0f rpm no pairado | v_i = %.2f m/s | λ_i = %.4f | κ = %.3f\n', rpm_h, vh, lam, kap);
fprintf('   relaxação do influxo τ = %.4f s (%.0f Hz) → nos doublets (1 a 3 Hz)\n', tau_inflow, 1/(2*pi*tau_inflow));
fprintf('   o influxo é QUASE-ESTÁTICO, então vale k_v = %.3f e não %.3f N/(m/s)\n', k_v, k_v_frozen);
fprintf('   download medido pelo excesso de empuxo: %.2f N (%.1f%% do peso)\n', D_dl, 100*DL_FRAC);
fprintf('   C_lp, C_mq, C_nr do voo de asa fixa: %.3f, %.2f, %.3f\n', Clp, Cmq, Cnr);
fprintf('  ======================================================================\n\n');
fprintf('  %-32s %10s %10s %10s\n', 'mecanismo', eixo{:});
for i = 1:5
    fprintf('  %-32s %10.3f %10.3f %10.3f\n', mech{i}, B(i,:));
end
fprintf('  %-32s %10.3f %10.3f %10.3f\n', 'SOMA (física, sem ajuste)', sum(B,1));
fprintf('  %-32s %10.3f %10.3f %10.3f\n', 'identificado por EEM', EEM);
fprintf('  %-32s %10.3f %10.3f %10.3f\n', 'identificado por OEM (oficial)', OEM);
fprintf('  %-32s %9.0f%% %9.0f%% %9.0f%%\n', 'soma / EEM', 100*sum(B,1)./EEM);
fprintf('  %-32s %9.0f%% %9.0f%% %9.0f%%\n', 'soma / OEM', 100*sum(B,1)./OEM);
fprintf('\n  Reynolds da asa: %.0f no pairado (V = %.1f m/s) contra %.0f em asa fixa (17 m/s)\n', ...
    rho*V_ENV*cbar/1.81e-5, V_ENV, rho*17*cbar/1.81e-5);

save(fullfile(paths.outputs,'damping_budget.mat'), 'B','mech','eixo','EEM','OEM', ...
     'k_v','k_h','vh','kap','tau_inflow','D_dl','Clp','Cmq','Cnr','V_ENV');

%% ---------------- figura
f = figure('Position',[40 40 1150 520],'Color','w'); try, f.Theme='light'; catch, end
cor = [0.20 0.40 0.65; 0.35 0.60 0.80; 0.55 0.75 0.90; 0.85 0.55 0.20; 0.45 0.70 0.35];
for i = 1:3
    subplot(1,3,i); hold on; grid on;
    hb = bar([1 2 3], [B(:,i)'; zeros(1,5); zeros(1,5)], 'stacked');
    for k = 1:5, hb(k).FaceColor = cor(k,:); end
    bar(2, EEM(i), 0.6, 'FaceColor',[0.75 0.75 0.75], 'EdgeColor','k');
    bar(3, OEM(i), 0.6, 'FaceColor',[0.45 0.45 0.45], 'EdgeColor','k');
    set(gca,'XTick',1:3,'XTickLabel',{'física','EEM','OEM'});
    ylabel('amortecimento [1/s]'); title(eixo{i});
    ylim([0 max(OEM(i),sum(B(:,i)))*1.25]);
    if i == 1
        legend(hb, mech, 'Location','northwest', 'FontSize',7.5);
    end
end
sgtitle('Balanço físico do amortecimento no pairado, contra o identificado por erro de equação e por erro de saída');
fn = fullfile(paths.images,'damping_budget.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',140);
fprintf('  Figura: %s\n', fn);
