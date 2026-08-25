% fw_vs_quad_damping.m — o amortecimento do multirrotor pode ser aerodinâmico?
% =========================================================================
% Junta as duas campanhas para responder, com dado de voo e não com estimativa
% de painéis, quanto do amortecimento identificado no multirrotor (c_p, c_q,
% c_r) pode vir da aerodinâmica da estrutura.
%
% Lógica: a aerodinâmica da estrutura escala com a velocidade,
%     L_p(V) = ¼ρSb²·C_lp·V / J_x    [1/s]      (idem M_q e N_r)
% então basta medir C_lp, C_mq e C_nr onde eles são grandes e bem observáveis,
% que é em voo à frente a 16 a 18 m/s (voos de asa fixa de dez/2024, em MANUAL,
% com doublets de superfície), e extrapolar para o envelope do multirrotor
% (V com p95 de 1,3 m/s no treino).
%
% As três fontes comparadas:
%   AVL          — painéis, o que era usado até agora como estimativa
%   voo asa fixa — identificado nos doublets (11_fixed_wing/fw_longitudinal.m
%                  e fw_lateral.m), mesma inércia do CAD usada no multirrotor
%   voo quad     — c_p, c_q, c_r da identificação oficial (outputs/P_identified.mat)
%
% Uso:  >> fw_vs_quad_damping
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
rho = p.rho; S = p.wing.S; b = p.wing.b; c = p.wing.c;
Jx = p.J.Jx; Jy = p.J.Jy; Jz = p.J.Jz;

L = load(fullfile(paths.outputs,'fw_longitudinal.mat'));
T = load(fullfile(paths.outputs,'fw_lateral.mat'));
gi = @(S_,nm) mean(S_.Pid(strcmp(S_.pn,nm),:));
gs = @(S_,nm)  std(S_.Pid(strcmp(S_.pn,nm),:));

FW  = [gi(T,'Clp'), gi(L,'Cmq'), gi(T,'Cnr')];        % identificado em voo
FWs = [gs(T,'Clp'), gs(L,'Cmq'), gs(T,'Cnr')];
AVL = [-0.4060, -8.9593, -0.0696];
Jv  = [Jx, Jy, Jz];  Lref = [b, c, b];                % braço da adimensionalização

% ganho por m/s:  ¼ρS·(b² ou c̄²)·C/J   [1/s por (m/s)]
kV = @(C) 0.25*rho*S*(Lref.^2).*C ./ Jv;
gFW = kV(FW);  gAVL = kV(AVL);

Q = load(fullfile(paths.outputs,'P_identified.mat'));
cq = Q.P_final(13:15)';                                % c_p, c_q, c_r do multirrotor
V_TR = 1.29;   V_VAL = 2.3;                            % p95 de treino e de validação
Vfw = 17.0;

nm = {'rolagem (C_lp)','arfagem (C_mq)','guinada (C_nr)'};
fprintf('\n  ================================================================================\n');
fprintf('   Amortecimento aerodinâmico da estrutura: AVL contra voo de asa fixa\n');
fprintf('  ================================================================================\n');
fprintf('  %-16s %10s %10s %8s | %s\n','eixo','AVL','voo','voo/AVL','1/s por m/s: AVL / voo');
for i = 1:3
    fprintf('  %-16s %10.4f %10.4f %8.2f | %8.4f / %8.4f\n', nm{i}, AVL(i), FW(i), FW(i)/AVL(i), gAVL(i), gFW(i));
end
fprintf('  (voo: média de 4 doublets, desvio %.3f, %.2f e %.3f)\n', FWs);

fprintf('\n  ================================================================================\n');
fprintf('   Quanto isso vale no envelope do multirrotor\n');
fprintf('  ================================================================================\n');
fprintf('  %-10s %12s %12s %12s %12s\n','eixo','c identif.','estrut. AVL','estrut. voo','voo / c identif.');
for i = 1:3
    fprintf('  %-10s %12.3f %12.3f %12.3f %11.1f%%\n', nm{i}(1:7), cq(i), ...
        abs(gAVL(i))*V_TR, abs(gFW(i))*V_TR, 100*abs(gFW(i))*V_TR/cq(i));
end
fprintf('  (em V = %.2f m/s, p95 do trecho de treino; c identificados em outputs/P_identified.mat)\n', V_TR);
fprintf('\n  Leitura direta, sem passar por inércia: em asa fixa a 17 m/s a constante de\n');
fprintf('  tempo de rolagem medida é %.2f s. No multirrotor pairado o modelo identificado\n', -1/(gFW(1)*Vfw));
fprintf('  tem %.2f s. Como o efeito aerodinâmico cresce com V, a %.2f m/s ele daria %.1f s,\n', ...
    1/cq(1), V_TR, -1/(gFW(1)*V_TR));
fprintf('  ou seja %.0f vezes mais lento que o observado.\n', (1/abs(gFW(1)*V_TR))/(1/cq(1)));

%% figura
f = figure('Position',[40 40 1250 460],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
Vg = linspace(0, 20, 200)';
lab = {'p','q','r'};
for i = 1:3
    nexttile; hold on; grid on;
    plot(Vg, abs(gAVL(i))*Vg, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',1.6);
    plot(Vg, abs(gFW(i))*Vg, '-',  'Color',[0.85 0.37 0.01], 'LineWidth',2);
    plot(Vfw, abs(gFW(i))*Vfw, 'o', 'MarkerSize',9, 'MarkerFaceColor',[0.85 0.37 0.01], 'MarkerEdgeColor','k');
    yline(cq(i), '-', sprintf('c_%s = %.2f (multirrotor)', lab{i}, cq(i)), ...
        'Color',[0 0.45 0.7], 'LineWidth',2, 'LabelHorizontalAlignment','left');
    xregion(0, V_VAL, 'FaceColor',[0 0.45 0.7], 'FaceAlpha',0.10);
    xlabel('V [m/s]'); ylabel('amortecimento [1/s]');
    title(nm{i});
    if i == 1
        legend({'AVL (painéis)','voo de asa fixa (identificado)','doublets a 17 m/s'}, ...
            'Location','northwest','FontSize',8);
    end
    text(0.03, 0.62, sprintf('envelope\ndo multirrotor\n(V \\leq %.1f m/s)', V_VAL), ...
        'Units','normalized','FontSize',8,'Color',[0 0.35 0.55]);
end
title(tl, 'Amortecimento aerodinâmico da estrutura, medido em asa fixa, contra o amortecimento identificado em multirrotor');
fn = fullfile(paths.images,'new_flights','fw_vs_quad_damping.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',140);
fprintf('\n  Figura: %s\n', fn);
