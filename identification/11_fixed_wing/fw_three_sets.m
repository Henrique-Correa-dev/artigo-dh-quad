% fw_three_sets.m — os três conjuntos de coeficientes de asa fixa do DH
% =========================================================================
% Sato, Ana e o identificado neste trabalho, lado a lado, em três formas:
%
%   (a) CRU, como cada um publicou. Os três usam a MESMA geometria de
%       referência (S 0,27 m², c̄ 0,226 m, b 1,20 m), então as derivadas
%       ESTÁTICAS (C_L0, C_Lα, C_mα, C_Yβ, C_lβ, C_nβ) já são comparáveis.
%
%   (b) DERIVADAS DE TAXA na mesma convenção. Não são comparáveis cruas,
%       porque cada um normaliza a taxa de um jeito:
%         Sato   p̂ = p·b/(2·V_ref),  V_ref = 12 m/s fixo
%         Ana    p̂ = p·b/V_ref,      V_ref = 15 m/s fixo, sem o fator 2
%         nosso  p̂ = p·b/(2·V),      V instantânea
%       Igualando o momento dimensional em V = 17 m/s, o coeficiente
%       equivalente na convenção usual (p·b/2V) é
%         Sato   C ← C·(V/V_ref)          = C·1,417
%         Ana    C ← C·(2V/V_ref)         = C·2,267
%
%   (c) DIMENSIONAL [1/s], já com a inércia de cada um. É o que a dinâmica
%       enxerga, e a única forma em que as três colunas significam a mesma
%       coisa física. Os dois usam I_xx 3,3 vezes o nosso.
%
% Uso:  >> fw_three_sets
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
S = proj.wing.S;  c = proj.wing.c;  b = proj.wing.b;  rho = proj.rho;
V = 17;  qb = 0.5*rho*V^2;

N = {'Sato','Ana (FINEP 2)','IDENTIFICADO'};
% massa e inércia
GEN = [2.200 0.14410 0.11550 0.25716; 2.300 0.14410 0.11550 0.25720; proj.m proj.J.Jx proj.J.Jy proj.J.Jz];
% fator de conversão da convenção de taxa para p·b/(2V) em V = 17
KV  = [V/12, 2*V/15, 1];
% estáticas: C_L0 C_D0 C_Lα C_m0 C_mα | C_Yβ C_lβ C_nβ
EST = [0.361632 0.0813815 6.98729 0.108861 -1.34314   0.609958 -0.1036037 0.127030; ...
       0.193615 0.1004    2.87836 0.083009 -0.881949 -0.398966 -0.277332  0.153474; ...
       0.337    0.05      2.938   0.012900 -0.434800 -0.196400 -0.022200  0.047200];
est_n = {'C_L0','C_D0','C_Lα','C_m0','C_mα','C_Yβ','C_lβ','C_nβ'};
% de taxa (cruas): C_Lq C_mq | C_Yp C_Yr C_lp C_lr C_np C_nr
TAX = [15.4571 -40.2966   5.00530 11.731150 -0.785502 0.433016 -0.326726 -0.781478; ...
       31.1023  -9.92909 -0.693329 0.0366157 -1.32871 1.00433  -0.108353 -0.461062; ...
        8.112   -5.540    0.112100 0.154000  -0.076900 0.142000 -0.031900 -0.203700];
tax_n = {'C_Lq','C_mq','C_Yp','C_Yr','C_lp','C_lr','C_np','C_nr'};

fprintf('\n  =============== MASSA E INÉRCIA\n');
fprintf('  %-16s %8s %9s %9s %9s\n','','m [kg]','J_x','J_y','J_z');
for k = 1:3, fprintf('  %-16s %8.3f %9.4f %9.4f %9.4f\n', N{k}, GEN(k,:)); end
fprintf('  %-16s %8.2f %9.2f %9.2f %9.2f   (Sato/Ana ÷ nosso)\n','razão', GEN(1,:)./GEN(3,:));

fprintf('\n  =============== (a) DERIVADAS ESTÁTICAS, diretamente comparáveis\n');
fprintf('  %-8s %12s %14s %14s %10s\n','coef',N{1},N{2},N{3},'Sato/Ana?');
for i = 1:numel(est_n)
    fprintf('  %-8s %12.4f %14.4f %14.4f %10s\n', est_n{i}, EST(1,i), EST(2,i), EST(3,i), sinal(EST(:,i)));
end

fprintf('\n  =============== (b) DERIVADAS DE TAXA na convenção p·b/(2V), em V = %.0f m/s\n', V);
fprintf('  (fatores de conversão: Sato ×%.3f, Ana ×%.3f)\n', KV(1), KV(2));
fprintf('  %-8s %12s %14s %14s | %9s %9s\n','coef',N{1},N{2},N{3},'S/nosso','A/nosso');
TAXC = TAX .* KV(:);
for i = 1:numel(tax_n)
    fprintf('  %-8s %12.4f %14.4f %14.4f | %9.1f %9.1f\n', tax_n{i}, TAXC(1,i), TAXC(2,i), TAXC(3,i), ...
        TAXC(1,i)/TAXC(3,i), TAXC(2,i)/TAXC(3,i));
end

fprintf('\n  =============== (c) DIMENSIONAL [1/s] em V = %.0f m/s, com a inércia de cada um\n', V);
fprintf('  %-14s %11s %11s %11s | %10s %10s %10s\n','', 'M_q','L_p','N_r','tau_arf','tau_rol','tau_gui');
DIM = nan(3,3);
for k = 1:3
    DIM(k,1) = qb*S*c/GEN(k,3) * TAXC(k,2) * c/(2*V);
    DIM(k,2) = qb*S*b/GEN(k,2) * TAXC(k,5) * b/(2*V);
    DIM(k,3) = qb*S*b/GEN(k,4) * TAXC(k,8) * b/(2*V);
    fprintf('  %-14s %11.2f %11.2f %11.2f | %10.3f %10.3f %10.3f\n', N{k}, DIM(k,:), -1./DIM(k,:));
end
fprintf('  %-14s %11.1f %11.1f %11.1f\n','Sato / nosso', DIM(1,:)./DIM(3,:));
fprintf('  %-14s %11.1f %11.1f %11.1f\n','Ana / nosso',  DIM(2,:)./DIM(3,:));

%% figura
CC = [0.55 0.68 0.82; 0.72 0.82 0.66; 0.85 0.37 0.01];
f = figure('Position',[30 30 1400 480],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hb = bar(categorical(est_n, est_n), EST'); for k=1:3, hb(k).FaceColor = CC(k,:); end
grid on; ylabel('valor'); title('(a) derivadas estáticas'); legend(N,'Location','southwest','FontSize',8);
nexttile; hb = bar(categorical(tax_n, tax_n), TAXC'); for k=1:3, hb(k).FaceColor = CC(k,:); end
grid on; ylabel('valor'); title(sprintf('(b) derivadas de taxa, convenção p·b/2V, V = %.0f m/s', V));
nexttile; hb = bar(categorical({'M_q','L_p','N_r'},{'M_q','L_p','N_r'}), DIM'); for k=1:3, hb(k).FaceColor = CC(k,:); end
grid on; ylabel('[1/s]'); title('(c) dimensional, com a inércia de cada um');
title(tl, sprintf('Três conjuntos de coeficientes de asa fixa do DH, em V = %.0f m/s', V), 'FontWeight','bold');
fn = fullfile(paths.images,'fw_three_sets.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',150);
dd = fullfile(getenv('HOME'),'Desktop','DH_asafixa_modelos');
if ~exist(dd,'dir'), mkdir(dd); end
copyfile(fn, fullfile(dd,'fw_three_sets.png'));
fprintf('\n  Figura: %s  (cópia na Mesa)\n', fn);

function s = sinal(v)
    if all(v > 0) || all(v < 0), s = 'mesmo sinal'; else, s = 'SINAL DIFERE'; end
end
