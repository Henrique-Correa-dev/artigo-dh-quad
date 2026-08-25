% fw_compare_models.m — modelos de asa fixa de referência contra o identificado
% =========================================================================
% Compara os modelos aerodinâmicos de asa fixa disponíveis para o DH com o que
% foi identificado nos doublets de dez/2024, simulando todos nas MESMAS janelas
% de validação, todas FORA dos trechos usados na identificação.
%
% MODELOS COMPARADOS (cada um com a sua geometria de referência própria)
%   AVL DH (dh_st)        S 0,270  c̄ 0,226  b 1,20   identification/2_model/avl
%   AVL VTOL ITA          S 0,350  c̄ 0,250  b 1,40   AVL Stuff/SIMULAÇÃO (workshop)
%   AVL M1, AVL M2        S 0,350  c̄ 0,250  b 1,40   AVL Stuff/coef_M1.m, coef_M2.m
%   FINEP_2               S 0,270  c̄ 0,226  b 1,20   relatório FINEP 2, OEM em voo do DH
%   Baby Shark NTNU       S 0,662  c̄ 0,242  b 2,50   Græsdal (2021), identificado em voo
%
% FINEP_2: modelo completo do Anexo B do relatório infos_Rela_FINEP_2.pdf, união
% do modelo 4 longitudinal com o modelo 2 látero-direcional, identificado por OEM
% em voo do próprio DH. A geometria de referência não aparece em texto no
% relatório (a Tabela 01 é imagem), e assume-se a do DH: o C_Lα deles, 2,878,
% cai a 2% do nosso, 2,938, o que só fecha se a normalização for a mesma.
%   Identificado (voo)    S 0,270  c̄ 0,226  b 1,20   fw_longitudinal.m, fw_lateral.m
%
% NORMALIZAÇÃO. Cada fonte adimensionaliza pela sua própria referência, então os
% coeficientes NÃO são comparáveis diretamente. Todos são convertidos para a
% geometria do DH pela igualdade da força ou do momento dimensional:
%   força  ∝ S·C            → C ← C·(S'/S)
%   momento estático ∝ S·ℓ·C → C ← C·(S'ℓ'/Sℓ)      ℓ = c̄ (long) ou b (lat)
%   momento de taxa  ∝ S·ℓ²·C → C ← C·(S'ℓ'²/Sℓ²)
% Massa e inércia são SEMPRE as do DH: o que se troca é só a aerodinâmica.
%
% DERIVADAS DE CONTROLE. São as NOSSAS em todos os modelos. O curso mecânico das
% superfícies não está no log, cada fonte usa uma normalização diferente, e o
% servo do DH tem a convenção de sinal invertida (ver fw_longitudinal.m). Manter
% o controle fixo isola o que está em comparação: as derivadas de ESTABILIDADE.
%
% Uso:  >> fw_compare_models
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
S = proj.wing.S;  c = proj.wing.c;  b = proj.wing.b;  rho = proj.rho;  g = proj.g;  m = proj.m;
Jx = proj.J.Jx; Jy = proj.J.Jy; Jz = proj.J.Jz; Jxz = proj.J.Jxz;
DIR = fullfile(getenv('HOME'),'Desktop','DH_asafixa_modelos');
if ~exist(DIR,'dir'), mkdir(DIR); end

% janelas de identificação (para conferir que a validação não as toca)
IDL = [306.5 310.5; 326.5 330.5; 344.0 348.0; 367.5 371.5];   % longitudinal
IDT = [303.5 307.5; 312.5 316.5; 322.5 326.5; 333.0 337.0];   % látero-direcional
% janelas de VALIDAÇÃO, escolhidas por atividade de comando fora das de cima
VLL = [312 318; 334 340; 352 358; 360 366; 376 382; 396 402];
VLT = [340 346; 356 362; 380 386; 408 414; 432 438; 446 452];

%% ---------------------------------------------------------------- modelos
% [nome, S', c', b', struct de coeficientes], já com a conversão feita adiante
LON = {};  LAT = {};
gDH  = [0.270 0.226 1.20];  gITA = [0.350 0.250 1.40];  gBS = [0.6617 0.242 2.50];

LON{end+1} = mk_lon('AVL DH (dh\_st)',   gDH,  0.350,  4.585,  8.112,  0.000, -0.7377,  -8.959);
LON{end+1} = mk_lon('AVL VTOL ITA',      gITA, 0.502,  4.635, 11.172, -0.148, -2.2349, -13.335);
LON{end+1} = mk_lon('AVL M1',            gITA, 1.0534, 4.4384, 16.479, -0.0445, -0.2036, -9.5431);
LON{end+1} = mk_lon('AVL M2',            gITA, 1.0533, 4.2169, 18.951, -0.0434, -0.1647, -11.440);
LON{end+1} = mk_lon('FINEP 2 (OEM, voo)', gDH, 0.1936, 2.8784, 31.102,  0.0830, -0.8819,  -9.929);
LON{end+1} = mk_lon('IDENTIFICADO (voo)',gDH,  0.337,  2.938,   8.112,  0.0129, -0.4348,  -5.540);
% O Baby Shark (Græsdal 2021) fica FORA da simulação: é a mesma família de
% aeronave, mas com S 2,5× e b 2,1× maiores, e m e J 6× e 13× maiores. Rodar a
% aerodinâmica dele com a inércia do DH mede a diferença de porte, não a
% qualidade do modelo. Entra nas tabelas de coeficientes e de derivadas
% dimensionais, que são as comparações válidas entre aeronaves diferentes.
BS_lon = mk_lon('Baby Shark NTNU',   gBS,  0.4606, 5.3253,  0.000,  0.0950, -1.4947, -13.140);
BS_lat = mk_lat('Baby Shark NTNU',   gBS,  -0.7310, -0.0354, -0.2419,  0.0759, -0.0752,  1.0778, 0.0000,  0.0953, -0.0823);

LAT{end+1} = mk_lat('AVL DH (dh\_st)',   gDH,  -0.1934, -0.0615, -0.4060,  0.0656, -0.0696,  0.1121, 0.1540,  0.1420, -0.0319);
LAT{end+1} = mk_lat('AVL VTOL ITA',      gITA, -0.2061, -0.0522, -0.4200,  0.1069, -0.1248,  0.1159, 0.2271,  0.1483, -0.0459);
LAT{end+1} = mk_lat('AVL M1 e M2',       gITA, -0.1222, -0.0408, -0.4954,  0.0451, -0.0437,  0.0694, 0.1046,  0.1349, -0.0302);
LAT{end+1} = mk_lat('FINEP 2 (OEM, voo)',gDH, -0.3990, -0.2773, -1.3287,  0.1535, -0.4611, -0.6933, 0.0366,  1.0043, -0.1084);
LAT{end+1} = mk_lat('IDENTIFICADO (voo)',gDH,  -0.1964, -0.0222, -0.0769,  0.0472, -0.2037,  0.1121, 0.1540,  0.1420, -0.0319);

% derivadas de controle: as NOSSAS, iguais em todos os modelos
UE = struct('CLde', 0.007, 'Cmde', 0.150);
UA = struct('CYda', 0.0562, 'Clda', -0.0213, 'Cnda', 0.0042);

% ---- tabelas de coeficientes, montadas ANTES da conversão -----------------
AL = [LON(:); {BS_lon}];  AT = [LAT(:); {BS_lat}];
fprintf('\n  ======== LONGITUDINAL: coeficientes crus (cada fonte na SUA referência)\n');
fprintf('  %-22s %6s %6s | %7s %7s %7s %8s\n', 'modelo','S','c̄','C_L0','C_Lα','C_mα','C_mq');
for k = 1:numel(AL)
    A = AL{k};
    fprintf('  %-22s %6.3f %6.3f | %7.3f %7.3f %7.3f %8.2f\n', lbl(A.nome), A.S, A.c, A.CL0, A.CLa, A.Cma, A.Cmq);
end
fprintf('\n           mesmos, convertidos para a geometria do DH (S %.3f, c̄ %.3f)\n', S, c);
for k = 1:numel(AL)
    A = conv_lon(AL{k}, S, c);
    fprintf('  %-22s %6s %6s | %7.3f %7.3f %7.3f %8.2f\n', lbl(A.nome), '', '', A.CL0, A.CLa, A.Cma, A.Cmq);
end
fprintf('\n  ======== LÁTERO-DIRECIONAL: coeficientes crus\n');
fprintf('  %-22s %6s %6s | %8s %8s %8s %8s %8s\n', 'modelo','S','b','C_Yβ','C_lβ','C_lp','C_nβ','C_nr');
for k = 1:numel(AT)
    B = AT{k};
    fprintf('  %-22s %6.3f %6.2f | %8.4f %8.4f %8.4f %8.4f %8.4f\n', lbl(B.nome), B.S, B.b, B.CYb, B.Clb, B.Clp, B.Cnb, B.Cnr);
end
fprintf('\n           mesmos, convertidos para a geometria do DH (S %.3f, b %.2f)\n', S, b);
for k = 1:numel(AT)
    B = conv_lat(AT{k}, S, b);
    fprintf('  %-22s %6s %6s | %8.4f %8.4f %8.4f %8.4f %8.4f\n', lbl(B.nome), '', '', B.CYb, B.Clb, B.Clp, B.Cnb, B.Cnr);
end
fprintf('\n  ======== DERIVADAS DIMENSIONAIS, cada modelo na SUA aeronave e V de cruzeiro\n');
fprintf('  (comparação válida entre aeronaves de porte diferente: dá a constante de tempo)\n');
fprintf('  %-22s %6s | %10s %12s | %10s %12s\n','modelo','V','M_q [1/s]','tau_pitch [s]','L_p [1/s]','tau_roll [s]');
AER = {'AVL DH (dh_st)',17.0,Jx,Jy; 'AVL VTOL ITA',17.0,Jx,Jy; 'AVL M1',17.0,Jx,Jy; ...
       'AVL M2',17.0,Jx,Jy; 'FINEP_2',17.0,Jx,Jy; 'IDENTIFICADO (voo)',17.0,Jx,Jy; ...
       'Baby Shark NTNU',21.0,NaN,1.0664};
LATOF = [1 2 3 3 4 5 6];    % qual entrada de AT corresponde a cada entrada de AL
for k = 1:numel(AL)
    A = AL{k};  B = AT{LATOF(k)};  V = AER{k,2};  Jxk = AER{k,3};  Jyk = AER{k,4};
    Mq = 0.25*rho*A.S*A.c^2*A.Cmq*V/Jyk;
    Lp = 0.25*rho*B.S*B.b^2*B.Clp*V/Jxk;
    fprintf('  %-22s %6.1f | %10.2f %12.3f | %10.2f %12.3f\n', lbl(A.nome), V, Mq, -1/Mq, Lp, -1/Lp);
end

% conversão de geometria (a partir daqui os modelos vivem na geometria do DH)
for k = 1:numel(LON), LON{k} = conv_lon(LON{k}, S, c); end
for k = 1:numel(LAT), LAT{k} = conv_lat(LAT{k}, S, b); end

%% ---------------------------------------------------------------- longitudinal
F = fw_load('long');
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
TIC = @(y,yh) sqrt(mean((y-yh).^2))/(sqrt(mean(y.^2))+sqrt(mean(yh.^2)));
nL = numel(LON);  RL = nan(size(VLL,1), nL, 2);  TL = nan(size(VLL,1), nL, 2);
fprintf('\n  ================ LONGITUDINAL, janelas de validação (identificação em %s)\n', mat2str(IDL(:,1)'));
for w = 1:size(VLL,1)
    tw = VLL(w,:);
    if any(tw(1) < IDL(:,2) & tw(2) > IDL(:,1)), warning('janela %d toca a identificação', w); end
    ii = F.t >= tw(1) & F.t <= tw(2);
    d = struct('t',F.t(ii),'dt',F.dt,'V',max(F.V(ii),5),'th',F.theta(ii),'ph',F.phi(ii), ...
               'ue',F.de(ii),'q',F.q(ii),'az',F.az(ii));
    Y = cell(nL,1);
    for k = 1:nL
        P = [LON{k}.CL0; LON{k}.CLa; UE.CLde; LON{k}.Cm0; LON{k}.Cma; LON{k}.Cmq; UE.Cmde];
        Y{k} = sim_sp(P, d, S, c, rho, g, m, Jy, LON{k}.CLq);
        RL(w,k,:) = [R2(d.q,Y{k}.q), R2(d.az,Y{k}.az)];
        TL(w,k,:) = [TIC(d.q,Y{k}.q), TIC(d.az,Y{k}.az)];
    end
    fig_lon(d, Y, LON, RL(w,:,:), tw, DIR, mean(d.V));
end
tabela('LONGITUDINAL', {'q','a_z'}, VLL, LON, RL);

%% ---------------------------------------------------------------- lateral
F = fw_load('lat');
gam = Jx*Jz - Jxz^2;
G = struct('g1',Jxz*(Jx-Jy+Jz)/gam, 'g2',(Jz*(Jz-Jy)+Jxz^2)/gam, 'g3',Jz/gam, ...
           'g4',Jxz/gam, 'g7',(Jx*(Jx-Jy)+Jxz^2)/gam, 'g8',Jx/gam);
nT = numel(LAT);  RT = nan(size(VLT,1), nT, 3);
fprintf('\n  ================ LÁTERO-DIRECIONAL, janelas de validação\n');
for w = 1:size(VLT,1)
    tw = VLT(w,:);
    if any(tw(1) < IDT(:,2) & tw(2) > IDT(:,1)), warning('janela %d toca a identificação', w); end
    ii = F.t >= tw(1) & F.t <= tw(2);
    d = struct('t',F.t(ii),'dt',F.dt,'V',max(F.V(ii),5),'th',F.theta(ii),'ph',F.phi(ii), ...
               'q',F.q(ii),'ua',F.da(ii),'p',F.p(ii),'r',F.r(ii),'ay',F.ay(ii));
    Y = cell(nT,1);
    for k = 1:nT
        P = [LAT{k}.CYb; UA.CYda; LAT{k}.Clb; LAT{k}.Clp; UA.Clda; LAT{k}.Cnb; LAT{k}.Cnr; UA.Cnda];
        FIX = struct('CYp',LAT{k}.CYp,'CYr',LAT{k}.CYr,'Clr',LAT{k}.Clr,'Cnp',LAT{k}.Cnp);
        Y{k} = sim_lat(P, d, S, b, rho, g, m, G, FIX);
        RT(w,k,:) = [R2(d.p,Y{k}.p), R2(d.r,Y{k}.r), R2(d.ay,Y{k}.ay)];
    end
    fig_lat(d, Y, LAT, RT(w,:,:), tw, DIR, mean(d.V));
end
tabela('LÁTERO-DIRECIONAL', {'p','r','a_y'}, VLT, LAT, RT);
fprintf('\n  Figuras em %s\n', DIR);

%% ================================================================= helpers
function M = mk_lon(nome, geo, CL0, CLa, CLq, Cm0, Cma, Cmq)
    M = struct('nome',nome,'S',geo(1),'c',geo(2),'b',geo(3), ...
               'CL0',CL0,'CLa',CLa,'CLq',CLq,'Cm0',Cm0,'Cma',Cma,'Cmq',Cmq);
end
function M = mk_lat(nome, geo, CYb, Clb, Clp, Cnb, Cnr, CYp, CYr, Clr, Cnp)
    M = struct('nome',nome,'S',geo(1),'c',geo(2),'b',geo(3), ...
               'CYb',CYb,'Clb',Clb,'Clp',Clp,'Cnb',Cnb,'Cnr',Cnr, ...
               'CYp',CYp,'CYr',CYr,'Clr',Clr,'Cnp',Cnp);
end
function M = conv_lon(M, S, c)
    kF = M.S/S;  kM = (M.S*M.c)/(S*c);  kM2 = (M.S*M.c^2)/(S*c^2);
    M.CL0 = M.CL0*kF;  M.CLa = M.CLa*kF;  M.CLq = M.CLq*kM;
    M.Cm0 = M.Cm0*kM;  M.Cma = M.Cma*kM;  M.Cmq = M.Cmq*kM2;
end
function M = conv_lat(M, S, b)
    kF = M.S/S;  kFb = (M.S*M.b)/(S*b);  kM = (M.S*M.b)/(S*b);  kM2 = (M.S*M.b^2)/(S*b^2);
    M.CYb = M.CYb*kF;   M.CYp = M.CYp*kFb;  M.CYr = M.CYr*kFb;
    M.Clb = M.Clb*kM;   M.Cnb = M.Cnb*kM;
    M.Clp = M.Clp*kM2;  M.Clr = M.Clr*kM2;  M.Cnp = M.Cnp*kM2;  M.Cnr = M.Cnr*kM2;
end
function tabela(titulo, canais, W, MODS, R)
    fprintf('\n  ---- %s: R² por janela\n', titulo);
    for ch = 1:numel(canais)
        fprintf('\n   %s:  %-11s', canais{ch}, 'janela');
        for k = 1:numel(MODS), fprintf(' %14s', trunca(strrep(MODS{k}.nome,'\',''),14)); end
        fprintf('\n');
        for w = 1:size(W,1)
            fprintf('        %4.0f–%-6.0f', W(w,1), W(w,2));
            for k = 1:numel(MODS), fprintf(' %14.3f', R(w,k,ch)); end
            fprintf('\n');
        end
        fprintf('        %-11s', 'MÉDIA');
        for k = 1:numel(MODS), fprintf(' %14.3f', mean(R(:,k,ch),'omitnan')); end
        fprintf('\n');
    end
end
function s = lbl(s), s = strrep(s,'\_','_'); end
function s = trunca(s, n), if numel(s) > n, s = s(1:n); end, end
function [C, LW] = cores(n)
%CORES  Paleta: modelos de referência em tom fraco e linha fina, o modelo
%   identificado neste trabalho (SEMPRE o último) em laranja forte e grosso.
    ref = [0.55 0.68 0.82; 0.72 0.82 0.66; 0.82 0.70 0.86; 0.95 0.85 0.55; 0.70 0.75 0.78];
    nref = n - 1;
    C = [ref(1:nref,:); 0.85 0.37 0.01];
    LW = [repmat(1.0, nref, 1); 2.0];
end
function fig_lon(d, Y, MODS, R, tw, DIR, Vm)
    n = numel(MODS);  [C, LW] = cores(n);
    f = figure('Position',[30 30 1250 830],'Color','w','Visible','off');
    try, f.Theme = 'light'; catch, end
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    nexttile; hold on; grid on;
    for k = 1:n, plot(d.t, Y{k}.q, '-', 'Color',C(k,:), 'LineWidth',LW(k)); end
    plot(d.t, d.q, 'k-', 'LineWidth',2.2);
    ylabel('q [rad/s]'); xlim(tw);
    lg = cell(n+1,1);
    for k = 1:n, lg{k} = sprintf('%s  (R² %.2f)', MODS{k}.nome, R(1,k,1)); end
    lg{n+1} = 'medido';
    legend(lg, 'Location','eastoutside','FontSize',8);
    nexttile; hold on; grid on;
    for k = 1:n, plot(d.t, Y{k}.az, '-', 'Color',C(k,:), 'LineWidth',LW(k)); end
    plot(d.t, d.az, 'k-', 'LineWidth',2.2);
    ylabel('a_z [m/s²]'); xlim(tw);
    lg2 = cell(n+1,1);
    for k = 1:n, lg2{k} = sprintf('R² %.2f', R(1,k,2)); end
    lg2{n+1} = 'medido';
    legend(lg2, 'Location','eastoutside','FontSize',8);
    nexttile; hold on; grid on;
    plot(d.t, d.ue, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1.4);
    ylabel('\delta_e normalizado'); xlabel('t [s]'); xlim(tw); ylim([-1 1]);
    title(tl, sprintf(['Asa fixa, LONGITUDINAL, janela de validação %.0f a %.0f s (V %.1f m/s, voo 17/12/2024)\n' ...
        'modelos de referência contra o identificado, derivadas de controle iguais em todos'], tw, Vm), 'FontWeight','bold');
    exportgraphics(f, fullfile(DIR, sprintf('long_%04.0f-%04.0f.png', tw)), 'BackgroundColor','white','Resolution',130);
    close(f);
end
function fig_lat(d, Y, MODS, R, tw, DIR, Vm)
    n = numel(MODS);  [C, LW] = cores(n);
    f = figure('Position',[30 30 1250 980],'Color','w','Visible','off');
    try, f.Theme = 'light'; catch, end
    tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
    ch = {'p','r','ay'};  rot = {'p [rad/s]','r [rad/s]','a_y [m/s²]'};
    for j = 1:3
        nexttile; hold on; grid on;
        for k = 1:n, plot(d.t, Y{k}.(ch{j}), '-', 'Color',C(k,:), 'LineWidth',LW(k)); end
        plot(d.t, d.(ch{j}), 'k-', 'LineWidth',2.2);
        ylabel(rot{j}); xlim(tw);
        lg = cell(n+1,1);
        for k = 1:n
            if j == 1, lg{k} = sprintf('%s  (R² %.2f)', MODS{k}.nome, R(1,k,j));
            else,      lg{k} = sprintf('R² %.2f', R(1,k,j)); end
        end
        lg{n+1} = 'medido';
        legend(lg, 'Location','eastoutside','FontSize',8);
    end
    nexttile; hold on; grid on;
    plot(d.t, d.ua, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1.4);
    ylabel('\delta_a normalizado'); xlabel('t [s]'); xlim(tw); ylim([-1 1]);
    title(tl, sprintf(['Asa fixa, LÁTERO-DIRECIONAL, janela de validação %.0f a %.0f s (V %.1f m/s, voo 17/12/2024)\n' ...
        'modelos de referência contra o identificado, derivadas de controle iguais em todos'], tw, Vm), 'FontWeight','bold');
    exportgraphics(f, fullfile(DIR, sprintf('lat_%04.0f-%04.0f.png', tw)), 'BackgroundColor','white','Resolution',130);
    close(f);
end
function y = sim_sp(P, d, S, c, rho, g, m, Jy, CLq)
    CL0=P(1); CLa=P(2); CLde=P(3);
    N = numel(d.t);  nsub = 5;  h = d.dt/nsub;
    al = zeros(N,1); q = zeros(N,1);
    q(1) = d.q(1);
    qb1 = 0.5*rho*d.V(1)^2;  CLt = -d.az(1)*m/(qb1*S);
    al(1) = max(min((CLt - CL0 - CLde*d.ue(1))/max(CLa,0.5), 0.4), -0.4);
    for k = 1:N-1
        a = al(k); qq = q(k);
        for s = 1:nsub
            f = (s-1)/nsub;  fh = f + 0.5/nsub;  f1 = f + 1/nsub;
            [k1a,k1q] = dyn_l(a,         qq,         k, f,  P, d, S, c, rho, g, m, Jy, CLq);
            [k2a,k2q] = dyn_l(a+h/2*k1a, qq+h/2*k1q, k, fh, P, d, S, c, rho, g, m, Jy, CLq);
            [k3a,k3q] = dyn_l(a+h/2*k2a, qq+h/2*k2q, k, fh, P, d, S, c, rho, g, m, Jy, CLq);
            [k4a,k4q] = dyn_l(a+h*k3a,   qq+h*k3q,   k, f1, P, d, S, c, rho, g, m, Jy, CLq);
            a  = a  + h/6*(k1a + 2*k2a + 2*k3a + k4a);
            qq = qq + h/6*(k1q + 2*k2q + 2*k3q + k4q);
        end
        al(k+1) = max(min(a,0.6),-0.6);  q(k+1) = max(min(qq,20),-20);
    end
    qbar = 0.5*rho*d.V.^2;
    y.q = q;  y.al = al;
    y.az = -(qbar*S/m).*(CL0 + CLa*al + CLq*(q.*c./(2*d.V)) + CLde*d.ue);
end
function [ad, qd] = dyn_l(al, q, k, f, P, d, S, c, rho, g, m, Jy, CLq)
    ip = @(x) x(k) + f*(x(min(k+1,numel(x))) - x(k));
    V = ip(d.V);  th = ip(d.th);  ph = ip(d.ph);  ue = ip(d.ue);
    qbar = 0.5*rho*V^2;  qhat = q*c/(2*V);
    CL = P(1) + P(2)*al + CLq*qhat + P(3)*ue;
    Cm = P(4) + P(5)*al + P(6)*qhat + P(7)*ue;
    ad = q + (g/V)*cos(th)*cos(ph) - (qbar*S/(m*V))*CL;
    qd = (qbar*S*c/Jy)*Cm;
end
function y = sim_lat(P, d, S, b, rho, g, m, G, FIX)
    N = numel(d.t);  nsub = 5;  h = d.dt/nsub;
    be = zeros(N,1); p = zeros(N,1); r = zeros(N,1);
    p(1) = d.p(1);  r(1) = d.r(1);
    qb1 = 0.5*rho*d.V(1)^2;
    be(1) = max(min(d.ay(1)*m/(qb1*S*P(1)), 0.3), -0.3);
    for k = 1:N-1
        bb = be(k); pp = p(k); rr = r(k);
        for s = 1:nsub
            f = (s-1)/nsub;  fh = f + 0.5/nsub;  f1 = f + 1/nsub;
            [a1,a2,a3] = dyn_t(bb,pp,rr, k,f,  P,d,S,b,rho,g,m,G,FIX);
            [b1,b2,b3] = dyn_t(bb+h/2*a1, pp+h/2*a2, rr+h/2*a3, k,fh, P,d,S,b,rho,g,m,G,FIX);
            [c1,c2,c3] = dyn_t(bb+h/2*b1, pp+h/2*b2, rr+h/2*b3, k,fh, P,d,S,b,rho,g,m,G,FIX);
            [e1,e2,e3] = dyn_t(bb+h*c1,   pp+h*c2,   rr+h*c3,   k,f1, P,d,S,b,rho,g,m,G,FIX);
            bb = bb + h/6*(a1+2*b1+2*c1+e1);
            pp = pp + h/6*(a2+2*b2+2*c2+e2);
            rr = rr + h/6*(a3+2*b3+2*c3+e3);
        end
        be(k+1) = max(min(bb,0.5),-0.5);  p(k+1) = max(min(pp,20),-20);  r(k+1) = max(min(rr,20),-20);
    end
    qbar = 0.5*rho*d.V.^2;
    y.p = p;  y.r = r;  y.be = be;
    y.ay = (qbar*S/m).*(P(1)*be + FIX.CYp*(p*b./(2*d.V)) + FIX.CYr*(r*b./(2*d.V)) + P(2)*d.ua);
end
function [bd, pd, rd] = dyn_t(be, p, r, k, f, P, d, S, b, rho, g, m, G, FIX)
    ip = @(x) x(k) + f*(x(min(k+1,numel(x))) - x(k));
    V = ip(d.V);  th = ip(d.th);  ph = ip(d.ph);  q = ip(d.q);  ua = ip(d.ua);
    qbar = 0.5*rho*V^2;  ph_ = p*b/(2*V);  rh_ = r*b/(2*V);
    CY = P(1)*be + FIX.CYp*ph_ + FIX.CYr*rh_ + P(2)*ua;
    Cl = P(3)*be + P(4)*ph_ + FIX.Clr*rh_ + P(5)*ua;
    Cn = P(6)*be + FIX.Cnp*ph_ + P(7)*rh_ + P(8)*ua;
    L = qbar*S*b*Cl;  Nn = qbar*S*b*Cn;
    bd = (qbar*S/(m*V))*CY - r + (g/V)*cos(th)*sin(ph);
    pd = G.g1*p*q - G.g2*q*r + G.g3*L + G.g4*Nn;
    rd = G.g7*p*q - G.g1*q*r + G.g4*L + G.g8*Nn;
end
