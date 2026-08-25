% fw_compare_exec.m — modelos executáveis do Sato e da Ana contra o identificado
% =========================================================================
% Compara os três modelos de asa fixa do DH rodando cada um COM A SUA PRÓPRIA
% ESTRUTURA (massa, tensor de inércia, normalização das taxas e velocidade de
% referência), nas mesmas janelas de validação, todas fora dos trechos usados
% na identificação de cada trabalho.
%
% FONTES
%   Sato   AVL Stuff/Modelo executável Sato/.../modelo_DH.m + coef_DH.m
%   Ana    AVL Stuff/Modelo executável Ana/modelo_DH/aerodynamics.m + Coef_DH.mat
%          (é o modelo do relatório infos_Rela_FINEP_2.pdf, Anexo B)
%   nosso  11_fixed_wing/fw_longitudinal.m, fw_lateral.m
%
% AS TRÊS DIFERENÇAS ESTRUTURAIS QUE PRECISAM SER RESPEITADAS
%
%  1. Inércia e massa. Os dois usam o mesmo conjunto, bem diferente do nosso:
%                     m       I_xx      I_yy      I_zz      I_xz
%       Sato         2,20    0,1441    0,1155    0,2572   −0,00167
%       Ana          2,30    0,1441    0,1155    0,2572   −0,0017
%       nosso (CAD)  1,993   0,0432    0,0844    0,1262   +0,00157
%     I_xx deles é 3,3 vezes o nosso. Comparar C_lp sem isso não faz sentido:
%     o que a dinâmica enxerga é C_lp/I_xx, não C_lp.
%
%  2. Normalização das taxas. Nenhum dos dois usa a velocidade instantânea:
%       Sato   q̄ = q·c̄/(2·V_ref)   com V_ref = 12 m/s fixo
%       Ana    q̄ = q·c̄/V_ref       com V_ref = 15 m/s fixo, SEM o fator 2
%       nosso  q̄ = q·c̄/(2·V)       com V do pitot, instantânea
%     O amortecimento efetivo deles cresce com V², e não com V.
%
%  3. Derivadas de controle. As nossas são por deflexão NORMALIZADA em [−1,1]
%     (o curso mecânico não está no log), as deles por radiano. Na comparação
%     principal o controle é o NOSSO em todos, para isolar as derivadas de
%     ESTABILIDADE, que é o que está em questão. A opção 'own_ctrl' roda cada
%     um com o seu próprio conjunto, com δ = u·δ_max e δ_max = 25°.
%
% Uso:  >> fw_compare_exec              controle padronizado
%       >> fw_compare_exec('own_ctrl', true)   cada um com o seu controle
% =========================================================================
function fw_compare_exec(varargin)
opt = struct('own_ctrl', false, 'dmax_deg', 25, 'wins_lon', [], 'wins_lat', [], 'tag', '');
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
S = proj.wing.S;  c = proj.wing.c;  b = proj.wing.b;  rho = proj.rho;  g = proj.g;
DIR = fullfile(getenv('HOME'),'Desktop','DH_asafixa_modelos');
if ~exist(DIR,'dir'), mkdir(DIR); end
dmax = deg2rad(opt.dmax_deg);

%% ---------------------------------------------------------------- modelos
M = {};
% --- Sato: coef_DH.m do Modelo_DH_Completo_SATO_REV00
M{end+1} = struct('nome','Sato', 'm',2.20, 'Jx',0.14410, 'Jy',0.11550, 'Jz',0.25716, 'Jxz',0.00167, ...
    'Vref',12, 'fat',2, ...                       % q̄ = q·c/(fat·Vref)
    'CL0',0.361632, 'CLa',6.98729, 'CLq',15.4571, 'Cm0',0.108861, 'Cma',-1.34314, 'Cmq',-40.2966, ...
    'CYb',0.609958, 'CYp',5.00530, 'CYr',11.73115, 'Clb',-0.1036037, 'Clp',-0.785502, 'Clr',0.433016, ...
    'Cnb',0.127030, 'Cnp',-0.326726, 'Cnr',-0.781478, ...
    'CLde',0, 'Cmde',3.46691, 'CYda',0.227152, 'Clda',-0.278348, 'Cnda',0.127621);
% --- Ana: Coef_DH.mat (= modelo do FINEP 2)
M{end+1} = struct('nome','Ana (FINEP 2)', 'm',2.30, 'Jx',0.14410, 'Jy',0.11550, 'Jz',0.25720, 'Jxz',0.00170, ...
    'Vref',15, 'fat',1, ...                       % q̄ = q·c/Vref, sem o 2
    'CL0',0.193615, 'CLa',2.87836, 'CLq',31.1023, 'Cm0',0.0830091, 'Cma',-0.881949, 'Cmq',-9.92909, ...
    'CYb',-0.398966, 'CYp',-0.693329, 'CYr',0.0366157, 'Clb',-0.277332, 'Clp',-1.32871, 'Clr',1.00433, ...
    'Cnb',0.153474, 'Cnp',-0.108353, 'Cnr',-0.461062, ...
    'CLde',0, 'Cmde',1.04662, 'CYda',0, 'Clda',-0.710799, 'Cnda',0.0395561);
% --- nosso
% Coeficientes longitudinais da identificação com 3 doublets (fw_longitudinal_3win),
% que deixa o doublet de 326,5 a 330,5 s livre para validação. Diferem em menos de
% 4% dos obtidos com os 4 doublets, o que já é um teste de robustez do conjunto.
M{end+1} = struct('nome','IDENTIFICADO (nosso)', 'm',proj.m, 'Jx',proj.J.Jx, 'Jy',proj.J.Jy, ...
    'Jz',proj.J.Jz, 'Jxz',proj.J.Jxz, 'Vref',NaN, 'fat',2, ...   % Vref NaN = usa V instantânea
    'CL0',0.351, 'CLa',2.836, 'CLq',8.112, 'Cm0',0.012, 'Cma',-0.420, 'Cmq',-5.516, ...
    'CYb',-0.1964, 'CYp',0.1121, 'CYr',0.1540, 'Clb',-0.0222, 'Clp',-0.0769, 'Clr',0.1420, ...
    'Cnb',0.0472, 'Cnp',-0.0319, 'Cnr',-0.2037, ...
    'CLde',0.005, 'Cmde',0.147, 'CYda',0.0562, 'Clda',-0.0213, 'Cnda',0.0042);
nM = numel(M);  OUR = nM;

% controle padronizado (o nosso, em deflexão normalizada) salvo se own_ctrl
if ~opt.own_ctrl
    for k = 1:nM-1
        M{k}.CLde = M{OUR}.CLde;  M{k}.Cmde = M{OUR}.Cmde;
        M{k}.CYda = M{OUR}.CYda;  M{k}.Clda = M{OUR}.Clda;  M{k}.Cnda = M{OUR}.Cnda;
        M{k}.uscale = 1;
    end
    M{OUR}.uscale = 1;
else
    for k = 1:nM-1, M{k}.uscale = dmax; end     % deflexão em radianos
    M{OUR}.uscale = 1;
end
for k = 1:nM, M{k} = gammas(M{k}); end

fprintf('\n  ======== ESTRUTURA DE CADA MODELO\n');
fprintf('  %-22s %6s %8s %8s %8s %9s %8s %6s\n','modelo','m','J_x','J_y','J_z','J_xz','V_ref','fator');
for k = 1:nM
    fprintf('  %-22s %6.3f %8.4f %8.4f %8.4f %9.5f %8s %6d\n', M{k}.nome, M{k}.m, ...
        M{k}.Jx, M{k}.Jy, M{k}.Jz, M{k}.Jxz, vstr(M{k}.Vref), M{k}.fat);
end
fprintf('\n  ======== AMORTECIMENTO EFETIVO, já com a normalização e a inércia de cada um, a V = 17 m/s\n');
fprintf('  %-22s %11s %12s | %11s %12s\n','modelo','M_q [1/s]','tau_arf [s]','L_p [1/s]','tau_rol [s]');
V0 = 17;  qb0 = 0.5*rho*V0^2;
for k = 1:nM
    Vr = M{k}.Vref;  if ~isfinite(Vr), Vr = V0; end
    Mq = qb0*S*c/M{k}.Jy * M{k}.Cmq * c/(M{k}.fat*Vr);
    Lp = qb0*S*b/M{k}.Jx * M{k}.Clp * b/(M{k}.fat*Vr);
    fprintf('  %-22s %11.2f %12.3f | %11.2f %12.3f\n', M{k}.nome, Mq, -1/Mq, Lp, -1/Lp);
end

%% ---------------------------------------------------------------- longitudinal
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
F = fw_load('long');
VLL = [312 318; 334 340; 352 358; 360 366; 376 382; 396 402];
if ~isempty(opt.wins_lon), VLL = opt.wins_lon; end
RL = nan(size(VLL,1), nM, 2);
for w = 1:size(VLL,1)
    ii = F.t >= VLL(w,1) & F.t <= VLL(w,2);
    d = struct('t',F.t(ii),'dt',F.dt,'V',max(F.V(ii),5),'th',F.theta(ii),'ph',F.phi(ii), ...
               'ue',F.de(ii),'q',F.q(ii),'az',F.az(ii));
    Y = cell(nM,1);
    for k = 1:nM
        Y{k} = sim_lon(M{k}, d, S, c, rho, g);
        RL(w,k,:) = [R2(d.q,Y{k}.q), R2(d.az,Y{k}.az)];
    end
    fig_lon(d, Y, M, RL(w,:,:), VLL(w,:), DIR, mean(d.V), opt.own_ctrl, opt.tag);
end
tab('LONGITUDINAL', {'q','a_z'}, VLL, M, RL);

%% ---------------------------------------------------------------- lateral
F = fw_load('lat');
VLT = [340 346; 356 362; 380 386; 408 414; 432 438; 446 452];
if ~isempty(opt.wins_lat), VLT = opt.wins_lat; end
RT = nan(size(VLT,1), nM, 3);
for w = 1:size(VLT,1)
    ii = F.t >= VLT(w,1) & F.t <= VLT(w,2);
    d = struct('t',F.t(ii),'dt',F.dt,'V',max(F.V(ii),5),'th',F.theta(ii),'ph',F.phi(ii), ...
               'q',F.q(ii),'ua',F.da(ii),'p',F.p(ii),'r',F.r(ii),'ay',F.ay(ii));
    Y = cell(nM,1);
    for k = 1:nM
        Y{k} = sim_lat(M{k}, d, S, b, rho, g);
        RT(w,k,:) = [R2(d.p,Y{k}.p), R2(d.r,Y{k}.r), R2(d.ay,Y{k}.ay)];
    end
    fig_lat(d, Y, M, RT(w,:,:), VLT(w,:), DIR, mean(d.V), opt.own_ctrl, opt.tag);
end
tab('LÁTERO-DIRECIONAL', {'p','r','a_y'}, VLT, M, RT);
fprintf('\n  Figuras em %s  (prefixo exec_)\n', DIR);
end

%% ================================================================= helpers
function s = vstr(v), if isfinite(v), s = sprintf('%.0f',v); else, s = 'V(t)'; end, end
function M = gammas(M)
    G = M.Jx*M.Jz - M.Jxz^2;
    M.g1 = M.Jxz*(M.Jx-M.Jy+M.Jz)/G;  M.g2 = (M.Jz*(M.Jz-M.Jy)+M.Jxz^2)/G;
    M.g3 = M.Jz/G;  M.g4 = M.Jxz/G;
    M.g7 = (M.Jx*(M.Jx-M.Jy)+M.Jxz^2)/G;  M.g8 = M.Jx/G;
end
function h = hatf(M, V)      % divisor da normalização de taxa
    if isfinite(M.Vref), h = M.fat*M.Vref; else, h = M.fat*V; end
end
function y = sim_lon(M, d, S, c, rho, g)
    N = numel(d.t);  nsub = 5;  h = d.dt/nsub;
    al = zeros(N,1);  q = zeros(N,1);  q(1) = d.q(1);
    qb1 = 0.5*rho*d.V(1)^2;  CLt = -d.az(1)*M.m/(qb1*S);
    al(1) = max(min((CLt - M.CL0 - M.CLde*M.uscale*d.ue(1))/max(M.CLa,0.5), 0.4), -0.4);
    for k = 1:N-1
        a = al(k);  qq = q(k);
        for s = 1:nsub
            f = (s-1)/nsub;  fh = f+0.5/nsub;  f1 = f+1/nsub;
            [k1a,k1q] = dl(a,         qq,         k,f,  M,d,S,c,rho,g);
            [k2a,k2q] = dl(a+h/2*k1a, qq+h/2*k1q, k,fh, M,d,S,c,rho,g);
            [k3a,k3q] = dl(a+h/2*k2a, qq+h/2*k2q, k,fh, M,d,S,c,rho,g);
            [k4a,k4q] = dl(a+h*k3a,   qq+h*k3q,   k,f1, M,d,S,c,rho,g);
            a = a + h/6*(k1a+2*k2a+2*k3a+k4a);  qq = qq + h/6*(k1q+2*k2q+2*k3q+k4q);
        end
        al(k+1) = max(min(a,0.6),-0.6);  q(k+1) = max(min(qq,20),-20);
    end
    qbar = 0.5*rho*d.V.^2;
    CL = M.CL0 + M.CLa*al + M.CLq*(q*c./hatf(M,d.V)) + M.CLde*M.uscale*d.ue;
    y.q = q;  y.al = al;  y.az = -(qbar*S/M.m).*CL;
end
function [ad, qd] = dl(al, q, k, f, M, d, S, c, rho, g)
    ip = @(x) x(k) + f*(x(min(k+1,numel(x))) - x(k));
    V = ip(d.V);  th = ip(d.th);  ph = ip(d.ph);  ue = ip(d.ue);
    qbar = 0.5*rho*V^2;  qh = q*c/hatf(M,V);
    CL = M.CL0 + M.CLa*al + M.CLq*qh + M.CLde*M.uscale*ue;
    Cm = M.Cm0 + M.Cma*al + M.Cmq*qh + M.Cmde*M.uscale*ue;
    ad = q + (g/V)*cos(th)*cos(ph) - (qbar*S/(M.m*V))*CL;
    qd = (qbar*S*c/M.Jy)*Cm;
end
function y = sim_lat(M, d, S, b, rho, g)
    N = numel(d.t);  nsub = 5;  h = d.dt/nsub;
    be = zeros(N,1); p = zeros(N,1); r = zeros(N,1);
    p(1) = d.p(1);  r(1) = d.r(1);
    qb1 = 0.5*rho*d.V(1)^2;
    be(1) = max(min(d.ay(1)*M.m/(qb1*S*M.CYb), 0.3), -0.3);
    for k = 1:N-1
        bb = be(k); pp = p(k); rr = r(k);
        for s = 1:nsub
            f = (s-1)/nsub;  fh = f+0.5/nsub;  f1 = f+1/nsub;
            [a1,a2,a3] = dt_(bb,pp,rr, k,f,  M,d,S,b,rho,g);
            [b1,b2,b3] = dt_(bb+h/2*a1, pp+h/2*a2, rr+h/2*a3, k,fh, M,d,S,b,rho,g);
            [c1,c2,c3] = dt_(bb+h/2*b1, pp+h/2*b2, rr+h/2*b3, k,fh, M,d,S,b,rho,g);
            [e1,e2,e3] = dt_(bb+h*c1,   pp+h*c2,   rr+h*c3,   k,f1, M,d,S,b,rho,g);
            bb = bb + h/6*(a1+2*b1+2*c1+e1);
            pp = pp + h/6*(a2+2*b2+2*c2+e2);
            rr = rr + h/6*(a3+2*b3+2*c3+e3);
        end
        be(k+1) = max(min(bb,0.5),-0.5);  p(k+1) = max(min(pp,20),-20);  r(k+1) = max(min(rr,20),-20);
    end
    qbar = 0.5*rho*d.V.^2;  hh = hatf(M, d.V);
    y.p = p;  y.r = r;  y.be = be;
    y.ay = (qbar*S/M.m).*(M.CYb*be + M.CYp*(p*b./hh) + M.CYr*(r*b./hh) + M.CYda*M.uscale*d.ua);
end
function [bd, pd, rd] = dt_(be, p, r, k, f, M, d, S, b, rho, g)
    ip = @(x) x(k) + f*(x(min(k+1,numel(x))) - x(k));
    V = ip(d.V);  th = ip(d.th);  ph = ip(d.ph);  q = ip(d.q);  ua = ip(d.ua)*M.uscale;
    qbar = 0.5*rho*V^2;  hh = hatf(M,V);  ph_ = p*b/hh;  rh_ = r*b/hh;
    CY = M.CYb*be + M.CYp*ph_ + M.CYr*rh_ + M.CYda*ua;
    Cl = M.Clb*be + M.Clp*ph_ + M.Clr*rh_ + M.Clda*ua;
    Cn = M.Cnb*be + M.Cnp*ph_ + M.Cnr*rh_ + M.Cnda*ua;
    L = qbar*S*b*Cl;  Nn = qbar*S*b*Cn;
    bd = (qbar*S/(M.m*V))*CY - r + (g/V)*cos(th)*sin(ph);
    pd = M.g1*p*q - M.g2*q*r + M.g3*L + M.g4*Nn;
    rd = M.g7*p*q - M.g1*q*r + M.g4*L + M.g8*Nn;
end
function tab(titulo, canais, W, M, R)
    fprintf('\n  ---- %s: R² por janela\n', titulo);
    for ch = 1:numel(canais)
        fprintf('\n   %s:  %-12s', canais{ch}, 'janela');
        for k = 1:numel(M), fprintf(' %22s', M{k}.nome); end
        fprintf('\n');
        for w = 1:size(W,1)
            fprintf('        %4.0f–%-7.0f', W(w,1), W(w,2));
            for k = 1:numel(M), fprintf(' %22.3f', R(w,k,ch)); end
            fprintf('\n');
        end
        fprintf('        %-12s', 'MÉDIA');
        for k = 1:numel(M), fprintf(' %22.3f', mean(R(:,k,ch),'omitnan')); end
        fprintf('\n');
    end
end
function [C, LW] = paleta(n)
    ref = [0.55 0.68 0.82; 0.72 0.82 0.66; 0.82 0.70 0.86];
    C = [ref(1:n-1,:); 0.85 0.37 0.01];  LW = [repmat(1.2,n-1,1); 2.0];
end
function sub = subtit(own)
    if own, sub = 'cada modelo com a SUA estrutura e o SEU controle';
    else,   sub = 'cada modelo com a SUA estrutura (m, J, normalização), controle padronizado'; end
end
function fig_lon(d, Y, M, R, tw, DIR, Vm, own, tag)
    if nargin < 9, tag = ''; end
    opt.tag = tag;
    n = numel(M);  [C, LW] = paleta(n);
    f = figure('Position',[30 30 1250 830],'Color','w','Visible','off');
    try, f.Theme='light'; catch, end
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    nexttile; hold on; grid on;
    for k = 1:n, plot(d.t, Y{k}.q, '-','Color',C(k,:),'LineWidth',LW(k)); end
    plot(d.t, d.q, 'k-','LineWidth',2.2);
    ylabel('q [rad/s]'); xlim(tw); ylim(lim_(d.q));
    lg = cell(n+1,1);
    for k = 1:n, lg{k} = sprintf('%s  (R² %.2f)', M{k}.nome, R(1,k,1)); end
    lg{n+1} = 'medido';  legend(lg,'Location','eastoutside','FontSize',8.5);
    nexttile; hold on; grid on;
    for k = 1:n, plot(d.t, Y{k}.az, '-','Color',C(k,:),'LineWidth',LW(k)); end
    plot(d.t, d.az, 'k-','LineWidth',2.2);
    ylabel('a_z [m/s²]'); xlim(tw); ylim(lim_(d.az));
    lg2 = cell(n+1,1);
    for k = 1:n, lg2{k} = sprintf('R² %.2f', R(1,k,2)); end
    lg2{n+1} = 'medido';  legend(lg2,'Location','eastoutside','FontSize',8.5);
    nexttile; hold on; grid on;
    plot(d.t, d.ue, '-','Color',[0.3 0.3 0.3],'LineWidth',1.4);
    ylabel('\delta_e normalizado'); xlabel('t [s]'); xlim(tw); ylim([-1 1]);
    title(tl, sprintf('Asa fixa LONGITUDINAL, validação %.0f a %.0f s (V %.1f m/s)\n%s', tw, Vm, subtit(own)), 'FontWeight','bold');
    exportgraphics(f, fullfile(DIR, sprintf(['exec_long%s_%04.0f-%04.0f.png'], opt.tag, tw)), 'BackgroundColor','white','Resolution',130);
    close(f);
end
function fig_lat(d, Y, M, R, tw, DIR, Vm, own, tag)
    if nargin < 9, tag = ''; end
    opt.tag = tag;
    n = numel(M);  [C, LW] = paleta(n);
    f = figure('Position',[30 30 1250 980],'Color','w','Visible','off');
    try, f.Theme='light'; catch, end
    tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
    ch = {'p','r','ay'};  rot = {'p [rad/s]','r [rad/s]','a_y [m/s²]'};
    for j = 1:3
        nexttile; hold on; grid on;
        for k = 1:n, plot(d.t, Y{k}.(ch{j}), '-','Color',C(k,:),'LineWidth',LW(k)); end
        plot(d.t, d.(ch{j}), 'k-','LineWidth',2.2);
        ylabel(rot{j}); xlim(tw); ylim(lim_(d.(ch{j})));
        lg = cell(n+1,1);
        for k = 1:n
            if j == 1, lg{k} = sprintf('%s  (R² %.2f)', M{k}.nome, R(1,k,j));
            else,      lg{k} = sprintf('R² %.2f', R(1,k,j)); end
        end
        lg{n+1} = 'medido';  legend(lg,'Location','eastoutside','FontSize',8.5);
    end
    nexttile; hold on; grid on;
    plot(d.t, d.ua, '-','Color',[0.3 0.3 0.3],'LineWidth',1.4);
    ylabel('\delta_a normalizado'); xlabel('t [s]'); xlim(tw); ylim([-1 1]);
    title(tl, sprintf('Asa fixa LÁTERO-DIRECIONAL, validação %.0f a %.0f s (V %.1f m/s)\n%s', tw, Vm, subtit(own)), 'FontWeight','bold');
    exportgraphics(f, fullfile(DIR, sprintf(['exec_lat%s_%04.0f-%04.0f.png'],  opt.tag, tw)), 'BackgroundColor','white','Resolution',130);
    close(f);
end
function L = lim_(y)
%LIM_  eixo pelo dado medido, com folga: modelo que dispara não estica a escala
    a = min(y); z = max(y); d = max(z-a, eps);
    L = [a - 0.9*d, z + 0.9*d];
end
