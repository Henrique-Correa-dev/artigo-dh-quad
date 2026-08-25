% fw_lateral.m — dinâmica látero-direcional do DH em asa fixa, dez/2024
% =========================================================================
% Identifica o modelo látero-direcional nos doublets de aileron do voo de
% 17/12/2024 (MANUAL, malha aberta).
%
% MODELO (β, p, r integrados; V, θ, φ e q vêm medidos, o que evita a deriva da
% atitude integrada e mantém o teste focado nas derivadas aerodinâmicas)
%   β̇ = (q̄S/(mV))·C_Y − r + (g/V)·cosθ·senφ
%   ṗ = Γ1·p·q − Γ2·q·r + Γ3·L + Γ4·N
%   ṙ = Γ7·p·q − Γ1·q·r + Γ4·L + Γ8·N
%   a_y = (q̄S/m)·C_Y
% com L = q̄Sb·C_l , N = q̄Sb·C_n e os Γ do tensor de inércia (Beard-McLain),
%   C_Y = C_Yβ·β + C_Yp·p̂ + C_Yr·r̂ + C_Yδa·u_a
%   C_l = C_lβ·β + C_lp·p̂ + C_lr·r̂ + C_lδa·u_a
%   C_n = C_nβ·β + C_np·p̂ + C_nr·r̂ + C_nδa·u_a
% e p̂ = pb/2V, r̂ = rb/2V (mesma convenção da AVL).
%
% β NÃO é medido (a mensagem AOA do log é toda zero e o EKF não estimou vento),
% então entra como estado, observado por a_y. Saídas: p, r e a_y.
%
% O leme fica no zero nesses doublets, então C_lδr, C_nδr e C_Yδr não são
% identificáveis aqui e ficam fora. As derivadas cruzadas pequenas (C_Yp, C_Yr,
% C_lr, C_np) ficam FIXAS nos valores da AVL.
%
% Deflexão de aileron NORMALIZADA em [-1,1]: o identificado é C_lδa·δ_max.
% ESCALAS ASSUMIDAS: m = 1,993 kg e o tensor de inércia do CAD.
%
% Uso:  >> fw_lateral
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();

WINS = [303.5 307.5; 312.5 316.5; 322.5 326.5; 333.0 337.0];
% Levantamento de eventos (scratchpad/fw_events.m) mostrou que os 4 doublets
% acima são só parte da manobra disponível: sobraram 77 eventos de aileron em
% asa fixa, entre eles 3 doublets de amplitude cheia. Hook para ampliar o
% conjunto sem editar o arquivo:  setappdata(0,'fw_lat_wins', W)
W_ = getappdata(0,'fw_lat_wins');  if ~isempty(W_), WINS = W_; end
S = proj.wing.S;  b = proj.wing.b;  rho = proj.rho;  g = proj.g;  m = proj.m;
Jx = proj.J.Jx; Jy = proj.J.Jy; Jz = proj.J.Jz; Jxz = proj.J.Jxz;
gam = Jx*Jz - Jxz^2;
G = struct('g1', Jxz*(Jx-Jy+Jz)/gam, 'g2', (Jz*(Jz-Jy)+Jxz^2)/gam, 'g3', Jz/gam, ...
           'g4', Jxz/gam, 'g7', (Jx*(Jx-Jy)+Jxz^2)/gam, 'g8', Jx/gam);
FIX = struct('CYp',0.1121, 'CYr',0.1540, 'Clr',0.1420, 'Cnp',-0.0319);   % AVL, fixas

DMAX = 25;                       % [°] curso mecânico suposto do aileron
AVL = struct('CYb',-0.1934, 'CYda',0.000273*DMAX, ...
             'Clb',-0.0615, 'Clp',-0.4060, 'Clda',-0.004085*DMAX, ...
             'Cnb', 0.0656, 'Cnr',-0.0696, 'Cnda',0.000073*DMAX);
pn = {'CYb','CYda','Clb','Clp','Clda','Cnb','Cnr','Cnda'};
P0 = cellfun(@(f) AVL.(f), pn)';
%              CYb    CYda   Clb    Clp    Clda   Cnb   Cnr    Cnda
lb = [-2.0;  -1.0;  -1.0;  -3.0;  -2.0;  -0.5;  -2.0;  -1.0];
ub = [ 0.0;   1.0;   0.0;   0.0;   2.0;   1.0;   0.0;   1.0];
% Sinais: amortecimentos negativos, C_lβ < 0 (diedro), C_nβ > 0 (cata-vento).
% Os de controle ficam com sinal LIVRE, pela mesma inversão de servo vista no
% profundor (deflexão positiva no log dá momento de nariz para cima).

F = fw_load('lat');
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

fprintf('\n  DH asa fixa, dez/2024 — látero-direcional em %d doublets de aileron\n', size(WINS,1));
fprintf('  m = %.3f kg (suposta), J_x = %.4f, J_z = %.4f kg·m² (CAD), S = %.3f m², b = %.2f m\n\n', ...
    m, Jx, Jz, S, b);

Pid = nan(numel(P0), size(WINS,1));  R2w = nan(size(WINS,1),3);  Vw = nan(size(WINS,1),1);
D = cell(size(WINS,1),1);
for w = 1:size(WINS,1)
    ii = F.t >= WINS(w,1) & F.t <= WINS(w,2);
    d = struct('t',F.t(ii), 'dt',F.dt, 'V',max(F.V(ii),5), 'th',F.theta(ii), 'ph',F.phi(ii), ...
               'q',F.q(ii), 'ua',F.da(ii), 'p',F.p(ii), 'r',F.r(ii), 'ay',F.ay(ii));
    wts = [1/var(d.p); 1/var(d.r); 1/var(d.ay)];
    cost = @(P) resid(P, d, wts, S, b, rho, g, m, G, FIX);
    opts = optimoptions('lsqnonlin','Display','off','MaxIterations',400, ...
        'MaxFunctionEvaluations',20000,'FunctionTolerance',1e-12,'StepTolerance',1e-12);
    [P, ~, res, ~, ~, ~, J] = lsqnonlin(cost, P0, lb, ub, opts);
    y = sim_lat(P, d, S, b, rho, g, m, G, FIX);
    Pid(:,w) = P;  Vw(w) = mean(d.V);
    R2w(w,:) = [R2(d.p,y.p), R2(d.r,y.r), R2(d.ay,y.ay)];
    sig2 = sum(res.^2)/max(numel(res)-numel(P),1);
    d.se = sqrt(max(diag(sig2*pinv(full(J.'*J))),0));  d.y = y;  d.P = P;  D{w} = d;
    fprintf('  janela %d (%.0f a %.0f s, V %.1f m/s):  R² p = %.3f | R² r = %.3f | R² a_y = %.3f\n', ...
        w, WINS(w,1), WINS(w,2), mean(d.V), R2w(w,:));
end

nW = size(WINS,1);
fprintf('\n  %-6s %9s |', 'coef','AVL');
for w = 1:nW, fprintf(' %8s', sprintf('jan.%d',w)); end
fprintf(' | %9s %9s %7s\n', 'média','desvio','disp%');
for i = 1:numel(pn)
    fprintf('  %-6s %9.4f |', pn{i}, P0(i));
    fprintf(' %8.4f', Pid(i,:));
    mu = mean(Pid(i,:));  sd = std(Pid(i,:));
    fprintf(' | %9.4f %9.4f %6.0f%%\n', mu, sd, 100*sd/max(abs(mu),1e-9));
end
fprintf('\n  CR%% por janela:\n');
for i = 1:numel(pn)
    cr = arrayfun(@(w) 100*D{w}.se(i)/max(abs(D{w}.P(i)),1e-9), 1:numel(D));
    fprintf('  %-6s %s\n', pn{i}, sprintf('%8.1f', cr));
end
Vm = mean(Vw);  qbar = 0.5*rho*Vm^2;
Lp = qbar*S*b/Jx * mean(Pid(4,:)) * b/(2*Vm);
Nr = qbar*S*b/Jz * mean(Pid(7,:)) * b/(2*Vm);
fprintf('\n  Em V = %.1f m/s:  L_p = %.2f 1/s (constante de rolagem %.2f s)   N_r = %.3f 1/s\n', ...
    Vm, Lp, -1/Lp, Nr);
save(fullfile(paths.outputs,'fw_lateral.mat'), 'Pid','R2w','WINS','pn','P0','Vw','DMAX');

%% figura
f = figure('Position',[30 30 1350 900],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(4, numel(D), 'TileSpacing','compact','Padding','compact');
for w = 1:numel(D)
    d = D{w};
    nexttile(w); hold on; grid on;
    plot(d.t, d.ua, '-','Color',[0 0.45 0.7],'LineWidth',1.4);
    ylabel('\delta_a normalizado'); xlim(WINS(w,:)); ylim([-1.1 1.1]);
    text(0.02, 1.09, sprintf('doublet %d  (V %.1f m/s)', w, mean(d.V)), 'Units','normalized','FontWeight','bold');
    ch = {'p','r','ay'};  lb_ = {'p [rad/s]','r [rad/s]','a_y [m/s²]'};
    for j = 1:3
        nexttile(j*numel(D)+w); hold on; grid on;
        plot(d.t, d.(ch{j}), 'k-','LineWidth',1.6);
        plot(d.t, d.y.(ch{j}), '-','Color',[0.85 0.37 0.01],'LineWidth',1.3);
        ylabel(lb_{j}); xlim(WINS(w,:));
        text(0.02, 0.06, sprintf('R² = %.3f', R2w(w,j)), 'Units','normalized','FontSize',9,'BackgroundColor',[1 1 1 0.7]);
        if j == 1 && w == 1, legend({'medido','modelo identificado'},'Location','northwest','FontSize',8); end
        if j == 3, xlabel('t [s]'); end
    end
end
title(tl, 'DH asa fixa dez/2024 — látero-direcional identificado em cada doublet de aileron');
fn = fullfile(paths.images,'new_flights','fw_lateral.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',130);
fprintf('  Figura: %s\n', fn);

%% ------------------------------------------------------------- funções
function e = resid(P, d, wts, S, b, rho, g, m, G, FIX)
    y = sim_lat(P, d, S, b, rho, g, m, G, FIX);
    e = [sqrt(wts(1))*(d.p - y.p); sqrt(wts(2))*(d.r - y.r); sqrt(wts(3))*(d.ay - y.ay)];
    e(~isfinite(e)) = 1e3;
end

function y = sim_lat(P, d, S, b, rho, g, m, G, FIX)
    N = numel(d.t);  nsub = 5;  h = d.dt/nsub;
    be = zeros(N,1); p = zeros(N,1); r = zeros(N,1);
    p(1) = d.p(1);  r(1) = d.r(1);
    qbar1 = 0.5*rho*d.V(1)^2;
    be(1) = max(min(d.ay(1)*m/(qbar1*S*P(1)), 0.3), -0.3);     % β inicial pelo a_y medido
    for k = 1:N-1
        bb = be(k); pp = p(k); rr = r(k);
        for s = 1:nsub
            f = (s-1)/nsub;  fh = f + 0.5/nsub;  f1 = f + 1/nsub;
            [d1,d2,d3] = dyn(bb,pp,rr, k,f,  P,d,S,b,rho,g,m,G,FIX);
            [e1,e2,e3] = dyn(bb+h/2*d1, pp+h/2*d2, rr+h/2*d3, k,fh, P,d,S,b,rho,g,m,G,FIX);
            [f1_,f2,f3] = dyn(bb+h/2*e1, pp+h/2*e2, rr+h/2*e3, k,fh, P,d,S,b,rho,g,m,G,FIX);
            [g1_,g2,g3] = dyn(bb+h*f1_, pp+h*f2, rr+h*f3, k,f1, P,d,S,b,rho,g,m,G,FIX);
            bb = bb + h/6*(d1+2*e1+2*f1_+g1_);
            pp = pp + h/6*(d2+2*e2+2*f2+g2);
            rr = rr + h/6*(d3+2*e3+2*f3+g3);
        end
        be(k+1) = max(min(bb,0.5),-0.5);  p(k+1) = max(min(pp,20),-20);  r(k+1) = max(min(rr,20),-20);
    end
    qbar = 0.5*rho*d.V.^2;
    CY = P(1)*be + FIX.CYp*(p*b./(2*d.V)) + FIX.CYr*(r*b./(2*d.V)) + P(2)*d.ua;
    y.p = p;  y.r = r;  y.be = be;  y.ay = (qbar*S/m).*CY;
end

function [bd, pd, rd] = dyn(be, p, r, k, f, P, d, S, b, rho, g, m, G, FIX)
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
