% fw_longitudinal.m — período curto do DH em asa fixa, dez/2024
% =========================================================================
% Identifica o modelo de PERÍODO CURTO nos doublets de profundor do voo de
% 17/12/2024, feito inteiro em MANUAL (sem estabilização, malha aberta pura).
%
% MODELO (eixos do corpo, V medida pelo pitot como entrada exógena)
%   ᾱ = q + (g/V)·cosθ·cosφ − (q̄S/(mV))·C_L
%   q̇ = (q̄S c̄/J_y)·C_m
%   a_z = −(q̄S/m)·C_L                        (força específica, sem gravidade)
% com q̄ = ½ρV² e
%   C_L = C_L0 + C_Lα·α + C_Lq·(q c̄/2V) + C_Lδe·u_e
%   C_m = C_m0 + C_mα·α + C_mq·(q c̄/2V) + C_mδe·u_e
%
% O ângulo de ataque NÃO é medido (a mensagem AOA do log é toda zero e o EKF
% não estimou vento), então α é ESTADO, observado por a_z. É a formulação
% padrão de ensaio em voo: saídas q e a_z, dois estados.
%
% u_e é a deflexão NORMALIZADA em [-1,1] lida do RCOU, porque o curso mecânico
% da superfície não está no log. Logo o que se identifica é o produto
% C_mδe·δ_max. Para comparar com a AVL, que dá C_mδe por grau, use
% C_mδe(AVL)·δ_max: com δ_max = 25° dá −0,562.
%
% ESCALAS ASSUMIDAS: m = 1,993 kg e J_y = 0,0844 kg·m² (CAD). A massa da
% configuração de asa fixa em dez/2024 não foi pesada. C_Lα escala com m e os
% C_m escalam com J_y, linearmente. O que o dado determina sem essa hipótese
% são as derivadas DIMENSIONAIS M_α [1/s²] e M_q [1/s], reportadas junto.
%
% Uso:  >> fw_longitudinal
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();

WINS = [306.5 310.5; 326.5 330.5; 344.0 348.0; 367.5 371.5];   % doublets de ±1
% Hook para reduzir o conjunto de identificação e liberar doublet para validação:
%   setappdata(0,'fw_lon_wins', W)
W_ = getappdata(0,'fw_lon_wins');  if ~isempty(W_), WINS = W_; end
TAG_ = getappdata(0,'fw_lon_tag'); if isempty(TAG_), TAG_ = ''; end
S = proj.wing.S;  c = proj.wing.c;  rho = proj.rho;  g = proj.g;
m = proj.m;  Jy = proj.J.Jy;
CLq_fix = 8.112;                 % AVL, mantido fixo (efeito pequeno e colinear com C_mq)

% Θ₀ da AVL (dh_st.txt), com o profundor convertido para deflexão normalizada
DMAX = 25;                       % [°] curso mecânico suposto do profundor
AVL = struct('CL0',0.35, 'CLa',4.585, 'CLde',0.010387*DMAX, ...
             'Cm0',0.0,  'Cma',-0.7377,'Cmq',-8.959, 'Cmde',-0.022483*DMAX);
pn = {'CL0','CLa','CLde','Cm0','Cma','Cmq','Cmde'};
P0 = cellfun(@(f) AVL.(f), pn)';
% SINAL DAS DERIVADAS DE CONTROLE: livre. No log, deflexão positiva do
% profundor produz q POSITIVO (nariz para cima), como se vê no doublet: δ_e vai
% a +1 e q sobe a +1,8 rad/s. Ou seja, a convenção do servo é invertida em
% relação à clássica (positivo = bordo de fuga para baixo = nariz para baixo),
% e C_mδe sai POSITIVO em relação ao sinal registrado. É a mesma inversão que
% aparece no relatório do FINEP, onde C_mδe foi reportado positivo.
lb = [-1.0;  1.0; -4.0;  -0.6; -6.0; -80.0; -6.0];
ub = [ 1.5; 20.0;  4.0;   0.6;  0.0;   0.0;  6.0];

F = fw_load('long');
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
TIC = @(y,yh) sqrt(mean((y-yh).^2)) / (sqrt(mean(y.^2)) + sqrt(mean(yh.^2)));

fprintf('\n  DH asa fixa, dez/2024 — período curto em %d doublets de profundor\n', size(WINS,1));
fprintf('  m = %.3f kg (suposta), J_y = %.4f kg·m² (CAD), S = %.3f m², c̄ = %.3f m\n\n', m, Jy, S, c);

Pid = nan(numel(P0), size(WINS,1));  R2w = nan(size(WINS,1), 2);  Vw = nan(size(WINS,1),1);
D = cell(size(WINS,1),1);
for w = 1:size(WINS,1)
    ii = F.t >= WINS(w,1) & F.t <= WINS(w,2);
    d = struct('t',F.t(ii), 'dt',F.dt, 'V',F.V(ii), 'th',F.theta(ii), 'ph',F.phi(ii), ...
               'ue',F.de(ii), 'q',F.q(ii), 'az',F.az(ii));
    d.V = max(d.V, 5);                       % guarda contra leitura ruim do pitot
    wq = 1/var(d.q);  wa = 1/var(d.az);
    cost = @(P) resid(P, d, wq, wa, S, c, rho, g, m, Jy, CLq_fix);
    opts = optimoptions('lsqnonlin','Display','off','MaxIterations',400, ...
        'MaxFunctionEvaluations',20000,'FunctionTolerance',1e-12,'StepTolerance',1e-12);
    [P, ~, res, ~, ~, ~, J] = lsqnonlin(cost, P0, lb, ub, opts);
    y = sim_sp(P, d, S, c, rho, g, m, Jy, CLq_fix);
    Pid(:,w) = P;  Vw(w) = mean(d.V);
    R2w(w,:) = [R2(d.q, y.q), R2(d.az, y.az)];
    % Cramér-Rao
    sig2 = sum(res.^2)/max(numel(res)-numel(P),1);
    se = sqrt(max(diag(sig2*pinv(full(J.'*J))),0));
    d.y = y;  d.P = P;  d.se = se;  D{w} = d;
    fprintf('  janela %d (%.0f a %.0f s, V %.1f m/s):  R² q = %.3f | R² a_z = %.3f\n', ...
        w, WINS(w,1), WINS(w,2), mean(d.V), R2w(w,:));
end

nW = size(WINS,1);
fprintf('\n  %-6s %9s |', 'coef', 'AVL');
for w = 1:nW, fprintf(' %9s', sprintf('jan.%d',w)); end
fprintf(' | %9s %9s\n', 'média', 'desvio');
for i = 1:numel(pn)
    fprintf('  %-6s %9.3f |', pn{i}, P0(i));
    fprintf(' %9.3f', Pid(i,:));
    fprintf(' | %9.3f %9.3f\n', mean(Pid(i,:)), std(Pid(i,:)));
end
fprintf('\n  CR%% por janela (erro-padrão relativo):\n');
for i = 1:numel(pn)
    cr = arrayfun(@(w) 100*D{w}.se(i)/max(abs(D{w}.P(i)),1e-9), 1:numel(D));
    fprintf('  %-6s %s\n', pn{i}, sprintf('%8.1f', cr));
end

% derivadas dimensionais (independentes da hipótese de J_y a menos da escala)
Vm = mean(Vw);  qbar = 0.5*rho*Vm^2;
Ma = qbar*S*c/Jy * mean(Pid(5,:));
Mq = qbar*S*c/Jy * mean(Pid(6,:)) * c/(2*Vm);
wn = sqrt(max(-Ma,0));   % aproximação de período curto ignorando Z_α
fprintf('\n  Em V = %.1f m/s:  M_α = %.1f 1/s²   M_q = %.2f 1/s   ω_n(aprox) = %.2f rad/s (%.2f Hz)\n', ...
    Vm, Ma, Mq, wn, wn/(2*pi));

save(fullfile(paths.outputs,['fw_longitudinal' TAG_ '.mat']), 'Pid','R2w','WINS','pn','P0','Vw','DMAX');

%% figura
f = figure('Position',[30 30 1350 760],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(3, numel(D), 'TileSpacing','compact','Padding','compact');
for w = 1:numel(D)
    d = D{w};
    nexttile(w); hold on; grid on;
    plot(d.t, d.ue, '-','Color',[0.85 0.37 0.01],'LineWidth',1.4);
    ylabel('\delta_e normalizado'); xlim(WINS(w,:)); ylim([-1.1 1.1]);
    text(0.02, 1.09, sprintf('doublet %d  (V %.1f m/s)', w, mean(d.V)), 'Units','normalized','FontWeight','bold');
    nexttile(numel(D)+w); hold on; grid on;
    plot(d.t, d.q, 'k-','LineWidth',1.6); plot(d.t, d.y.q, '-','Color',[0 0.45 0.7],'LineWidth',1.3);
    ylabel('q [rad/s]'); xlim(WINS(w,:));
    text(0.02, 0.06, sprintf('R² = %.3f', R2w(w,1)), 'Units','normalized','FontSize',9,'BackgroundColor',[1 1 1 0.7]);
    if w == 1, legend({'medido','modelo identificado'},'Location','northwest','FontSize',8); end
    nexttile(2*numel(D)+w); hold on; grid on;
    plot(d.t, d.az, 'k-','LineWidth',1.6); plot(d.t, d.y.az, '-','Color',[0 0.45 0.7],'LineWidth',1.3);
    ylabel('a_z [m/s²]'); xlabel('t [s]'); xlim(WINS(w,:));
    text(0.02, 0.06, sprintf('R² = %.3f', R2w(w,2)), 'Units','normalized','FontSize',9,'BackgroundColor',[1 1 1 0.7]);
end
title(tl, 'DH asa fixa dez/2024 — período curto identificado em cada doublet de profundor');
fn = fullfile(paths.images,'new_flights','fw_longitudinal.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',130);
fprintf('  Figura: %s\n', fn);

%% ------------------------------------------------------------- funções
function e = resid(P, d, wq, wa, S, c, rho, g, m, Jy, CLq)
    y = sim_sp(P, d, S, c, rho, g, m, Jy, CLq);
    e = [sqrt(wq)*(d.q - y.q); sqrt(wa)*(d.az - y.az)];
    e(~isfinite(e)) = 1e3;
end

function y = sim_sp(P, d, S, c, rho, g, m, Jy, CLq)
% Integra o período curto com RK4 e sub-passos, V, θ e φ vindos do dado.
    CL0=P(1); CLa=P(2); CLde=P(3); Cm0=P(4); Cma=P(5); Cmq=P(6); Cmde=P(7);
    N = numel(d.t);  nsub = 5;  h = d.dt/nsub;
    al = zeros(N,1); q = zeros(N,1);
    q(1) = d.q(1);
    al(1) = trim_alpha(P, d, 1, S, rho, g, m);      % α inicial pelo a_z medido
    for k = 1:N-1
        a = al(k); qq = q(k);
        for s = 1:nsub
            f = (s-1)/nsub;  fh = f + 0.5/nsub;  f1 = f + 1/nsub;
            [k1a,k1q] = dyn(a,           qq,           k, f,  P, d, S, c, rho, g, m, Jy, CLq);
            [k2a,k2q] = dyn(a+h/2*k1a,   qq+h/2*k1q,   k, fh, P, d, S, c, rho, g, m, Jy, CLq);
            [k3a,k3q] = dyn(a+h/2*k2a,   qq+h/2*k2q,   k, fh, P, d, S, c, rho, g, m, Jy, CLq);
            [k4a,k4q] = dyn(a+h*k3a,     qq+h*k3q,     k, f1, P, d, S, c, rho, g, m, Jy, CLq);
            a  = a  + h/6*(k1a + 2*k2a + 2*k3a + k4a);
            qq = qq + h/6*(k1q + 2*k2q + 2*k3q + k4q);
        end
        al(k+1) = max(min(a, 0.6), -0.6);  q(k+1) = max(min(qq, 20), -20);
    end
    qbar = 0.5*rho*d.V.^2;
    CL = CL0 + CLa*al + CLq*(q.*c./(2*d.V)) + CLde*d.ue;
    y.q = q;  y.al = al;  y.az = -(qbar*S/m).*CL;
end

function [ad, qd] = dyn(al, q, k, f, P, d, S, c, rho, g, m, Jy, CLq)
    ip = @(x) x(k) + f*(x(min(k+1,numel(x))) - x(k));
    V = ip(d.V);  th = ip(d.th);  ph = ip(d.ph);  ue = ip(d.ue);
    qbar = 0.5*rho*V^2;
    qhat = q*c/(2*V);
    CL = P(1) + P(2)*al + CLq*qhat + P(3)*ue;
    Cm = P(4) + P(5)*al + P(6)*qhat + P(7)*ue;
    ad = q + (g/V)*cos(th)*cos(ph) - (qbar*S/(m*V))*CL;
    qd = (qbar*S*c/Jy)*Cm;
end

function a0 = trim_alpha(P, d, k, S, rho, g, m)
% α inicial coerente com o a_z medido: −a_z·m/(q̄S) = C_L → α = (C_L − C_L0 − C_Lδe·u)/C_Lα
    qbar = 0.5*rho*d.V(k)^2;
    CL = -d.az(k)*m/(qbar*S);
    a0 = (CL - P(1) - P(3)*d.ue(k)) / max(P(2), 0.5);
    a0 = max(min(a0, 0.4), -0.4);
end
