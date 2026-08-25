% prior_damping.m — Valores A PRIORI dos arrastos angulares cp, cq, cr
% =========================================================================
% IDEIA
%   Os coeficientes cp, cq, cr do modelo rotacional (ṗ = ... − cp·p, etc.) são
%   derivadas de amortecimento do voo pairado. A fonte física dominante, que
%   se consegue derivar em forma fechada, é a variação do empuxo de cada rotor
%   com a velocidade axial induzida pela rotação do corpo (hipótese de influxo):
%
%     rotor i sobe/desce a  V_c,i = q·x_i − p·y_i   →   ΔT_i = −k_v·V_c,i
%     ΔL = −4·k_v·ly²·p ,  ΔM = −4·k_v·lx²·q         (configuração simétrica)
%     cp_0 = 4·k_v·ly²/Jx ,  cq_0 = 4·k_v·lx²/Jy
%
%   Guinada: rotor anda de lado a r·d, sofre força H (arrasto de rotor no plano):
%     ΔN = −4·k_h·d²·r  →  cr_0 = 4·k_h·d²/Jz
%
%   k_v pela teoria de elemento de pá + quantidade de movimento (derivada Z_w de
%   Padfield 2007, Bramwell 2001, Leishman 2006), reescrita como fração da
%   escala natural T0/vh:
%     k_v = (T0/vh)·κ ,   κ = a·σ/(16·λ_i + a·σ)
%   k_h por elemento de pá (arrasto de perfil + parcela induzida), Leishman 2006:
%     k_h = ρ·A·(ΩR)·(σ/4)·(Cd0 + a·θ75·λ_i)
%
%   Além disso, estima a ordem de grandeza dos termos aerodinâmicos da ESTRUTURA
%   (asa + empenagem, geometria de DeLucena2025) que crescem com a velocidade
%   de avanço V, para delimitar o envelope em que a hipótese "baixa velocidade"
%   vale:  amortecimento ∝ V (Cl_p, Cm_q, Cn_r) e termos estáticos ∝ V² (Cm_α, ...).
%
% FONTES DOS DADOS
%   parameters.m  : massa, braços, inércias CAD, asa, bancada (RPM de pairado)
%   DeLucena2025  : hélice slow flyer 1045, empenagem (0,45×0,17 m; 2 verticais
%                   de ~0,12 m), CG a 1/3 da corda, afilamento 0,8, cm do perfil
%   Hipóteses     : corda média da pá 25 mm, a = 5/rad, Cd0 = 0,03 (Re ~ 5e4–1e5)
%
% Uso:  >> prior_damping
% =========================================================================
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % raiz identification/ (setup_paths)
paths = setup_paths(); p = parameters();
rho = 1.225;  g = p.g;  m = p.m;

%% ---------------- ROTOR: hélice 1045 (10 x 4,5 in), 2 pás ----------------
R        = 0.127;                 % raio [m]
A        = pi*R^2;                % área do disco [m²]
Nb       = 2;                     % número de pás
c_blade  = 0.025;                 % corda média da pá [m]      (hipótese)
a        = 5.0;                   % inclinação de sustentação da seção [1/rad] (hipótese)
Cd0      = 0.03;                  % arrasto de perfil da seção          (hipótese)
pitch_in = 4.5;                   % passo geométrico [in]
sigma    = Nb*c_blade/(pi*R);                       % solidez
theta75  = atan(pitch_in/(0.75*2*R/0.0254*pi));     % ângulo de passo a 75 % R [rad]

T0    = m*g/4;                                       % empuxo por rotor no pairado [N]
rpm_h = interp1(p.bench.T_grams*1e-3*g, p.bench.RPM, T0, 'linear');   % rpm de pairado (bancada)
Omega = rpm_h*2*pi/60;  OmegaR = Omega*R;
vh    = sqrt(T0/(2*rho*A));                          % velocidade induzida no pairado [m/s]
lam_i = vh/OmegaR;                                   % razão de influxo
kappa = a*sigma/(16*lam_i + a*sigma);
k_v   = (T0/vh)*kappa;                               % N/(m/s): empuxo perdido por m/s de subida
k_h   = rho*A*OmegaR*(sigma/4)*(Cd0 + a*theta75*lam_i);   % N/(m/s): força H por m/s lateral

%% ---------------- GEOMETRIA E INÉRCIA (a priori = CAD) ----------------
% Em parameters.m: Lx_* = braço LATERAL (y), Ly_* = braço LONGITUDINAL (x).
ly = mean([p.arms.Lx_r, p.arms.Lx_l]);
lx = mean([p.arms.Ly_f, p.arms.Ly_r]);
d2 = lx^2 + ly^2;
Jx = p.J.Jx;  Jy = p.J.Jy;  Jz = p.J.Jz;

Lp0 = 4*k_v*ly^2;   Mq0 = 4*k_v*lx^2;   Nr0 = 4*k_h*d2;    % N·m·s
cp0 = Lp0/Jx;       cq0 = Mq0/Jy;       cr0 = Nr0/Jz;      % 1/s

%% ---------------- ESTRUTURA: asa + empenagem (DeLucena2025) ----------------
S = p.wing.S;  b = p.wing.b;  cbar = p.wing.c;  AR = b^2/S;
% Derivadas do AVL 3.40, modelo avl/dh.avl (asa + empenagem em H, 472 vórtices),
% eixos de estabilidade, ponto trimado CL = 0,525 (α = 2,03°, V = 15 m/s).
% Saída completa em avl/dh_st.txt. Substituem as estimativas por fórmula usadas
% antes; as fórmulas ficam no comentário ao lado, para comparação.
CLa  =  4.585;      % (fórmula de asa: 4,57)
Cma  = -0.7377;     % (fórmula: -0,95 | XFLR5 do relatório FINEP: -0,783)
Clp  = -0.4060;     % (fórmula: -0,72)
Cmq  = -8.959;      % (fórmula: -12,3)
Cnr  = -0.0696;     % (fórmula: -0,147)
Clb  = -0.0615;     % (típico assumido antes: -0,08)
Cnb  =  0.0656;     % (típico assumido antes: +0,07)
Cm0  = -0.08;       % USA-35B (Xflr5, DeLucena2025)
CL0  = 0.5*AR/(AR+2);                                % sustentação a α = 0 (perfil cambado)
alpha_max = deg2rad(10);  beta_max = deg2rad(10);

qbar = @(V) 0.5*rho*V.^2;
Lp_w = @(V) rho*V*S*b^2   *abs(Clp)/4;              % N·m·s
Mq_w = @(V) rho*V*S*cbar^2*abs(Cmq)/2;
Nr_w = @(V) rho*V*S*b^2   *abs(Cnr)/4;
M_st = @(V) qbar(V)*S*cbar*(abs(Cm0)+abs(Cma)*alpha_max);   % N·m (estático, arfagem)
L_st = @(V) qbar(V)*S*b   *abs(Clb)*beta_max;
N_st = @(V) qbar(V)*S*b   *abs(Cnb)*beta_max;
Lift = @(V) qbar(V)*S*CL0;                           % N (sustentação a α = 0)

%% ---------------- IDENTIFICADO (Θ_OEM usado no modelo linear) ----------------
lm = load(fullfile(paths.outputs,'linear_model.mat'));  P = lm.P(:);
cp = P(13); cq = P(14); cr = P(15);
Lp_id = cp*P(1);  Mq_id = cq*P(2);  Nr_id = cr*P(3);

%% ---------------- RELATÓRIO ----------------
fprintf('\n=============== ROTOR (hélice 1045, pairado) ===============\n');
fprintf('  T0 = %.2f N/rotor | rpm_h = %.0f | ΩR = %.1f m/s | vh = %.2f m/s | λ_i = %.3f\n', T0, rpm_h, OmegaR, vh, lam_i);
fprintf('  σ = %.3f | θ75 = %.1f° | κ = %.2f\n', sigma, rad2deg(theta75), kappa);
fprintf('  k_v = %.3f N/(m/s)   (escala T0/vh = %.2f)\n', k_v, T0/vh);
fprintf('  k_h = %.4f N/(m/s)\n', k_h);
fprintf('  amortecimento vertical (translação): 4k_v/m = %.2f 1/s\n', 4*k_v/m);

fprintf('\n=============== A PRIORI × IDENTIFICADO ===============\n');
fprintf('  braços: ly = %.3f m, lx = %.3f m, d = %.3f m | J CAD = [%.4f %.4f %.4f]\n', ly, lx, sqrt(d2), Jx, Jy, Jz);
fprintf('  %-6s %10s %12s %10s\n', '', 'a priori', 'identificado', 'razão');
fprintf('  %-6s %10.2f %12.2f %10.1f   [1/s]\n', 'cp', cp0, cp, cp/cp0);
fprintf('  %-6s %10.2f %12.2f %10.1f   [1/s]\n', 'cq', cq0, cq, cq/cq0);
fprintf('  %-6s %10.3f %12.3f %10.1f   [1/s]\n', 'cr', cr0, cr, cr/cr0);
fprintf('  razão cp/cq: a priori %.2f | identificado %.2f\n', cp0/cq0, cp/cq);
fprintf('  k_v implicado pelo identificado: rolagem %.2f, arfagem %.2f N/(m/s)\n', Lp_id/(4*ly^2), Mq_id/(4*lx^2));

fprintf('\n=============== ESTRUTURA (asa + empenagem) ===============\n');
fprintf('  AR = %.2f | CLα = %.2f | Cl_p = %.2f | Cm_q = %.1f | Cn_r = %.3f | Cm_α = %.2f /rad | V_H = %.2f\n', ...
    AR, CLa, Clp, Cmq, Cnr, Cma, VH);
fprintf('  amortecimento por m/s de V: L_p %.3f | M_q %.3f | N_r %.3f  [N·m·s por m/s]\n', Lp_w(1), Mq_w(1), Nr_w(1));
V = [1 1.3 2 4];
fprintf('\n  Amortecimento [N·m·s]      rotor    ');  fprintf('  V=%.1f ', V);  fprintf('   identificado\n');
fprintf('  L_p                     %8.3f  ', Lp0);  fprintf('%7.3f ', Lp_w(V));  fprintf('   %8.3f\n', Lp_id);
fprintf('  M_q                     %8.3f  ', Mq0);  fprintf('%7.3f ', Mq_w(V));  fprintf('   %8.3f\n', Mq_id);
fprintf('  N_r                     %8.3f  ', Nr0);  fprintf('%7.3f ', Nr_w(V));  fprintf('   %8.3f\n', Nr_id);
fprintf('  V em que estrutura = identificado: rolagem %.1f m/s | arfagem %.1f m/s | guinada %.1f m/s\n', ...
    Lp_id/Lp_w(1), Mq_id/Mq_w(1), Nr_id/Nr_w(1));
fprintf('  fração explicada (rotor + estrutura) a V=1,3 m/s: L_p %.0f%% | M_q %.0f%% | N_r %.0f%%\n', ...
    100*(Lp0+Lp_w(1.3))/Lp_id, 100*(Mq0+Mq_w(1.3))/Mq_id, 100*(Nr0+Nr_w(1.3))/Nr_id);

fprintf('\n  Termos estáticos (|α|,|β| ≤ 10°) e sustentação da asa\n');
fprintf('  V [m/s]   q̄ [Pa]   M_est [N·m]  L_est [N·m]  N_est [N·m]  Lift [N] (%% mg)\n');
for v = V
    fprintf('  %5.1f   %7.2f   %9.3f   %9.3f   %9.3f   %6.2f (%4.1f%%)\n', ...
        v, qbar(v), M_st(v), L_st(v), N_st(v), Lift(v), 100*Lift(v)/(m*g));
end
fprintf('  referência: momento de amortecimento a 0,5 rad/s ≈ %.2f (rol) %.2f (arf) N·m | controle ≈ 1 N·m\n', ...
    0.5*Lp_id, 0.5*Mq_id);

%% ---------------- SAÍDA PARA USO EM Θ0 ----------------
prior = struct('k_v',k_v,'k_h',k_h,'cp0',cp0,'cq0',cq0,'cr0',cr0, ...
               'Lp_per_V',Lp_w(1),'Mq_per_V',Mq_w(1),'Nr_per_V',Nr_w(1));
save(fullfile(paths.outputs,'prior_damping.mat'),'prior');
fprintf('\n  salvo: %s\n', fullfile(paths.outputs,'prior_damping.mat'));
