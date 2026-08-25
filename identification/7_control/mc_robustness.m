% mc_robustness.m — Robustez do controlador a INCERTEZA PARAMÉTRICA (Monte Carlo)
% =========================================================================
% IDEIA
%   O controlador foi projetado sobre Θ_OEM (nominal). Mas cada parâmetro tem
%   uma incerteza quantificada pelo limite de Cramér-Rao (CR%, Tabela 5.3).
%   Pergunta: se a aeronave REAL tiver parâmetros dentro dessa incerteza,
%   o controlador (fixo, projetado no nominal) continua funcionando bem?
%
% COMO
%   1. Sorteia N plantas: Θ_k = Θ_OEM .* (1 + CR%/100 .* z_k),  z_k ~ N(0,1)
%      → cada parâmetro varia com desvio-padrão igual ao seu CR% (1σ).
%   2. Para cada Θ_k, roda a MALHA FECHADA NÃO LINEAR (cl_loop) com o MESMO
%      controlador K (não re-sintoniza), no cenário de degraus do projeto
%      (altitude 2 m, velocidade 3 m/s, guinada 30°).
%   3. Mede em cada rodada: sobressinal, tempo de acomodação (2%), erro de
%      regime, rolagem máxima, % de saturação, e se ficou estável.
%   4. Reporta média/desvio/pior caso e a taxa de rodadas instáveis.
%
% O QUE ISSO RESPONDE À BANCA
%   "análises de robustez frente a incertezas paramétricas" — com a incerteza
%   que o PRÓPRIO trabalho mediu (CR%), não uma incerteza arbitrária.
%
% Uso:  >> design_control;  mc_robustness
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---------------- MODELO NOMINAL + CONTROLADOR (igual ao sim_control_scenarios) ----------------
lm   = load(fullfile(paths.outputs,'linear_model.mat'));
G    = load(fullfile(paths.control,'control_gains.mat'));
proj = parameters(); [fT,fQ] = motor_models(); dyn = vtol_dynamics('get_handles');
M0 = struct('A',lm.A,'B',lm.B,'u0',lm.u0,'P',lm.P,'K',G.K,'dyn',dyn, ...
            'bridge',forces_to_pwm(lm,'nonlinear'),'fT',fT,'fQ',fQ, ...
            'constants',struct('m',proj.m,'g',proj.g),'dt',0.01,'PWM_LIM',[1200 2000]);
P_nom = lm.P(:);          % Θ_OEM  [Jx Jy Jz Jxz  kT1..4  kQ1..4  cp cq cr]

%% ---------------- INCERTEZA: CR% do run oficial_2026_final (summary.mat), ORDEM de P ----------------
%  Somente os 15 parâmetros LIVRES da identificação são perturbados; os
%  coeficientes aerodinâmicos (fixos, medidos) permanecem nos valores nominais.
sm = load(fullfile(paths.outputs,'runs','oficial_2026_final','summary.mat')).summary;
np_ = numel(P_nom);
sigma = zeros(np_,1);  sigma(1:15) = sm.CR_pct(1:15)/100;   % 1σ relativo

%% ---------------- CENÁRIO: degraus do projeto (Seção "Simulação no Modelo Linear") ----------------
T_end = 25; t = (0:M0.dt:T_end)'; h0 = 0;
s = struct('name','mc','h0',h0,'t',t, ...
    'sp_h',   h0 + 2*(t>=1), ...                 % altitude: 0 → 2 m em t=1 s
    'sp_u',   3*(t>=8), ...                      % velocidade: 0 → 3 m/s em t=8 s
    'sp_psi', deg2rad(30)*(t>=16));              % guinada: 0 → 30° em t=16 s

%% ---------------- MONTE CARLO ----------------
N = 200;  rng(1);                                % reprodutível
Z = zeros(N,15);   % só os livres entram no sorteio
metr = struct('OS_h',nan(N,1),'OS_u',nan(N,1),'ts_h',nan(N,1),'ts_u',nan(N,1), ...
              'e_h',nan(N,1),'e_u',nan(N,1),'e_psi',nan(N,1),'phi_max',nan(N,1), ...
              'sat_pct',nan(N,1),'stable',false(N,1));

% referência: planta nominal
R0 = cl_loop(true, s, M0);
fprintf('\n  Nominal: OS_h=%.1f%%  OS_u=%.1f%%  ts_h=%.2fs  phi_max=%.2f°  sat=%.1f%%\n', ...
    overshoot(R0.H,2,t,1), overshoot(R0.X(:,7),3,t,8), settle(R0.H(t<16),2,t(t<16),1,0.05), ...
    rad2deg(max(abs(R0.X(:,4)))), 100*mean(R0.Usat));

fprintf('\n  Rodando %d plantas perturbadas (1σ = CR%%)...\n', N);
for k = 1:N
    z = zeros(np_,1);  z(1:15) = randn(15,1);  Z(k,:) = z(1:15)';
    P_k = P_nom .* (1 + sigma .* z);
    % proteção física: inércias e coeficientes positivos (Jxz pode ser ~0/±)
    P_k([1:3 5:15]) = max(P_k([1:3 5:15]), 0.05*abs(P_nom([1:3 5:15])));
    Mk = M0;  Mk.P = P_k;
    % ponte PWM também depende de kT/kQ da planta "real"? NÃO: a ponte usa o
    % modelo que o CONTROLADOR acredita (nominal) — é isso que torna o teste
    % de robustez honesto (controlador aloca com Θ_OEM, planta responde com Θ_k).
    try
        R = cl_loop(true, s, Mk);
        phi = rad2deg(R.X(:,4));
        metr.OS_h(k)    = overshoot(R.H,2,t,1);
        metr.OS_u(k)    = overshoot(R.X(:,7),3,t,8);
        mh = t<16;  metr.ts_h(k) = settle(R.H(mh),2,t(mh),1,0.05);   % até o degrau de guinada, banda 5%
        mu = t<16;  metr.ts_u(k) = settle(R.X(mu,7),3,t(mu),8,0.05);
        metr.e_h(k)     = abs(R.H(end)-2);
        metr.e_u(k)     = abs(R.X(end,7)-3);
        metr.e_psi(k)   = rad2deg(abs(R.X(end,6)-deg2rad(30)));
        metr.phi_max(k) = max(abs(phi));
        metr.sat_pct(k) = 100*mean(R.Usat);
        metr.stable(k)  = all(isfinite(R.X(:))) && max(abs(phi))<45 && ...
                          all(abs(R.H)<20) && metr.e_h(k)<0.5;
    catch
        metr.stable(k) = false;   % divergência numérica = instável
    end
    if mod(k,50)==0, fprintf('    %d/%d\n', k, N); end
end

%% ---------------- RELATÓRIO ----------------
ok = metr.stable;
fprintf('\n  ============ ROBUSTEZ PARAMÉTRICA (N=%d, 1σ = CR%%) ============\n', N);
fprintf('  Rodadas estáveis: %d/%d (%.1f%%)\n', nnz(ok), N, 100*mean(ok));
rep = @(nm,v,u) fprintf('  %-22s média %7.2f | desvio %6.2f | pior %7.2f  %s\n', ...
                        nm, mean(v(ok)), std(v(ok)), max(v(ok)), u);
rep('sobressinal altitude', metr.OS_h,   '%');
rep('sobressinal velocid.', metr.OS_u,   '%');
rep('t_s altitude (5%)',    metr.ts_h,   's');
rep('t_s velocidade (5%)',  metr.ts_u,   's');
rep('erro regime h',        metr.e_h,    'm');
rep('erro regime u',        metr.e_u,    'm/s');
rep('erro regime psi',      metr.e_psi,  'deg');
rep('rolagem máxima',       metr.phi_max,'deg');
rep('saturação PWM',        metr.sat_pct,'%');


%% ---------------- DIAGNÓSTICO: o que caracteriza as rodadas instáveis? ----------------
nm = {'Jx','Jy','Jz','Jxz','kT1','kT2','kT3','kT4','kQ1','kQ2','kQ3','kQ4','cp','cq','cr'};
if any(~ok)
    dz = mean(Z(~ok,:),1) - mean(Z(ok,:),1);     % deslocamento médio (em σ) nas instáveis
    [~,ord] = sort(abs(dz),'descend');
    fprintf('\n  Parâmetros mais associados à instabilidade (desvio médio em sigma, instáveis - estáveis):\n');
    for i = ord(1:5), fprintf('    %-4s %+5.2f sigma\n', nm{i}, dz(i)); end
end

%% ---------------- FIGURA: dispersão das métricas ----------------
f = figure('Position',[40 40 1100 380]); set(f,'Color','w','DefaultAxesFontSize',12);
try, f.Theme='light'; catch, end
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
nexttile; histogram(metr.OS_u(ok),20,'FaceColor',[.3 .3 .3]); grid on;
  xlabel('sobressinal de velocidade [%]'); ylabel('rodadas'); title('(a)','FontWeight','bold');
nexttile; histogram(metr.ts_h(ok),20,'FaceColor',[.3 .3 .3]); grid on;
  xlabel('t_s altitude [s]'); title('(b)','FontWeight','bold');
nexttile; histogram(metr.phi_max(ok),20,'FaceColor',[.3 .3 .3]); grid on;
  xlabel('rolagem máxima [°]'); title('(c)','FontWeight','bold');
sgtitle(sprintf('Monte Carlo paramétrico (N=%d, 1\\sigma = CR%%): dispersão das métricas de malha fechada', N));
exportgraphics(f, fullfile(paths.images,'mc_robustness.png'), 'BackgroundColor','white','Resolution',200);
fprintf('\n  Figura salva: %s\n', fullfile(paths.images,'mc_robustness.png'));

%% ======================= HELPERS =======================
function os = overshoot(y, yref, t, t0)
    % sobressinal (%) após o degrau em t0, relativo à amplitude do degrau
    m = t>=t0; y0 = y(find(m,1)); os = 100*max(0, (max(y(m))-yref)/(yref-y0+eps));
end
function ts = settle(y, yref, t, t0, tol)
    % tempo de acomodação: último instante fora da faixa ±tol·|degrau|
    m = t>=t0; y0 = y(find(m,1)); band = tol*abs(yref-y0);
    idx = find(abs(y-yref)>band & m, 1, 'last');
    if isempty(idx), ts = 0; else, ts = t(idx)-t0; end
end
