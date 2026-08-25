%COMPARE_V4_VS_VALIDATE  Compara saídas de quad_model_v4.slx vs sim_window (.m)
%
% Roda os DOIS modelos com EXATAMENTE os mesmos parâmetros e janela:
%   1) quad_model_v4.slx via sim()  → simOut.*_sim
%   2) sim_window('full', ...)      → res.*
%
% Plota overlay medido / sim_v4 / sim_dotm e calcula:
%   - R² por sinal (formato igual identify_plant.m → print_R2)
%   - max |diff| entre v4 e .m  (se modelo está sincronizado, diff ~ 0)
%
% Pré-requisitos:
%   - Workspace populado (setup_quad_v4.m roda na entrada).
%   - quad_model_v4.slx atualizado (com Accelerometer Model + p_dot/q_dot/r_dot).

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% =====================================================================
%  CONFIG (mesmo do setup_quad_v4 — herdam o mesmo P_J + janela)
%  =====================================================================
% Pra mudar log/janela/P, edite o topo do setup_quad_v4.m

%% =====================================================================
%  1. Rodar setup + simular o Simulink
%  =====================================================================
fprintf('═════════════════════════════════════════════════════════════\n');
fprintf('  PASSO 1: Setup workspace e simular quad_model_v4.slx\n');
fprintf('═════════════════════════════════════════════════════════════\n');

run(fullfile(paths.simulink, 'setup_quad_v4.m'));

fprintf('\nSimulando quad_model_v4 (t_sim=%.2fs)...\n', t_sim);
tic;
simOut = sim('quad_model_v4', 'StopTime', num2str(t_sim));
fprintf('Simulink OK (%.2fs)\n', toc);

% Extrair timeseries (ToWorkspace blocks → Timeseries no simOut)
slx = struct();
slx.time = ref.time;
for sig = {'p_sim','q_sim','r_sim','phi_sim','theta_sim','psi_sim','accX_sim','accY_sim','accZ_sim'}
    if isprop(simOut, sig{1})
        ts = simOut.(sig{1});
        % Reamostrar pra grade do ref (interp linear)
        slx.(sig{1}) = interp1(ts.Time, ts.Data, ref.time, 'linear', 'extrap');
    else
        warning('Simulink output %s não encontrado.', sig{1});
        slx.(sig{1}) = nan(size(ref.time));
    end
end

%% =====================================================================
%  2. Rodar sim_window do .m com mesmo P_J e janela
%  =====================================================================
fprintf('\n═════════════════════════════════════════════════════════════\n');
fprintf('  PASSO 2: Simular sim_window (.m) com mesmo P_J\n');
fprintf('═════════════════════════════════════════════════════════════\n');

constants.m = mass;
constants.g = g_acc;

% Reconstruir inputs no formato esperado pelo sim_window
time_m  = ref.time;
pwm_m   = ref.pwm;
pqr_m   = ref.pqr;
att_m   = ref.att_deg;
acc_m   = ref.acc;

P_J = ref.P_J;

tic;
mres = sim_window('full', P_J, time_m, pwm_m, pqr_m, att_m, constants);
fprintf('sim_window OK (%.2fs)\n', toc);

% Adapt: graus
mres.phi   = mres.phi(:);
mres.theta = mres.theta(:);
mres.psi   = mres.psi(:);

%% =====================================================================
%  2b. Simular o MODELO LINEAR de HOVER (fonte única: linear_model.mat)
%      Entrada em FORÇAS:  δẋ = A·δx + B·δv,  δv = [T,M] - v0,  x = x0 + δx
%      Carrega o MESMO A,B,x0,u0=[mg;0;0;0] do projeto (gerado por linearize.m),
%      via load_linear_model — sem recalcular, sem divergir.
%
%      As forças [T,M](t) vêm do PWM logado pela ALOCAÇÃO (forces handle +
%      tabelas) — a mesma alocação do modelo NL.
%
%      NOTA: é o modelo de HOVER. Sobre esta janela manobrada ele diverge — em
%      especial no yaw, pois a força da janela foge do trim e o polo fraco do
%      yaw integra o offset. É o LIMITE DE VALIDADE do linear, não um erro.
%  =====================================================================
lin = struct();
[func_T_l, func_Q_l] = motor_models();
N  = numel(time_m);
dt = time_m(2) - time_m(1);

lm = load_linear_model(paths);     % avisa se stale vs P_identified.mat
A_lin = lm.A;  B_lin = lm.B;  x0_lin = lm.x0;  u0_lin = lm.u0;   % u0=[mg;0;0;0]
fprintf('\nModelo linear de hover carregado de linear_model.mat\n');

% ALOCAÇÃO: PWM logado → [T, Mx, My, Mz](t)  (mesma do modelo NL)
dyn_l = vtol_dynamics('get_handles');
[Tt_l, Mxt_l, Myt_l, Mzt_l] = dyn_l.forces(pwm_m, P_J, func_T_l, func_Q_l);
V_lin = [Tt_l, Mxt_l, Myt_l, Mzt_l];                   % Nx4

% IC: estados medidos em t=0 (mesma IC dos modelos NL)
x_ic = [pqr_m(1,:)'; deg2rad(att_m(1,:))'; 0; 0; 0];   % 9x1
dx   = x_ic - x0_lin;
dv   = V_lin - u0_lin';                                % Nx4 (broadcast)

% Discretização EXATA (ZOH) em vez de Euler explícito — o Euler a dt=0.1 com
% polos rápidos (Dp~5.7) atrasa/distorce. Ad=e^{A·dt}, Bd via exponencial aumentada.
Mexp = expm([A_lin, B_lin; zeros(4,9), zeros(4,4)] * dt);
Ad = Mexp(1:9,1:9);  Bd = Mexp(1:9,10:13);
DX = zeros(N, 9);  DX(1,:) = dx';
for k = 1:N-1
    DX(k+1,:) = (Ad * DX(k,:)' + Bd * dv(k,:)')';        % hover: c≈0, omitido
end
X = DX + x0_lin';                                       % estado completo

lin.p = X(:,1);  lin.q = X(:,2);  lin.r = X(:,3);
lin.phi = rad2deg(X(:,4));  lin.theta = rad2deg(X(:,5));  lin.psi = rad2deg(X(:,6));
lin.u = X(:,7);  lin.v = X(:,8);  lin.w = X(:,9);

% Acelerômetro: mesmo modelo de sensor aplicado aos estados lineares
k_T_l = P_J(5:8);
proj_l = parameters();
r_imu_l = proj_l.imu_offset;
T_total_l = zeros(N,1);
for k = 1:N
    for j = 1:4
        T_total_l(k) = T_total_l(k) + k_T_l(j) * func_T_l(pwm_m(k,j));
    end
end
p_dot_l = gradient(lin.p, dt);  q_dot_l = gradient(lin.q, dt);  r_dot_l = gradient(lin.r, dt);
[lin.accX, lin.accY, lin.accZ] = accelerometer_model( ...
    lin.p, lin.q, lin.r, lin.u, lin.v, lin.w, T_total_l/mass, ...
    p_dot_l, q_dot_l, r_dot_l, r_imu_l);
lin.ok = true;
fprintf('Linear (hover) OK.\n');

%% =====================================================================
%  3. Tabela de R² (formato idêntico ao identify_plant.m → print_R2)
%  =====================================================================
fprintf('\n═════════════════════════════════════════════════════════════\n');
fprintf('  PASSO 3: R² (mesmo formato do identify_plant)\n');
fprintf('═════════════════════════════════════════════════════════════\n');

R2 = @(y_e, y_s) 1 - sum((y_e - y_s).^2) / max(sum((y_e - mean(y_e)).^2), 1e-12);
t_win = [ref.time_abs(1), ref.time_abs(end)];

% v4.slx vs medido
fprintf('\n  [v4.slx]   Validação (%g-%gs):\n', t_win(1), t_win(2));
fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
    R2(pqr_m(:,1), slx.p_sim), R2(pqr_m(:,2), slx.q_sim), R2(pqr_m(:,3), slx.r_sim));
fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
    R2(acc_m(:,1), slx.accX_sim), R2(acc_m(:,2), slx.accY_sim), R2(acc_m(:,3), slx.accZ_sim));

% .m vs medido
fprintf('\n  [.m]       Validação (%g-%gs):\n', t_win(1), t_win(2));
fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
    R2(pqr_m(:,1), mres.p), R2(pqr_m(:,2), mres.q), R2(pqr_m(:,3), mres.r));
fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
    R2(acc_m(:,1), mres.accX), R2(acc_m(:,2), mres.accY), R2(acc_m(:,3), mres.accZ));

% linear vs medido (válido só perto do trim — divergência em manobra é esperada)
if lin.ok
    fprintf('\n  [linear]   Validação (%g-%gs):\n', t_win(1), t_win(2));
    fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
        R2(pqr_m(:,1), lin.p), R2(pqr_m(:,2), lin.q), R2(pqr_m(:,3), lin.r));
    fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
        R2(acc_m(:,1), lin.accX), R2(acc_m(:,2), lin.accY), R2(acc_m(:,3), lin.accZ));
end

% v4 vs .m (diff direto — se modelo está sync, diff ~ 0)
fprintf('\n  [v4 vs .m] max(|diff|) (idealmente ~0):\n');
fprintf('    p=%.4e q=%.4e r=%.4e\n', ...
    max(abs(slx.p_sim - mres.p)), max(abs(slx.q_sim - mres.q)), max(abs(slx.r_sim - mres.r)));
fprintf('    phi=%.4e theta=%.4e psi=%.4e\n', ...
    max(abs(slx.phi_sim - mres.phi)), max(abs(slx.theta_sim - mres.theta)), max(abs(slx.psi_sim - mres.psi)));
fprintf('    AccX=%.4e AccY=%.4e AccZ=%.4e\n', ...
    max(abs(slx.accX_sim - mres.accX)), max(abs(slx.accY_sim - mres.accY)), max(abs(slx.accZ_sim - mres.accZ)));

%% =====================================================================
%  4. Plot 3x3 com 3 cores (medido, sim_v4, sim_.m)
%  =====================================================================
fprintf('\n═════════════════════════════════════════════════════════════\n');
fprintf('  PASSO 4: Plots overlay (medido vs v4 vs .m)\n');
fprintf('═════════════════════════════════════════════════════════════\n');

t = ref.time;
fig = figure('Name', 'compare_v4_vs_validate', 'Position', [80 50 1400 900]);

% 5ª coluna = modelo linear (NaN se não rodou)
if lin.ok
    lin_cols = {lin.p, lin.q, lin.r, lin.accX, lin.accY, lin.accZ};
else
    nanv = nan(size(t));
    lin_cols = {nanv, nanv, nanv, nanv, nanv, nanv};
end

panels = {
    'p (rad/s)',     pqr_m(:,1), slx.p_sim,     mres.p,     lin_cols{1};
    'q (rad/s)',     pqr_m(:,2), slx.q_sim,     mres.q,     lin_cols{2};
    'r (rad/s)',     pqr_m(:,3), slx.r_sim,     mres.r,     lin_cols{3};
    'AccX (m/s²)',   acc_m(:,1), slx.accX_sim,  mres.accX,  lin_cols{4};
    'AccY (m/s²)',   acc_m(:,2), slx.accY_sim,  mres.accY,  lin_cols{5};
    'AccZ (m/s²)',   acc_m(:,3), slx.accZ_sim,  mres.accZ,  lin_cols{6};
};

for i = 1:6
    subplot(2, 3, i);
    plot(t, panels{i,2}, 'b-',  'LineWidth', 1.0, 'DisplayName', 'medido'); hold on;
    plot(t, panels{i,3}, 'r--', 'LineWidth', 1.3, 'DisplayName', 'sim v4');
    plot(t, panels{i,4}, 'g:',  'LineWidth', 1.3, 'DisplayName', 'sim .m');
    plot(t, panels{i,5}, 'm-.', 'LineWidth', 1.2, 'DisplayName', 'linear');
    hold off;
    title(panels{i,1});
    grid on;
    if i == 1, legend('Location','best'); end
    xlabel('t (s)');
end
sgtitle(sprintf('compare v4.slx vs sim\\_window vs linear  —  janela [%g, %g]s', t_win(1), t_win(2)));

save_path = fullfile(paths.images, 'compare_v4_vs_validate.png');
saveas(fig, save_path);
fprintf('\nFigura salva: %s\n', save_path);
fprintf('═════════════════════════════════════════════════════════════\n');
