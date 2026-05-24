function compare_models(mode, varargin)
%COMPARE_MODELS  Utilitário unificado de comparação dos modelos da planta VTOL.
%
% Substitui os 3 scripts antigos:
%   - compare_v4_vs_script.m
%   - validate_scenarios.m
%   - validate_lin_vs_nl.m
%
% Uso:
%   compare_models('log')        % v4.slx vs vtol_dynamics.m vs log medido
%                                %  Usa janela definida em setup_quad_v4.m
%                                %  Plot 2x3: p,q,r + acc_x,acc_y,acc_z
%
%   compare_models('scenarios')  % v4.slx vs vtol_dynamics.m em 5 cenários PWM
%                                %  sintéticos (hover, step, sin, sweep, random)
%                                %  Não usa o log. Tabela de erros + plot do pior
%
%   compare_models('linear')     % Linear vs NL.m vs NL.slx
%                                %  fsolve para achar trim PWM real
%                                %  Aplica perturbações pequenas, compara
%
%   compare_models('all')        % Roda os 3 modos em sequência
%
% NOTAS:
%   - Os 3 modos usam função vtol_dynamics.m (modo 17 estados com lag).
%   - O .slx alvo é sempre quad_model_v4.
%   - Modelos de motor: spline Akima da bancada (mesma usada no .slx).
%   - Massas usam os mesmos valores dos scripts originais (ainda 1.6011 kg;
%     atualizar para 2.20 quando task #72 propagar para v4.slx também).

if nargin < 1, mode = 'log'; end
% Adicionar raiz do projeto ao path (sobe um nível de 5_validation/)
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();  % adiciona todas subpastas ativas ao path

switch lower(mode)
    case 'log',       compare_vs_log();
    case 'scenarios', compare_scenarios();
    case 'linear',    compare_linear();
    case 'all'
        compare_vs_log();
        compare_scenarios();
        compare_linear();
    otherwise
        error('compare_models:badMode', ...
            'Modo desconhecido: "%s". Use: log | scenarios | linear | all', mode);
end
end


%% ========================================================================
%  MODO 1: compare_vs_log
%  v4.slx vs vtol_dynamics.m vs log experimental (medido pelo IMU)
%  ========================================================================
function compare_vs_log()
fprintf('\n========== MODO: log (v4.slx vs vtol_dynamics.m vs medido) ==========\n');

bdclose all;
evalin('base', 'clearvars');
evalin('base', 'run setup_quad_v4');
gv = @(n) evalin('base', n);

P_J     = gv('P_J');
t_sim   = gv('t_sim');
ref     = gv('ref');
phi0    = gv('phi0');
theta0  = gv('theta0');
psi0    = gv('psi0');
tau_m   = gv('tau_motor');

% --- 1) Rodar v4.slx capturando estados via outports temporários ---
[state_v4, sub_specs] = run_v4_with_extra_outports(t_sim, ref.time, true);

% --- 2) Rodar vtol_dynamics.m (17 estados, mesmo tau do .slx) ---
[func_T, func_Q] = build_motor_models();
constants_sim = build_constants(tau_m);

y0 = build_y0_17(P_J, ref.pwm(1,:), phi0, theta0, psi0, func_T, func_Q);

ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_s, y_s] = ode45(@(t,y) vtol_dynamics(t, y, P_J, ref.time, ref.pwm, ...
    func_T, func_Q, constants_sim), ref.time, y0, ode_opts);
y_script = interp1(t_s, y_s, ref.time, 'linear', 'extrap');

state_script = unpack_states_17(y_script);
% Adicionar derivadas u_dot,v_dot,w_dot via vtol_dynamics
[udot, vdot, wdot] = compute_dots(P_J, ref.time, ref.pwm, y_script, ...
    func_T, func_Q, constants_sim);
state_script.udot = udot;
state_script.vdot = vdot;
state_script.wdot = wdot;

% --- 3) Converter w_dot em aceleração específica AccZ (mesma convenção IMU) ---
acc_v4     = build_acc_from_state(state_v4,     constants_sim.g);
acc_script = build_acc_from_state(state_script, constants_sim.g);
acc_meas   = struct('x', ref.acc(:,1), 'y', ref.acc(:,2), 'z', ref.acc(:,3));

% --- 4) Tabela de erros (v4 vs script) ---
print_error_table(state_v4, state_script, acc_v4, acc_script);

% --- 5) Plot 3x3: pqr + atitude + aceleracoes ---
sim_dict = struct( ...
    'v4',     pack_sim(state_v4, acc_v4), ...
    'script', pack_sim(state_script, acc_script));
meas_dict = pack_meas(ref.pqr, struct2acc(acc_meas), ref.att_deg);

title_str = sprintf('Comparacao v4 vs script vs medido | janela %g-%g s', ...
    ref.time_abs(1), ref.time_abs(end));

plot_3x3_full(ref.time, sim_dict, meas_dict, title_str, ...
    fullfile(setup_paths().images, 'compare_v4_vs_script.png'));
end


%% ========================================================================
%  MODO 2: compare_scenarios
%  v4.slx vs vtol_dynamics.m em 5 cenários PWM sintéticos (sem log)
%  ========================================================================
function compare_scenarios()
fprintf('\n========== MODO: scenarios (v4 vs script em cenarios PWM) ==========\n');

bdclose all; clear; clc;

% --- Parâmetros (P0 - chute inicial) ---
P_J = default_P0();
[func_T, func_Q] = build_motor_models();
tau_m = 0.05;
constants = build_constants(tau_m);
phi0 = 0; theta0 = 0; psi0 = 0;

% --- Cenários ---
dt = 0.01; T_end = 3;
t_vec = (0:dt:T_end)';
N = numel(t_vec);
scenarios = build_scenarios(t_vec, N);

% --- Setup v4.slx com outports temporários ---
load_system('quad_model_v4.slx');
model = 'quad_model_v4';
set_param(model, 'SaveOutput','on', 'SaveFormat','Dataset');

sub_specs = {
    'Rotational Dynamics',    1, 'p';
    'Rotational Dynamics',    2, 'q';
    'Rotational Dynamics',    3, 'r';
    'pqr2euler1',             1, 'phi';
    'pqr2euler1',             2, 'theta';
    'pqr2euler1',             3, 'psi';
    'Translational Dynamics', 1, 'u';
    'Translational Dynamics', 2, 'v';
    'Translational Dynamics', 3, 'w';
};
added = add_outports(model, sub_specs);

% Workspace: setup que o quad_model_v4 lê
P_estimated = P_J_to_simulink(P_J); %#ok<NASGU>
assignin('base', 'P_estimated', P_estimated);
assignin('base', 'tau_motor', tau_m);   % MotorLag T1..T4, Q1..Q4 lêem isso

% --- Loop principal ---
results = struct();
fprintf('\n%-40s  %-7s  %-7s  %-7s  %-7s\n', ...
    'Cenario', 'max|p|', 'max|q|', 'max|r|', 'max|w|');
fprintf('%s\n', repmat('-', 1, 80));

z_ts = timeseries(zeros(N,1), t_vec);

for s = 1:numel(scenarios)
    sc = scenarios(s);
    assign_pwm_workspace(sc.pwm, t_vec);
    for nm = {'p_out','q_out','r_out','phi_out','theta_out', ...
              'psi_out','u_dot_out','v_dot_out','w_dot_out'}
        assignin('base', nm{1}, z_ts);
    end

    pwm0 = sc.pwm(1,:);
    T_eff_init = zeros(4,1); Q_eff_init = zeros(4,1);
    for i = 1:4
        T_eff_init(i) = func_T(pwm0(i));
        Q_eff_init(i) = func_Q(pwm0(i));
    end
    assignin('base', 'T_eff_init', T_eff_init);
    assignin('base', 'Q_eff_init', Q_eff_init);

    set_param(model, 'StopTime', num2str(T_end));
    out  = sim(model);
    yout = out.yout;

    state_v4 = extract_v4_state(yout, sub_specs, t_vec);

    % vtol_dynamics no script (mesmas ICs do v4: p=0, q=-0.1, r=-0.1)
    y0_9  = [0; -0.1; -0.1; phi0; theta0; psi0; 0; 0; 0];
    k_T = P_J(5:8); k_Q = P_J(9:12);
    y0_17 = [y0_9; T_eff_init .* k_T; Q_eff_init .* k_Q];
    ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t_s, y_s] = ode45(@(t,y) vtol_dynamics(t, y, P_J, t_vec, sc.pwm, ...
        func_T, func_Q, constants), t_vec, y0_17, ode_opts);
    y_sc = interp1(t_s, y_s, t_vec, 'linear','extrap');
    state_sc = unpack_states_17(y_sc);

    % Erros relativos (% da amplitude do v4)
    fields = {'p','q','r','phi','theta','psi','u','v','w'};
    err_pct = zeros(1, numel(fields));
    for k = 1:numel(fields)
        n = fields{k};
        d = state_v4.(n) - state_sc.(n);
        rng = max(state_v4.(n)) - min(state_v4.(n));
        if rng < 1e-9, rng = max(abs(state_v4.(n))) + 1e-9; end
        err_pct(k) = 100 * max(abs(d)) / rng;
    end

    fprintf('%-40s  %6.3f%%  %6.3f%%  %6.3f%%  %6.3f%%\n', sc.name, ...
        err_pct(1), err_pct(2), err_pct(3), err_pct(9));

    results(s).name     = sc.name;
    results(s).err_pct  = err_pct;
    results(s).state_v4 = state_v4;
    results(s).state_sc = state_sc;
    results(s).t        = t_vec;
end

% Cleanup
remove_outports(added);
bdclose all;

% Resumo completo
print_scenarios_summary(results);

% Plot do pior cenário
plot_worst_scenario(results, ...
    fullfile(setup_paths().images,'validate_scenarios_worst.png'));
end


%% ========================================================================
%  MODO 3: compare_linear
%  Linear vs NL.m vs NL.slx (com fsolve para trim)
%  ========================================================================
function compare_linear()
fprintf('\n========== MODO: linear (linear vs NL.m vs NL.slx) ==========\n');

bdclose all; clear; clc;

% Carrega modelo linear pré-computado (linearize.m gera em outputs/)
lm_path = fullfile(setup_paths().outputs, 'linear_model.mat');
if ~exist(lm_path, 'file')
    error(['%s nao encontrado. Rode linearize.m primeiro ', ...
        'para gerar A, B, x0, u0, f0.'], lm_path);
end
lm = load(lm_path);

% Pegar P_J consistente com o ponto de linearização
P_J = default_P0();
[func_T, func_Q] = build_motor_models();
constants = build_constants(0.05);

% Re-trim numérico para garantir consistência atual (PWM_f ≠ PWM_r por dx_cg)
fprintf('\nBuscando trim numérico (fsolve)...\n');
[pwm_trim, success] = find_trim_pwm(P_J, func_T, func_Q, constants);
if ~success
    warning('Trim numérico não convergiu. Usando PWM uniforme como fallback.');
    pwm_trim = [1500 1500 1500 1500];
end
fprintf('PWM trim encontrado: %s us\n', mat2str(round(pwm_trim)));

% Cenários de perturbação a partir do trim
dt = 0.01; T_end = 2.0;
t_vec = (0:dt:T_end)';
N = numel(t_vec);

perturb_scenarios = build_perturbation_scenarios(pwm_trim, t_vec, N);

% Roda cada cenário em 3 modelos
results = struct();
for s = 1:numel(perturb_scenarios)
    sc = perturb_scenarios(s);
    fprintf('\nCenario: %s\n', sc.name);

    % NL.m (script)
    y0 = build_y0_17_at_trim(P_J, pwm_trim, func_T, func_Q);
    [t_s, y_s] = ode45(@(t,y) vtol_dynamics(t, y, P_J, t_vec, sc.pwm, ...
        func_T, func_Q, constants), t_vec, y0, odeset('RelTol',1e-8,'AbsTol',1e-10));
    y_nlm = interp1(t_s, y_s, t_vec, 'linear', 'extrap');

    % Linear
    y_lin = simulate_linear(lm, t_vec, sc.pwm, pwm_trim);

    % NL.slx
    y_slx = simulate_v4_for_perturbation(sc.pwm, t_vec, P_J, func_T, func_Q);

    results(s).name  = sc.name;
    results(s).t     = t_vec;
    results(s).y_nlm = y_nlm;
    results(s).y_lin = y_lin;
    results(s).y_slx = y_slx;
end

plot_linear_comparison(results, fullfile(setup_paths().images,'validate_lin_vs_nl.png'));
end


%% ========================================================================
%  HELPERS COMPARTILHADOS
%  ========================================================================

function [pwm_b, T_b, Q_b] = bench_table()
% Tabela de bancada (consistente com o que setup_quad_v4 usa)
pwm_b = [1000; 1200; 1400; 1600; 1800; 2000];
% Thrust em kg-força (converter via g)
thrust_g = [0; 143; 328; 532; 784; 843];
T_b = thrust_g * 9.80665 / 1000;     % N
Q_b = [0; 0.034; 0.070; 0.115; 0.171; 0.176];  % Nm
end

function [func_T, func_Q] = build_motor_models()
% Delega para motor_models.m (2_model/) — usa tabela centralizada
[func_T, func_Q] = motor_models();
end

function c = build_constants(tau_motor)
% Massa e gravidade vêm de parameters.m (fonte única — slide oficial)
p = parameters();
c = struct('m', p.m, 'g', p.g, 'tau_motor', tau_motor);
% ATENÇÃO: o v4.slx ainda usa massa 1.6011 no bloco Massa. Para máxima
% consistência com .slx, manter v4.slx atualizado também (task #72 parcial).
end

function P_J = default_P0()
% Delega para parameters.m (2_model/) — fonte única de verdade
p = parameters();
P_J = p.P0_J;
end

function y0 = build_y0_17(P_J, pwm0, phi0, theta0, psi0, func_T, func_Q)
% 17 estados: [p q r phi theta psi u v w T_eff(1:4) Q_eff(1:4)]
%
% ICs de p,q,r: [0; -0.1; -0.1] para CASAR com os integradores
% hardcoded do quad_model_v4.slx. Mudar aqui se v4.slx for atualizado.
T_eff = zeros(4,1); Q_eff = zeros(4,1);
k_T = P_J(5:8); k_Q = P_J(9:12);
for i = 1:4
    T_eff(i) = k_T(i) * func_T(pwm0(i));
    Q_eff(i) = k_Q(i) * func_Q(pwm0(i));
end
y0 = [0; -0.1; -0.1; phi0; theta0; psi0; 0; 0; 0; T_eff; Q_eff];
end

function y0 = build_y0_17_at_trim(P_J, pwm_trim, func_T, func_Q)
y0 = build_y0_17(P_J, pwm_trim, 0, 0, 0, func_T, func_Q);
end

function s = unpack_states_17(y)
s = struct('p', y(:,1), 'q', y(:,2), 'r', y(:,3), ...
           'phi', y(:,4), 'theta', y(:,5), 'psi', y(:,6), ...
           'u', y(:,7), 'v', y(:,8), 'w', y(:,9));
end

function [udot, vdot, wdot] = compute_dots(P_J, t_vec, pwm, y, func_T, func_Q, constants)
N = numel(t_vec);
udot = zeros(N,1); vdot = zeros(N,1); wdot = zeros(N,1);
for i = 1:N
    dy = vtol_dynamics(t_vec(i), y(i,:)', P_J, t_vec, pwm, ...
        func_T, func_Q, constants);
    udot(i) = dy(7); vdot(i) = dy(8); wdot(i) = dy(9);
end
end

function acc = build_acc_from_state(state, g)
% Converte u_dot/v_dot/w_dot integrados → aceleração específica IMU
acc = struct();
acc.x = state.udot;
acc.y = state.vdot;
acc.z = state.wdot - (state.q .* state.u - state.p .* state.v) ...
        - g * cos(state.phi) .* cos(state.theta);
end

function added = add_outports(model, sub_specs)
added = cell(1, size(sub_specs,1));
for k = 1:size(sub_specs,1)
    op = sprintf('%s/OUT_%s', model, sub_specs{k,3});
    add_block('simulink/Sinks/Out1', op, 'Position', [1100 1000+30*k 1130 1020+30*k]);
    add_line(model, sprintf('%s/%d', sub_specs{k,1}, sub_specs{k,2}), ...
        sprintf('OUT_%s/1', sub_specs{k,3}), 'autorouting','on');
    added{k} = op;
end
end

function remove_outports(added)
for k = 1:numel(added)
    try, delete_block(added{k}); catch, end %#ok<NOSEM>
end
end

function state = extract_v4_state(yout, sub_specs, t_target)
state = struct();
for k = 1:yout.numElements
    el = yout{k};
    bp = el.BlockPath.getBlock(1);
    for j = 1:size(sub_specs,1)
        if endsWith(bp, ['/OUT_' sub_specs{j,3}])
            state.(sub_specs{j,3}) = interp1(el.Values.Time, ...
                el.Values.Data, t_target, 'linear', 'extrap');
        end
    end
end
end

function [state, sub_specs] = run_v4_with_extra_outports(t_sim, t_target, include_dots)
load_system('quad_model_v4.slx');
model = 'quad_model_v4';
set_param(model, 'StopTime', num2str(t_sim));
set_param(model, 'SaveOutput','on', 'SaveFormat','Dataset');

sub_specs = {
    'Rotational Dynamics',    1, 'p';
    'Rotational Dynamics',    2, 'q';
    'Rotational Dynamics',    3, 'r';
    'pqr2euler1',             1, 'phi';
    'pqr2euler1',             2, 'theta';
    'pqr2euler1',             3, 'psi';
    'Translational Dynamics', 1, 'u';
    'Translational Dynamics', 2, 'v';
    'Translational Dynamics', 3, 'w';
};
if include_dots
    sub_specs = [sub_specs; {
        'Translational Dynamics', 4, 'udot';
        'Translational Dynamics', 5, 'vdot';
        'Translational Dynamics', 6, 'wdot'}];
end

added  = add_outports(model, sub_specs);
simOut = sim(model);
state  = extract_v4_state(simOut.yout, sub_specs, t_target);

remove_outports(added);
bdclose all;
end

function assign_pwm_workspace(pwm, t_vec)
for i = 1:4
    assignin('base', sprintf('pwm%d_in', i), timeseries(pwm(:,i), t_vec));
end
end

function scenarios = build_scenarios(t_vec, N)
scenarios = struct();
scenarios(1).name = 'Hover constante (PWM=1620)';
scenarios(1).pwm  = 1620 * ones(N, 4);

pwm = 1620 * ones(N, 4); pwm(t_vec >= 0.5, 1) = 1700;
scenarios(2).name = 'Step motor 1: 1620 -> 1700 em t=0.5s';
scenarios(2).pwm  = pwm;

pwm = zeros(N, 4);
for i = 1:4, pwm(:, i) = 1620 + 100 * sin(2*pi*1*t_vec + (i-1)*pi/2); end
scenarios(3).name = 'Sin defasadas (A=100us, f=1Hz)';
scenarios(3).pwm  = pwm;

pwm = zeros(N, 4);
for i = 1:4, pwm(:, i) = 1300 + 600 * t_vec/t_vec(end) + (i-1)*20; end
pwm = max(1000, min(2000, pwm));
scenarios(4).name = 'Sweep 1300->1900';
scenarios(4).pwm  = pwm;

rng(42);
raw = 1400 + 450*rand(N, 4);
for i = 1:4, raw(:,i) = movmean(raw(:,i), 20); end
scenarios(5).name = 'Random filtrado [1400, 1850]';
scenarios(5).pwm  = raw;
end

function perturb = build_perturbation_scenarios(pwm_trim, t_vec, N)
% Pequenas perturbações ±30 us a partir do trim
perturb = struct();
perturb(1).name = 'Roll step (+30us motor 1, -30us motor 2)';
pwm = repmat(pwm_trim, N, 1);
pwm(t_vec >= 0.5, 1) = pwm_trim(1) + 30;
pwm(t_vec >= 0.5, 2) = pwm_trim(2) - 30;
perturb(1).pwm = pwm;

perturb(2).name = 'Pitch step (+30us motores frente, -30us trás)';
pwm = repmat(pwm_trim, N, 1);
pwm(t_vec >= 0.5, [1 3]) = pwm(t_vec >= 0.5, [1 3]) + 30;
pwm(t_vec >= 0.5, [2 4]) = pwm(t_vec >= 0.5, [2 4]) - 30;
perturb(2).pwm = pwm;

perturb(3).name = 'Yaw step (+30us CW, -30us CCW)';
pwm = repmat(pwm_trim, N, 1);
pwm(t_vec >= 0.5, [1 2]) = pwm(t_vec >= 0.5, [1 2]) + 30;
pwm(t_vec >= 0.5, [3 4]) = pwm(t_vec >= 0.5, [3 4]) - 30;
perturb(3).pwm = pwm;

perturb(4).name = 'Throttle bump (+50us em todos)';
pwm = repmat(pwm_trim, N, 1);
pwm(t_vec >= 0.5, :) = pwm(t_vec >= 0.5, :) + 50;
perturb(4).pwm = pwm;
end

function [pwm_trim, success] = find_trim_pwm(P_J, func_T, func_Q, constants)
% Acha PWM em hover steady. Por causa de dx_cg/dy_cg, PWM_front ≠ PWM_rear.
m = constants.m; g = constants.g;
total_thrust_needed = m * g;

% Chute inicial: PWM uniforme tal que soma ~= mg
chute_pwm = 1500;
chute = [chute_pwm chute_pwm chute_pwm chute_pwm];

opts = optimoptions('fsolve','Display','off','TolFun',1e-8);
try
    [pwm_trim, ~, exitflag] = fsolve(@(pwm) trim_residuals(pwm, P_J, ...
        func_T, func_Q, total_thrust_needed), chute, opts);
    success = (exitflag > 0);
catch ME
    fprintf('fsolve falhou: %s\n', ME.message);
    pwm_trim = chute;
    success = false;
end
end

function r = trim_residuals(pwm, P_J, func_T, func_Q, T_needed)
% Resíduos: [Mx, My, Mz, T_total - mg] em hover. dx_cg, dy_cg de P_J
k_T = P_J(5:8); k_Q = P_J(9:12);
dx_cg = P_J(19); dy_cg = P_J(20);

T = zeros(4,1); Q = zeros(4,1);
for i = 1:4
    T(i) = k_T(i) * func_T(pwm(i));
    Q(i) = k_Q(i) * func_Q(pwm(i));
end

Lx_r = 0.232 - dy_cg;  Lx_l = 0.232 + dy_cg;
Ly_f = 0.311185 - dx_cg;  Ly_r = 0.342865 + dx_cg;

Mx = -(Lx_r*T(1) - Lx_l*T(2) - Lx_l*T(3) + Lx_r*T(4));
My =  (Ly_f*T(1) - Ly_r*T(2) + Ly_f*T(3) - Ly_r*T(4));
Mz =   Q(1) + Q(2) - Q(3) - Q(4);
T_total = sum(T);

r = [Mx; My; Mz; T_total - T_needed];
end

function y = simulate_linear(lm, t_vec, pwm, pwm_trim)
% Simula modelo linear estado-espaço: dx_dot = A·dx + B·du
%   x = x0 + dx;  u = u0 + du
du = pwm - pwm_trim;  % perturbação em PWM
N = numel(t_vec); dt = t_vec(2)-t_vec(1);
nx = size(lm.A, 1);

dx = zeros(N, nx);
dx(1,:) = 0;
for k = 1:N-1
    % Euler explícito (simples, suficiente pra perturbações pequenas em janela curta)
    dx(k+1,:) = dx(k,:) + dt * (lm.A * dx(k,:)' + lm.B * du(k,:)')';
end
y = dx + lm.x0(:)';  % adiciona x0 do trim
end

function y = simulate_v4_for_perturbation(pwm, t_vec, P_J, func_T, func_Q)
% Roda o quad_model_v4.slx para a perturbação dada, retorna estados.
load_system('quad_model_v4.slx');
model = 'quad_model_v4';
set_param(model, 'SaveOutput','on', 'SaveFormat','Dataset');
set_param(model, 'StopTime', num2str(t_vec(end)));

sub_specs = {
    'Rotational Dynamics',    1, 'p';
    'Rotational Dynamics',    2, 'q';
    'Rotational Dynamics',    3, 'r';
    'pqr2euler1',             1, 'phi';
    'pqr2euler1',             2, 'theta';
    'pqr2euler1',             3, 'psi';
    'Translational Dynamics', 1, 'u';
    'Translational Dynamics', 2, 'v';
    'Translational Dynamics', 3, 'w';
};
added = add_outports(model, sub_specs);

assign_pwm_workspace(pwm, t_vec);
N = numel(t_vec);
z_ts = timeseries(zeros(N,1), t_vec);
for nm = {'p_out','q_out','r_out','phi_out','theta_out','psi_out', ...
          'u_dot_out','v_dot_out','w_dot_out'}
    assignin('base', nm{1}, z_ts);
end

T_eff_init = zeros(4,1); Q_eff_init = zeros(4,1);
for i = 1:4
    T_eff_init(i) = func_T(pwm(1,i));
    Q_eff_init(i) = func_Q(pwm(1,i));
end
assignin('base', 'T_eff_init', T_eff_init);
assignin('base', 'Q_eff_init', Q_eff_init);

P_estimated = P_J_to_simulink(P_J); %#ok<NASGU>
assignin('base', 'P_estimated', P_estimated);
assignin('base', 'tau_motor', 0.05);   % MotorLag T1..T4, Q1..Q4 lêem isso

simOut = sim(model);
state  = extract_v4_state(simOut.yout, sub_specs, t_vec);

remove_outports(added);
bdclose all;

y = [state.p, state.q, state.r, state.phi, state.theta, state.psi, ...
     state.u, state.v, state.w];
end


%% ========================================================================
%  PLOTTING HELPERS
%  ========================================================================

function print_error_table(state_v4, state_script, acc_v4, acc_script)
fprintf('\n=========================================================\n');
fprintf(' v4.slx vs vtol_dynamics.m | max|err| e %% relativo a amplitude\n');
fprintf('=========================================================\n');
for nm = {'p','q','r'}
    n = nm{1};
    err = state_v4.(n) - state_script.(n);
    rng = max(state_v4.(n)) - min(state_v4.(n));
    if rng < 1e-9, rng = 1e-9; end
    fprintf('  %-6s: max|err|=%9.5g   %.3f %%\n', n, max(abs(err)), 100*max(abs(err))/rng);
end
for nm = {'x','y','z'}
    n = nm{1};
    err = acc_v4.(n) - acc_script.(n);
    rng = max(acc_v4.(n)) - min(acc_v4.(n));
    if rng < 1e-9, rng = 1e-9; end
    fprintf('  acc_%s: max|err|=%9.5g   %.3f %%\n', n, max(abs(err)), 100*max(abs(err))/rng);
end
end

function sim = pack_sim(state, acc)
% Atitude vai em GRAUS pra plot ser legível
sim = struct('p', state.p, 'q', state.q, 'r', state.r, ...
             'phi',   rad2deg(state.phi), ...
             'theta', rad2deg(state.theta), ...
             'psi',   rad2deg(state.psi), ...
             'ax', acc.x, 'ay', acc.y, 'az', acc.z);
end

function meas = pack_meas(pqr, acc, att_deg)
% att_deg: matriz Nx3 [roll, pitch, yaw] em graus (vem do EKF do log)
meas = struct('p', pqr(:,1), 'q', pqr(:,2), 'r', pqr(:,3), ...
              'phi',   att_deg(:,1), ...
              'theta', att_deg(:,2), ...
              'psi',   att_deg(:,3), ...
              'ax', acc.x, 'ay', acc.y, 'az', acc.z);
end

function acc = struct2acc(s)
acc = struct('x', s.x, 'y', s.y, 'z', s.z);
end

function plot_3x3_full(t, sim_dict, meas, title_str, save_path)
% Plot 3x3:
%   linha 1: p, q, r          [rad/s]
%   linha 2: phi, theta, psi  [deg]   ◄── atitude
%   linha 3: ax, ay, az       [m/s²]

fig = figure('Position',[100 100 1600 900], 'Color','w');
labels = {'p','q','r','phi','theta','psi','ax','ay','az'};
units  = {'rad/s','rad/s','rad/s','deg','deg','deg','m/s^2','m/s^2','m/s^2'};

names = fieldnames(sim_dict);
colors = {'b-','r--','m:'};

for k = 1:9
    subplot(3,3,k); hold on; grid on;
    if ~isempty(meas)
        plot(t, meas.(labels{k}), 'Color',[0 0.9 0.3], ...
            'LineWidth',1.4, 'DisplayName','medido');
    end
    for j = 1:numel(names)
        plot(t, sim_dict.(names{j}).(labels{k}), colors{j}, ...
            'LineWidth',1.5, 'DisplayName', names{j});
    end
    title(labels{k}, 'FontWeight','bold');
    xlabel('t [s]'); ylabel(units{k});
    if k == 1, legend('Location','best'); end
end
sgtitle(title_str, 'FontWeight','bold');

if ~isempty(save_path)
    out_dir = fileparts(save_path);
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    saveas(fig, save_path);
    fprintf('Figura salva: %s\n', save_path);
end
end

function print_scenarios_summary(results)
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf(' RESUMO COMPLETO: erros maximos (%% da amplitude do v4) por estado\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('%-40s', 'Cenario');
for nm = {'p','q','r','phi','theta','psi','u','v','w'}
    fprintf(' %-5s', nm{1});
end
fprintf('\n');
for s = 1:numel(results)
    fprintf('%-40s', results(s).name);
    for k = 1:9
        fprintf(' %4.2f%%', results(s).err_pct(k));
    end
    fprintf('\n');
end

mean_errs = arrayfun(@(r) mean(r.err_pct), results);
[~, idx_worst] = max(mean_errs);
fprintf('\nPior cenario: %s (erro medio = %.2f%%)\n', ...
    results(idx_worst).name, mean_errs(idx_worst));
end

function plot_worst_scenario(results, save_path)
mean_errs = arrayfun(@(r) mean(r.err_pct), results);
[~, idx_worst] = max(mean_errs);
r = results(idx_worst);

fig = figure('Position',[100 100 1400 600], 'Color','w');
labels = {'p','q','r','phi','theta','psi','u','v','w'};
units  = {'rad/s','rad/s','rad/s','rad','rad','rad','m/s','m/s','m/s'};
for k = 1:9
    subplot(3,3,k); hold on; grid on;
    plot(r.t, r.state_v4.(labels{k}), 'b-', 'LineWidth',1.6, 'DisplayName','v4.slx');
    plot(r.t, r.state_sc.(labels{k}), 'r--', 'LineWidth',1.0, 'DisplayName','vtol\_dynamics.m');
    title(labels{k}); xlabel('t [s]'); ylabel(units{k});
    if k == 1, legend('Location','best'); end
end
sgtitle(sprintf('Pior cenario: %s', r.name));

if ~isempty(save_path)
    out_dir = fileparts(save_path);
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    saveas(fig, save_path);
    fprintf('Figura salva: %s\n', save_path);
end
end

function plot_linear_comparison(results, save_path)
n = numel(results);
fig = figure('Position',[100 100 1600 350*n], 'Color','w');
state_names = {'p','q','r','phi','theta','psi','u','v','w'};
state_idx   = [1 2 3 4 5 6 7 8 9];

for s = 1:n
    r = results(s);
    for k = 1:9
        subplot(n, 9, (s-1)*9 + k); hold on; grid on;
        plot(r.t, r.y_nlm(:, state_idx(k)), 'b-',  'LineWidth',1.5, 'DisplayName','NL .m');
        plot(r.t, r.y_lin(:, state_idx(k)), 'r--', 'LineWidth',1.2, 'DisplayName','Linear');
        plot(r.t, r.y_slx(:, state_idx(k)), 'm:',  'LineWidth',1.2, 'DisplayName','NL .slx');
        if s == 1, title(state_names{k}); end
        if k == 1, ylabel(r.name, 'FontWeight','bold'); end
        if s == n, xlabel('t [s]'); end
        if s == 1 && k == 1, legend('Location','best'); end
    end
end
sgtitle('Comparacao Linear vs NL.m vs NL.slx (perturbacoes do trim)');

if ~isempty(save_path)
    out_dir = fileparts(save_path);
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    saveas(fig, save_path);
    fprintf('Figura salva: %s\n', save_path);
end
end
