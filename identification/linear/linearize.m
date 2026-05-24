% linearize.m - Linearizacao do modelo nao-linear VTOL
%
% Metodo: Jacobiano numerico via diferencas finitas centrais
%   A(i,j) = df_i/dx_j ~ [f(x0+h*ej) - f(x0-h*ej)] / (2h)
%   B(i,j) = df_i/du_j ~ [f(x0,u0+h*ej) - f(x0,u0-h*ej)] / (2h)
%
% Ponto de trim: hover (estados nulos, PWM onde T_total = m*g)
%
% Saida: salva A, B, x0, u0, eigenvalues em linear_model.mat
%
% Uso:
%   >> linearize          % calcula e salva
%   >> simulate           % carrega e compara NL vs linear

%% ========================================================================
%  1. CONFIGURACAO (identica ao simulate.m)
%  ========================================================================

% ---------- Selecao de parametros ----------
P = [0.063244; 0.250554; 0.116192; 0.001571;   ... % Jx, Jy, Jz, Jxz
     0.55; 0.45; 1.0; 0.75;                    ... % k_T1..k_T4
     0.55; 0.45; 1.0; 0.75;                    ... % k_Q1..k_Q4
     10; 5; 0.5;                                ... % Dp, Dq, Dr
     0.7; 1.4; 0.3;                            ... % Bp, Bq, Br
     0.0; 0.0;                                  ... % dx_cg, dy_cg
     -4.0; -4.0; -0.1; -0.5];                     % Xu_m, Yv_m, Zw_m, Bz

% Descomente para usar parametros identificados:
%id = load(fullfile(fileparts(mfilename('fullpath')), 'P_identified.mat'));
%P = id.P_final;

constants.m = 1.6011;
constants.g = 9.81;

%% ========================================================================
%  2. MODELOS DE REFERENCIA DOS MOTORES
%  ========================================================================
pwm_values_exp   = [1000; 1200; 1400; 1600; 1800; 2000];
thrust_grams_exp = [0; 143; 328; 532; 784; 843];
torque_Nm_exp    = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];
poly_degree = 3;

func_T_ref = create_thrust_model(pwm_values_exp, thrust_grams_exp, poly_degree);
func_Q_ref = create_torque_model(pwm_values_exp, torque_Nm_exp, poly_degree);

%% ========================================================================
%  3. PONTO DE TRIM (HOVER)
%  ========================================================================
% Estado de trim: tudo zero (hover nivelado, sem velocidade)
x0 = zeros(9, 1);

% Encontrar PWM_hover onde empuxo total = peso
% T_total = sum(k_T(i) * T_ref(PWM_h)) = m*g
k_T = P(5:8);
m = constants.m;
g = constants.g;
peso = m * g;

% Resolver para PWM_hover (mesma PWM nos 4 motores)
sum_kT = sum(k_T);
T_ref_hover = peso / sum_kT;  % empuxo de referencia por motor no hover

% Inverter T_ref: encontrar PWM tal que T_ref(PWM) = T_ref_hover
PWM_hover = fzero(@(pw) func_T_ref(pw) - T_ref_hover, 1500);
u0 = PWM_hover * ones(4, 1);

fprintf('==========================================================\n');
fprintf('  LINEARIZACAO DO MODELO VTOL\n');
fprintf('==========================================================\n');
fprintf('  Metodo: Jacobiano numerico (diferencas finitas centrais)\n');
fprintf('  Trim: hover nivelado\n');
fprintf('  PWM hover = %.1f us\n', PWM_hover);
fprintf('  T_ref(PWM_h) = %.4f N/motor\n', func_T_ref(PWM_hover));
fprintf('  T_total = %.4f N  (peso = %.4f N)\n', sum_kT*func_T_ref(PWM_hover), peso);

% Verificar residuo no trim: f(x0, u0) deveria ser ~0 (exceto biases)
f0 = eval_dynamics_direct(x0, u0, P, func_T_ref, func_Q_ref, constants);
fprintf('\n  Residuo f(x0, u0):\n');
state_names = {'p_dot', 'q_dot', 'r_dot', 'phi_dot', 'theta_dot', 'psi_dot', 'u_dot', 'v_dot', 'w_dot'};
for i = 1:9
    fprintf('    %10s = %+.6f\n', state_names{i}, f0(i));
end

%% ========================================================================
%  4. JACOBIANO NUMERICO (DIFERENCAS FINITAS CENTRAIS)
%  ========================================================================
nx = 9;  % numero de estados
nu = 4;  % numero de entradas (PWM1..4)

% Perturbacoes (escala adequada para cada variavel)
h_x = [1e-4; 1e-4; 1e-4;     ... % p, q, r       (rad/s)
       1e-4; 1e-4; 1e-4;     ... % phi, theta, psi (rad)
       1e-3; 1e-3; 1e-3];        % u, v, w        (m/s)

h_u = 0.5 * ones(4, 1);          % PWM (us)

fprintf('\n  Calculando Jacobiano A (9x9)...\n');
A = zeros(nx, nx);
for j = 1:nx
    x_plus  = x0;  x_plus(j)  = x0(j) + h_x(j);
    x_minus = x0;  x_minus(j) = x0(j) - h_x(j);

    f_plus  = eval_dynamics_direct(x_plus,  u0, P, func_T_ref, func_Q_ref, constants);
    f_minus = eval_dynamics_direct(x_minus, u0, P, func_T_ref, func_Q_ref, constants);

    A(:, j) = (f_plus - f_minus) / (2 * h_x(j));
end

fprintf('  Calculando Jacobiano B (9x4)...\n');
B = zeros(nx, nu);
for j = 1:nu
    u_plus  = u0;  u_plus(j)  = u0(j) + h_u(j);
    u_minus = u0;  u_minus(j) = u0(j) - h_u(j);

    f_plus  = eval_dynamics_direct(x0, u_plus,  P, func_T_ref, func_Q_ref, constants);
    f_minus = eval_dynamics_direct(x0, u_minus, P, func_T_ref, func_Q_ref, constants);

    B(:, j) = (f_plus - f_minus) / (2 * h_u(j));
end

%% ========================================================================
%  5. ANALISE DE ESTABILIDADE
%  ========================================================================
eig_A = eig(A);

fprintf('\n  Matriz A (9x9):\n');
disp(A);
fprintf('  Matriz B (9x4):\n');
disp(B);

fprintf('  Autovalores de A:\n');
fprintf('  %-5s  %12s  %12s  %12s\n', '#', 'Real', 'Imag', '|lambda|');
fprintf('  %-5s  %12s  %12s  %12s\n', '---', '----------', '----------', '----------');
for i = 1:length(eig_A)
    fprintf('  %-5d  %+12.6f  %+12.6f  %12.6f\n', ...
        i, real(eig_A(i)), imag(eig_A(i)), abs(eig_A(i)));
end

% Modos naturais
fprintf('\n  Modos naturais:\n');
for i = 1:length(eig_A)
    re = real(eig_A(i));
    im = imag(eig_A(i));
    if abs(im) < 1e-10
        if re < 0
            fprintf('    lambda_%d: modo real estavel (tau = %.3f s)\n', i, -1/re);
        elseif re > 0
            fprintf('    lambda_%d: modo real INSTAVEL (tau = %.3f s)\n', i, 1/re);
        else
            fprintf('    lambda_%d: modo neutro\n', i);
        end
    elseif im > 0
        wn = abs(eig_A(i));
        zeta = -re / wn;
        fprintf('    lambda_%d,%d: modo oscilatorio (wn=%.3f rad/s, zeta=%.3f)\n', ...
            i, i+1, wn, zeta);
    end
end

%% ========================================================================
%  6. SALVAR MODELO LINEAR
%  ========================================================================
save_path = fullfile(fileparts(mfilename('fullpath')), 'linear_model.mat');
save(save_path, 'A', 'B', 'x0', 'u0', 'f0', 'eig_A', 'P', 'PWM_hover');
fprintf('\n  Modelo linear salvo em: %s\n', save_path);
fprintf('==========================================================\n');

%% ========================================================================
%  FUNCAO LOCAL: Avaliar dinamica NL em (x, u) sem interpolacao temporal
%  ========================================================================
function dxdt = eval_dynamics_direct(x, u, P, func_T_ref, func_Q_ref, constants)
    % Cria serie temporal "fake" constante para reutilizar vtol_dynamics.m
    fake_time = [0; 1];
    fake_pwm  = [u(:)'; u(:)'];  % mesma entrada em t=0 e t=1
    dxdt = vtol_dynamics(0.5, x, P, fake_time, fake_pwm, func_T_ref, func_Q_ref, constants);
end
