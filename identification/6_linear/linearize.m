% linearize.m — Linearização do modelo VTOL no HOVER para projeto de controle (LQR/PID)
%
% Estados   x = [p q r phi theta psi u v w]    (9)
% Entradas  v = [T Mx My Mz]                    (4)  ← forças generalizadas
% Saídas    y = x                               (full-state — LQR usa todos)
%
% Modelo linear:  dx/dt = A·δx + B·δv     (δx = x - x0, δv = v - v0)
%   A = ∂f/∂x |_(x0,v0)   (9×9)   — Jacobiano numérico (diferenças centrais)
%   B = ∂f/∂v |_(x0,v0)   (9×4)   — CONSTANTE (corpo rígido linear em [T,M])
%
% POR QUE [T,M] (e não PWM): a não-linearidade dos motores (curvas fT/fQ) fica
% isolada na ALOCAÇÃO (mixer), fora da dinâmica. Assim B é geometria/inércia
% pura e vale no envelope inteiro (não só num PWM de trim). Ganhos físicos:
% transfer functions tipo θ̈/My = (1/Jy)/s, ẅ/T = -(1/m)/s.
%
% TRIM: no espaço de forças, o trim de hover é TRIVIAL: v0 = [m·g; 0; 0; 0]
%       (empuxo = peso, momentos nulos) → c = f(x0,v0) ≈ 0 exato.
%       O PWM de trim (pwm_trim) ainda é resolvido por fsolve, pois é necessário
%       pra ALOCAÇÃO (mixer) e pra realismo — mas NÃO é a entrada do modelo linear.
%
% ALOCAÇÃO (mixer): M_alloc = ∂[T,Mx,My,Mz]/∂PWM |_trim  (4×4).
%   Forward (sim): PWM → [T,M] via M_alloc (localmente) ou forces (exato).
%   Inverse (controle): [T,M] desejado → PWM via M_mix = inv(M_alloc).
%   Construído a partir das SUAS tabelas fT, fQ (inclui a linha de yaw via fQ).
%
% Saída salva em outputs/linear_model.mat:
%   A, B, C, D, x0, u0(=v0=[mg;0;0;0]), sys, eig_A, P, pwm_trim,
%   M_alloc, M_mix, trim_residual, rank_Co
%
% Uso:  >> linearize

%% ========================================================================
%  0. SETUP
%  ========================================================================
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ========================================================================
%  1. PARÂMETROS (P_J 15 elementos)
%  ========================================================================
proj = parameters();

% Fonte do P: usa identificado se existir, senão P0
P = proj.P0_J;
P_file = fullfile(paths.outputs, 'P_identified.mat');
if exist(P_file, 'file')
    id = load(P_file);
    P = id.P_final;
    fprintf('linearize: usando P_final de P_identified.mat\n');
else
    fprintf('linearize: usando P0 (chute inicial de parameters.m)\n');
end

constants.m = proj.m;
constants.g = proj.g;

[func_T_ref, func_Q_ref] = motor_models();

% Lag de motor 1ª ordem 1/(τs+1) nas forças generalizadas [T,Mx,My,Mz].
%  Como os 4 motores compartilham τ e [T,M] são combinações LINEARES das
%  trações, o lag aparece igual nas 4 entradas → 4 estados de atuador em série
%  (modelo aumentado 13×13). Use >0 p/ projeto de controle (a dinâmica do motor
%  reduz a margem de fase); 0 = só a planta rígida 9×9 (default, não altera nada).
LAG_TAU = 0;     % s   (ex.: 0.05 físico ESC/motor)

%% ========================================================================
%  2. PONTO DE OPERAÇÃO (hover)
%     - Estado: x0 = 0 (nivelado, parado)
%     - Força:  v0 = [m·g; 0; 0; 0]  (trim trivial no espaço de forças)
%     - PWM de trim por fsolve (só pra alocação/realismo)
%  ========================================================================
x0 = zeros(9, 1);
peso = constants.m * constants.g;
v0 = [peso; 0; 0; 0];

trim_fun = @(pwm) trim_residuals(pwm, P, func_T_ref, func_Q_ref, proj.arms, peso);
opts_fs = optimoptions('fsolve', 'Display', 'off', ...
                       'FunctionTolerance', 1e-10, 'StepTolerance', 1e-12);
pwm_chute = 1500 * ones(1,4);
[pwm_trim, res_trim, exitflag] = fsolve(trim_fun, pwm_chute, opts_fs);
pwm_trim = pwm_trim(:);

fprintf('==========================================================\n');
fprintf('  LINEARIZAÇÃO VTOL NO HOVER — entrada [T, Mx, My, Mz]\n');
fprintf('==========================================================\n');
fprintf('  Trim de força v0 = [T=%.2f N, Mx=0, My=0, Mz=0]  (T=m·g)\n', peso);
fprintf('  PWM trim (p/ alocação) fsolve ef=%d |res|=%.2e: [%.1f %.1f %.1f %.1f] us\n', ...
        exitflag, norm(res_trim), pwm_trim);

%% ========================================================================
%  3. LINEARIZAÇÃO em torno de (x0, v0) (fonte única: linearize_at.m)
%     A = ∂f/∂x, B = ∂f/∂v, c = f(x0,v0) (deve ser ~0 no equilíbrio)
%  ========================================================================
nx = 9;  nv = 4;
fprintf('\n  Linearizando o corpo rígido em [T,M] (linearize_at)...\n');
[A, B, c] = linearize_at(x0, v0, P, constants);

state_names = {'p_dot','q_dot','r_dot','phi_dot','theta_dot','psi_dot', ...
               'u_dot','v_dot','w_dot'};
fprintf('  Resíduo c = f(x0,v0) (deve ser ~0):  |c| = %.2e\n', norm(c));

%% ========================================================================
%  4. ALOCAÇÃO (mixer): M_alloc = ∂[T,Mx,My,Mz]/∂PWM |_trim   (4×4)
%     Construída das tabelas fT, fQ. M_mix = inv(M_alloc) (controle → PWM).
%  ========================================================================
dyn = vtol_dynamics('get_handles');
forces = dyn.forces;
h_pwm = 0.5;   % us
M_alloc = zeros(4, 4);
for j = 1:4
    pp = pwm_trim'; pp(j) = pp(j) + h_pwm;
    pm = pwm_trim'; pm(j) = pm(j) - h_pwm;
    [Tp,Mxp,Myp,Mzp] = forces(pp, P, func_T_ref, func_Q_ref);
    [Tm,Mxm,Mym,Mzm] = forces(pm, P, func_T_ref, func_Q_ref);
    M_alloc(:,j) = ([Tp;Mxp;Myp;Mzp] - [Tm;Mxm;Mym;Mzm]) / (2*h_pwm);
end
M_mix = inv(M_alloc);
fprintf('\n  --- Alocação (mixer) ---\n');
fprintf('    cond(M_alloc) = %.1f', cond(M_alloc));
if cond(M_alloc) > 1e3
    fprintf('   (mal-condicionado — eixo de yaw fraco: pouco Mz por muito PWM)\n');
else
    fprintf('\n');
end

%% ========================================================================
%  5. SAÍDAS (full-state feedback — LQR mede todos os estados)
%  ========================================================================
C = eye(nx);
D = zeros(nx, nv);
u0 = v0;   % entrada de trim do modelo linear = força de hover
state_lbl = {'p','q','r','phi','theta','psi','u','v','w'};
input_lbl = {'T','Mx','My','Mz'};
sys = ss(A, B, C, D, ...
         'StateName', state_lbl, 'InputName', input_lbl, 'OutputName', state_lbl);

%% ========================================================================
%  5b. MODELO AUMENTADO COM LAG DE MOTOR (opcional, p/ controle)
%      x_aug = [x(9); v_lag(4)];  v_lag = 1/(τs+1)·v_cmd (lag comum nos 4 motores)
%      dv_lag/dt = (v_cmd - v_lag)/τ ;  a planta vê v_lag (não v_cmd).
%  ========================================================================
if LAG_TAU > 0
    Aact  = -(1/LAG_TAU) * eye(nv);
    Bact  =  (1/LAG_TAU) * eye(nv);
    A_lag = [A,             B;
             zeros(nv, nx), Aact];
    B_lag = [zeros(nx, nv);
             Bact];
    C_lag = [C, zeros(nx, nv)];     % saídas = mesmos 9 estados físicos
    D_lag = zeros(nx, nv);
    state_lbl_lag = [state_lbl, {'T_lag','Mx_lag','My_lag','Mz_lag'}];
    sys_lag = ss(A_lag, B_lag, C_lag, D_lag, ...
        'StateName', state_lbl_lag, 'InputName', input_lbl, 'OutputName', state_lbl);
    eig_A_lag = eig(A_lag);
    rank_Co_lag = rank(ctrb(A_lag, B_lag));
    fprintf('\n  --- Modelo AUMENTADO com lag de motor (τ=%.3f s) ---\n', LAG_TAU);
    fprintf('    %d estados (9 planta + 4 atuador) | %d polos de motor em -1/τ=%.1f\n', ...
        nx+nv, nv, -1/LAG_TAU);
    fprintf('    rank(ctrb_aug) = %d / %d\n', rank_Co_lag, nx+nv);
end

%% ========================================================================
%  6. ANÁLISE DE CONTROLABILIDADE (essencial pra LQR)
%  ========================================================================
Co = ctrb(A, B);
rank_Co = rank(Co);
fprintf('\n  --- Controlabilidade ---\n');
fprintf('    rank(ctrb) = %d / %d', rank_Co, nx);
if rank_Co == nx
    fprintf('  → SISTEMA CONTROLÁVEL\n');
else
    fprintf('  → NÃO totalmente controlável (%d modos não-controláveis)\n', nx - rank_Co);
    fprintf('    (esperado: psi e posição horizontal costumam ser marginais)\n');
end

%% ========================================================================
%  7. ANÁLISE MODAL (autovalores / polos da planta)
%  ========================================================================
eig_A = eig(A);
fprintf('\n  --- Autovalores de A (polos da planta) ---\n');
fprintf('  %-5s  %12s  %12s  %12s\n', '#', 'Real', 'Imag', '|lambda|');
for i = 1:nx
    fprintf('  %-5d  %+12.5f  %+12.5f  %12.5f\n', ...
        i, real(eig_A(i)), imag(eig_A(i)), abs(eig_A(i)));
end

fprintf('\n  --- Modos naturais ---\n');
done = false(nx,1);
for i = 1:nx
    if done(i), continue; end
    re = real(eig_A(i)); im = imag(eig_A(i));
    if abs(im) < 1e-9
        if re < -1e-9
            fprintf('    lambda=%+.4f: real estavel   (tau=%.3f s)\n', re, -1/re);
        elseif re > 1e-9
            fprintf('    lambda=%+.4f: real INSTAVEL  (tau=%.3f s)\n', re, 1/re);
        else
            fprintf('    lambda=%+.4f: integrador/neutro\n', re);
        end
    else
        wn = abs(eig_A(i)); zeta = -re/wn;
        fprintf('    lambda=%+.4f+-%.4fj: oscilatorio (wn=%.3f rad/s, zeta=%.3f)\n', ...
                re, abs(im), wn, zeta);
        done(i) = true;  % marca o par conjugado
    end
end

%% ========================================================================
%  8. SALVAR
%  ========================================================================
trim_residual = c;
save_path = fullfile(paths.outputs, 'linear_model.mat');
save(save_path, 'A', 'B', 'C', 'D', 'x0', 'u0', 'sys', ...
     'eig_A', 'P', 'pwm_trim', 'M_alloc', 'M_mix', 'trim_residual', 'rank_Co');
if LAG_TAU > 0
    save(save_path, 'A_lag', 'B_lag', 'C_lag', 'D_lag', 'sys_lag', ...
         'eig_A_lag', 'LAG_TAU', '-append');
    fprintf('  + modelo aumentado c/ lag: A_lag B_lag C_lag D_lag sys_lag eig_A_lag LAG_TAU\n');
end
fprintf('\n  Modelo linear salvo: %s\n', save_path);
fprintf('  Entrada: [T, Mx, My, Mz] | u0=[mg,0,0,0] | + alocação M_alloc/M_mix\n');
fprintf('  Vars: A B C D x0 u0 sys eig_A P pwm_trim M_alloc M_mix trim_residual rank_Co\n');
fprintf('==========================================================\n');


%% ========================================================================
%  FUNÇÕES LOCAIS
%  ========================================================================
function r = trim_residuals(pwm, P, func_T_ref, func_Q_ref, arms, peso)
%TRIM_RESIDUALS  Resíduo de equilíbrio em hover: [Mx; My; Mz; T_total - peso].
%   Zera quando os 4 PWM produzem momento nulo e empuxo = peso.
    k_T = P(5:8);  k_Q = P(9:12);
    T = zeros(4,1); Q = zeros(4,1);
    for i = 1:4
        T(i) = k_T(i) * func_T_ref(pwm(i));
        Q(i) = k_Q(i) * func_Q_ref(pwm(i));
    end
    % Momentos (ArduPilot QuadX — mesma fórmula de vtol_dynamics.m)
    Mx = -(arms.Lx_r*T(1) - arms.Lx_l*T(2) - arms.Lx_l*T(3) + arms.Lx_r*T(4));
    My =   arms.Ly_f*T(1) - arms.Ly_r*T(2) + arms.Ly_f*T(3) - arms.Ly_r*T(4);
    Mz =   Q(1) + Q(2) - Q(3) - Q(4);
    r = [Mx; My; Mz; sum(T) - peso];
end
