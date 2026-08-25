function bridge = forces_to_pwm(lm, mode)
%FORCES_TO_PWM  Constrói a PONTE (control allocation) u=[T,Mx,My,Mz] → PWM.
%
% É o elo entre o CONTROLADOR (que comanda forças/momentos [T,M]) e a PLANTA
% NÃO-LINEAR (que recebe PWM). Retorna um handle  pwm = bridge(u_cmd).
%
% Dois modos (correspondem aos blocos do diagrama do modelo não-linear):
%
%   'linear'    (padrão p/ robustez):
%       PWM = pwm_trim + M_mix·(u_cmd - u0)
%       M_mix = inv(M_alloc) — inverso EXATO da linearização do mixer no hover.
%       A não-linearidade fT⁻¹ já está embutida na inclinação de M_alloc.
%       Garante casamento NL ↔ linear perto do hover.
%
%   'nonlinear' (fiel ao diagrama: M_mix + fT⁻¹ em DOIS passos):
%       1) Alocação em ESPAÇO DE EMPUXO:  [T1..T4] = G_t \ [T,Mx,My,Mz]
%          (G_t monta T,Mx,My exatos + Mz≈Σ±κ_i·T_i, κ_i = Q_i/T_i no trim)
%       2) Motor inverso fT⁻¹:  PWM_i = fT_i⁻¹(T_i)  (inversão da spline Akima)
%
% Em ambos: PWM saturado em [1000, 2000] µs.
%
% ENTRADAS
%   lm    struct do linear_model.mat (usa M_mix, M_alloc, pwm_trim, u0, P)
%   mode  'linear' (padrão) | 'nonlinear'
%
% SAÍDA
%   bridge  function handle:  pwm(4x1) = bridge(u_cmd(4x1))

    if nargin < 2 || isempty(mode), mode = 'linear'; end
    PWM_MIN = 1000;  PWM_MAX = 2000;

    pwm_trim = lm.pwm_trim(:);
    u0       = lm.u0(:);

    switch lower(mode)
        case 'linear'
            M_mix = lm.M_mix;
            bridge = @(u_cmd) clamp_pwm(pwm_trim + M_mix*(u_cmd(:) - u0), PWM_MIN, PWM_MAX);

        case 'nonlinear'
            proj = parameters();
            P = lm.P;  arms = proj.arms;
            k_T = P(5:8);  k_Q = P(9:12);
            [fT, fQ] = motor_models();

            % κ_i = Q_i/T_i avaliado no PWM de trim (razão contra-torque/empuxo)
            kappa = zeros(4,1);
            for i = 1:4
                Ti = k_T(i)*fT(pwm_trim(i));
                Qi = k_Q(i)*fQ(pwm_trim(i));
                kappa(i) = Qi / max(Ti, eps);
            end

            % G_t: [T;Mx;My;Mz] = G_t · [T1;T2;T3;T4]   (mesma convenção ArduPilot QuadX)
            G_t = [ 1,           1,           1,           1; ...
                   -arms.Lx_r,   arms.Lx_l,   arms.Lx_l,  -arms.Lx_r; ...
                    arms.Ly_f,  -arms.Ly_r,   arms.Ly_f,  -arms.Ly_r; ...
                    kappa(1),    kappa(2),   -kappa(3),   -kappa(4) ];

            % Tabela inversa fT⁻¹ por motor (spline é monotônica em [1000,2000])
            pwm_grid = (PWM_MIN:1:PWM_MAX)';
            T_grid   = fT(pwm_grid);                 % empuxo de 1 motor de referência (N)

            bridge = @(u_cmd) alloc_nl(u_cmd(:), G_t, k_T, T_grid, pwm_grid, PWM_MIN, PWM_MAX);

        otherwise
            error('forces_to_pwm:badMode', 'mode deve ser "linear" ou "nonlinear".');
    end
end


function pwm = alloc_nl(u_cmd, G_t, k_T, T_grid, pwm_grid, lo, hi)
    T_motor = G_t \ u_cmd;                  % empuxo desejado por motor (N)
    pwm = zeros(4,1);
    for i = 1:4
        T_ref = max(T_motor(i), 0) / k_T(i); % empuxo equivalente do motor de referência
        pwm(i) = interp1(T_grid, pwm_grid, T_ref, 'linear', 'extrap');
    end
    pwm = clamp_pwm(pwm, lo, hi);
end

function p = clamp_pwm(p, lo, hi)
    p = min(max(p, lo), hi);
end
