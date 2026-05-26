function dydt = vtol_dynamics(t, y, P, pwm_time, pwm_signals, func_T_ref, func_Q_ref, constants)
    % Funcao ODE da dinamica do VTOL com dispatch por tamanho de estado.
    %
    % Modos:
    %   length(y) == 3   -> rotacional puro                    [p; q; r]
    %   length(y) == 9   -> completo SEM lag de motor          [p..r; phi..psi; u..w]
    %   length(y) == 17  -> completo COM lag de motor          + [T_eff_1..4; Q_eff_1..4]
    %
    % Modo dispatch: vtol_dynamics('get_handles') retorna struct com handles
    %   .trans_dot      -> derivadas translacionais [ud, vd, wd]
    %   .specific_force -> forca especifica IMU [fx, fy, fz]
    %
    % P (vector de 22 elementos):
    %   P(1:4)   = [Jx, Jy, Jz, Jxz]
    %   P(5:8)   = k_T1..k_T4
    %   P(9:12)  = k_Q1..k_Q4
    %   P(13:15) = Dp, Dq, Dr
    %   P(16:18) = Bp, Bq, Br
    %   P(19:22) = Xu_m, Yv_m, Zw_m, Bz       (opcional, defaults -4,-4,-0.1,-0.5)
    %
    % NOTA: dx_cg, dy_cg removidos (CG oficial = onde dx=dy=0). Braços fixos
    % no CAD: Lx=0.232, Ly_f=0.311185, Ly_r=0.342865.
    %
    % constants: struct com .m, .g; OPCIONAL .tau_motor (escalar ou 4x1).
    %            Se length(y)==17 e tau_motor nao especificado, default 0.05 s.

    %% Dispatch para handles
    if nargin == 1 && ischar(t) && strcmp(t, 'get_handles')
        dydt = struct('trans_dot',      @trans_dot_local, ...
                       'specific_force', @specific_force_local);
        return;
    end

    n_states = length(y);

    %% Parametros (comum aos 3 modos)
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
    k_T = P(5:8);
    k_Q = P(9:12);
    Dp = P(13); Dq = P(14); Dr = P(15);
    Bp = P(16); Bq = P(17); Br = P(18);
    if length(P) >= 22
        Xu_m = P(19); Yv_m = P(20); Zw_m = P(21); Bz_param = P(22);
    else
        Xu_m = -4.0; Yv_m = -4.0; Zw_m = -0.1; Bz_param = -0.5;
    end

    %% Constantes G (corpo rigido)
    gamma0 = Jx*Jz - Jxz^2;
    G1 = Jxz*(Jx - Jy + Jz) / gamma0;
    G2 = (Jz*(Jz - Jy) + Jxz^2) / gamma0;
    G3 = Jz / gamma0;
    G4 = Jxz / gamma0;
    G5 = (Jz - Jx) / Jy;
    G6 = Jxz / Jy;
    G7 = (Jx*(Jx - Jy) + Jxz^2) / gamma0;
    G8 = Jx / gamma0;
    invJy = 1 / Jy;

    %% Braços (do CAD — CG oficial está no ponto onde dx_cg=dy_cg=0)
    Lx_r = 0.232;
    Lx_l = 0.232;
    Ly_f = 0.311185;
    Ly_r = 0.342865;

    p = y(1); q = y(2); r = y(3);

    %% PWM interpolado no instante t
    current_pwm = zeros(1,4);
    for i = 1:4
        current_pwm(i) = interp1(pwm_time, pwm_signals(:,i), t, 'linear', 'extrap');
    end

    %% T_target / Q_target (estaticos, do polinomio com k_T/k_Q)
    Tmr_target = zeros(4,1);
    Qmr_target = zeros(4,1);
    for i = 1:4
        Tmr_target(i) = k_T(i) * func_T_ref(current_pwm(i));
        Qmr_target(i) = k_Q(i) * func_Q_ref(current_pwm(i));
    end

    %% Lag de motor (so se length(y)==17)
    use_lag = (n_states == 17);
    if use_lag
        T_eff = y(10:13);
        Q_eff = y(14:17);

        if nargin >= 8 && ~isempty(constants) && isfield(constants, 'tau_motor')
            tau_motor = constants.tau_motor;
        else
            tau_motor = 0.05;
        end
        if isscalar(tau_motor), tau_motor = tau_motor * ones(4, 1); end

        dT_eff = (Tmr_target - T_eff) ./ tau_motor;
        dQ_eff = (Qmr_target - Q_eff) ./ tau_motor;

        Tmr = T_eff;   % entrada dos momentos = empuxo lagado
        Qmr = Q_eff;
    else
        Tmr = Tmr_target;   % sem lag: motor instantaneo
        Qmr = Qmr_target;
    end

    %% Momentos
    Mx = -(Lx_r*Tmr(1) - Lx_l*Tmr(2) - Lx_l*Tmr(3) + Lx_r*Tmr(4));
    My =   Ly_f*Tmr(1) - Ly_r*Tmr(2) + Ly_f*Tmr(3) - Ly_r*Tmr(4);
    Mz =   Qmr(1) + Qmr(2) - Qmr(3) - Qmr(4);

    %% Aceleracoes angulares
    p_dot = G1*p*q - G2*q*r + G3*Mx + G4*Mz - Dp*p + Bp;
    q_dot = G5*p*r - G6*(p^2 - r^2) + invJy*My - Dq*q + Bq;
    r_dot = G7*p*q - G1*q*r + G4*Mx + G8*Mz - Dr*r + Br;

    if n_states == 3
        dydt = [p_dot; q_dot; r_dot];
        return;
    end

    %% Modo 9 ou 17 estados: cinematica de atitude + translacional
    phi   = y(4);
    theta = y(5);
    % psi   = y(6);   % nao usado nas equacoes (so sai como integral)
    u = y(7); v = y(8); w = y(9);

    if nargin >= 8 && ~isempty(constants)
        m_body = constants.m;
        g_acc  = constants.g;
    else
        m_body = 1.6011;
        g_acc  = 9.81;
    end

    cos_theta = cos(theta);
    if abs(cos_theta) < 1e-7
        cos_theta = 1e-7 * sign(cos_theta + 1e-12);
    end
    sin_theta = sin(theta);
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    tan_theta = sin_theta / cos_theta;

    phi_dot   = p + (q*sin_phi + r*cos_phi) * tan_theta;
    theta_dot = q*cos_phi - r*sin_phi;
    psi_dot   = (q*sin_phi + r*cos_phi) / cos_theta;

    T_total = sum(Tmr);
    gx_body = -g_acc * sin_theta;
    gy_body =  g_acc * sin_phi * cos_theta;
    gz_body =  g_acc * cos_phi * cos_theta;

    [u_dot, v_dot, w_dot] = trans_dot_local(p, q, r, u, v, w, ...
        gx_body, gy_body, gz_body, T_total/m_body, Xu_m, Yv_m, Zw_m, Bz_param);

    if use_lag
        dydt = [p_dot; q_dot; r_dot; ...
                phi_dot; theta_dot; psi_dot; ...
                u_dot; v_dot; w_dot; ...
                dT_eff; dQ_eff];
    else
        dydt = [p_dot; q_dot; r_dot; ...
                phi_dot; theta_dot; psi_dot; ...
                u_dot; v_dot; w_dot];
    end
end

% =========================================================================
%  Subfuncoes translacionais
% =========================================================================
function [ud, vd, wd] = trans_dot_local(p, q, r, u, v, w, gx, gy, gz, T_m, Xu, Yv, Zw, Bz)
    ud = r.*v - q.*w + gx + Xu*u;
    vd = p.*w - r.*u + gy + Yv*v;
    wd = q.*u - p.*v - T_m + gz + Zw*w + Bz;
end

function [fx, fy, fz] = specific_force_local(p, q, r, u, v, w, gx, gy, gz, T_m, Xu, Yv, Zw, Bz)
    fx = r.*v - q.*w + gx + Xu.*u;
    fy = p.*w - r.*u + gy + Yv.*v;
    fz = -T_m + Zw.*w + Bz;
end
