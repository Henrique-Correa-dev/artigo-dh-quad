function dydt = vtol_dynamics(t, y, P, pwm_time, pwm_signals, func_T_ref, func_Q_ref, constants)
    % VTOL_DYNAMICS  Dinâmica do drone — modular via subfunções.
    %
    % DISPATCH ESPECIAL:
    %   dyn_h = vtol_dynamics('get_handles')
    %       Retorna struct com handles pras subfunções:
    %         .moments   = @moments_local(Tmr, Qmr, Lx_r, Lx_l, Ly_f, Ly_r)
    %         .trans_dot = @trans_dot_local(p,q,r, u,v,w, gx,gy,gz, T_m)
    %       Útil pra sim_window, identify_plant: chamam as mesmas fórmulas
    %       sem duplicar (única fonte de verdade).
    %
    % MODOS NORMAIS (dispatch por tamanho de estado):
    %   length(y) == 3   -> rotacional puro             [p; q; r]
    %   length(y) == 9   -> completo SEM lag motor      [p..r; phi..psi; u..w]
    %   length(y) == 17  -> completo COM lag motor      + [T_eff_1..4; Q_eff_1..4]
    %
    % P (vector de 15 elementos — PURAMENTE ROTACIONAL + MOTOR):
    %   P(1:4)   = [Jx, Jy, Jz, Jxz]
    %   P(5:8)   = k_T1..k_T4
    %   P(9:12)  = k_Q1..k_Q4
    %   P(13:15) = Dp, Dq, Dr
    %
    % Parâmetros removidos do modelo (por design):
    %   dx_cg, dy_cg : CG oficial = onde dx=dy=0 (Lx, Ly tratam offset)
    %   Bp, Bq, Br   : CG offset capturado via Lx_r/Lx_l/Ly_f/Ly_r assimétricos
    %   Xu, Yv, Zw   : drag translacional não-identificável sem GPS velocity
    %   Bz           : bias de sensor — vai pra accelerometer_model.m
    %
    % constants: struct com .m, .g; OPCIONAL .tau_motor (escalar ou 4x1).
    %            Se length(y)==17 e tau_motor nao especificado, default 0.05 s.

    %% =========================================================================
    %  Dispatch especial: retorna handles pras subfunções
    %  =========================================================================
    if ischar(t) && strcmp(t, 'get_handles')
        dydt = struct();
        dydt.moments   = @moments_local;
        dydt.trans_dot = @trans_dot_local;
        return;
    end

    n_states = length(y);

    %% =========================================================================
    %  Parâmetros e constantes
    %  =========================================================================
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
    k_T = P(5:8);
    k_Q = P(9:12);
    Dp = P(13); Dq = P(14); Dr = P(15);

    % Constantes G do corpo rígido (Beard-McLain)
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

    % Braços lidos de parameters() — fonte única de verdade
    proj_p = parameters();
    Lx_r = proj_p.arms.Lx_r;
    Lx_l = proj_p.arms.Lx_l;
    Ly_f = proj_p.arms.Ly_f;
    Ly_r = proj_p.arms.Ly_r;

    p = y(1); q = y(2); r = y(3);

    %% =========================================================================
    %  PWM interpolado e empuxo/torque dos motores
    %  =========================================================================
    current_pwm = zeros(1,4);
    for i = 1:4
        current_pwm(i) = interp1(pwm_time, pwm_signals(:,i), t, 'linear', 'extrap');
    end

    % T_target / Q_target estáticos (sem lag) — escalados por k_T, k_Q
    Tmr_target = zeros(4,1);
    Qmr_target = zeros(4,1);
    for i = 1:4
        Tmr_target(i) = k_T(i) * func_T_ref(current_pwm(i));
        Qmr_target(i) = k_Q(i) * func_Q_ref(current_pwm(i));
    end

    % Lag de motor (só se length(y)==17)
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
        Tmr = T_eff;
        Qmr = Q_eff;
    else
        Tmr = Tmr_target;
        Qmr = Qmr_target;
    end

    %% =========================================================================
    %  Momentos no body frame (via subfunção — fonte única)
    %  =========================================================================
    [Mx, My, Mz] = moments_local(Tmr, Qmr, Lx_r, Lx_l, Ly_f, Ly_r);

    %% =========================================================================
    %  Acelerações angulares (rotacional)
    %  =========================================================================
    p_dot = G1*p*q - G2*q*r + G3*Mx + G4*Mz - Dp*p;
    q_dot = G5*p*r - G6*(p^2 - r^2) + invJy*My - Dq*q;
    r_dot = G7*p*q - G1*q*r + G4*Mx + G8*Mz - Dr*r;

    if n_states == 3
        dydt = [p_dot; q_dot; r_dot];
        return;
    end

    %% =========================================================================
    %  Modo 9 ou 17 estados: cinemática de atitude + translacional
    %  =========================================================================
    phi   = y(4);
    theta = y(5);
    % psi   = y(6);   % não usado nas equações (só sai como integral)
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

    %% Cinemática de Euler (taxa de atitude)
    phi_dot   = p + (q*sin_phi + r*cos_phi) * tan_theta;
    theta_dot = q*cos_phi - r*sin_phi;
    psi_dot   = (q*sin_phi + r*cos_phi) / cos_theta;

    %% Gravidade no body frame
    gx_body = -g_acc * sin_theta;
    gy_body =  g_acc * sin_phi * cos_theta;
    gz_body =  g_acc * cos_phi * cos_theta;

    %% Acelerações translacionais (via subfunção — fonte única)
    T_total = sum(Tmr);
    T_m     = T_total / m_body;

    [u_dot, v_dot, w_dot] = trans_dot_local(p, q, r, u, v, w, ...
                                            gx_body, gy_body, gz_body, T_m);

    %% Monta vetor de saída
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

%% =========================================================================
%  SUBFUNÇÃO: Momentos no body frame (ArduPilot QuadX padrão)
%    M1=FR, M2=RL, M3=FL, M4=RR   |   M1+M2 CCW, M3+M4 CW
%  =========================================================================
function [Mx, My, Mz] = moments_local(Tmr, Qmr, Lx_r, Lx_l, Ly_f, Ly_r)
    % Aceita Tmr, Qmr como vetor 4x1, 1x4 ou Nx4 (vetorizado).
    if size(Tmr, 2) == 4   % formato Nx4 (vetorizado)
        Mx = -(Lx_r*Tmr(:,1) - Lx_l*Tmr(:,2) - Lx_l*Tmr(:,3) + Lx_r*Tmr(:,4));
        My =   Ly_f*Tmr(:,1) - Ly_r*Tmr(:,2) + Ly_f*Tmr(:,3) - Ly_r*Tmr(:,4);
        Mz =   Qmr(:,1) + Qmr(:,2) - Qmr(:,3) - Qmr(:,4);
    else                    % formato escalar (4x1 ou 1x4)
        Mx = -(Lx_r*Tmr(1) - Lx_l*Tmr(2) - Lx_l*Tmr(3) + Lx_r*Tmr(4));
        My =   Ly_f*Tmr(1) - Ly_r*Tmr(2) + Ly_f*Tmr(3) - Ly_r*Tmr(4);
        Mz =   Qmr(1) + Qmr(2) - Qmr(3) - Qmr(4);
    end
end

%% =========================================================================
%  SUBFUNÇÃO: Acelerações translacionais (sem drag)
%    Drag (Xu, Yv, Zw) removido — não-identificável sem GPS velocity.
%    Aceita escalares OU vetores Nx1 (broadcast por .*).
%  =========================================================================
function [u_dot, v_dot, w_dot] = trans_dot_local(p, q, r, u, v, w, gx, gy, gz, T_m)
    u_dot = r.*v - q.*w + gx;
    v_dot = p.*w - r.*u + gy;
    w_dot = q.*u - p.*v - T_m + gz;
end
