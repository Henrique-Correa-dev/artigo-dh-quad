function dydt = vtol_dynamics(t, y, P, pwm_time, pwm_signals, func_T_ref, func_Q_ref, constants)
    % Função ODE da dinâmica do VTOL
    %
    % Modo 3 estados: y = [p; q; r]                              -> rotacional
    % Modo 9 estados: y = [p; q; r; phi; theta; psi; u; v; w]    -> completa (acoplado)
    %
    % P: vetor de 20 ou 24 parâmetros identificáveis
    %   P(1:4)   = [Jx, Jy, Jz, Jxz] (momentos de inércia)
    %   P(5:8)   = k_T1..k_T4
    %   P(9:12)  = k_Q1..k_Q4
    %   P(13:15) = Dp, Dq, Dr
    %   P(16:18) = Bp, Bq, Br
    %   P(19:20) = [dx_cg, dy_cg] (offset CG vs CAD)
    %   P(21:24) = [Xu_m, Yv_m, Zw_m, Bz] (opcional, defaults se ausente)
    % constants: struct com .m, .g

    n_states = length(y);

    % Extrair CG offsets
    if length(P) >= 20
        dx_cg = P(19); dy_cg = P(20);
    else
        dx_cg = 0; dy_cg = 0;
    end

    % Extrair parâmetros translacionais (defaults do Simulink se não fornecidos)
    if length(P) >= 24
        Xu_m = P(21); Yv_m = P(22); Zw_m = P(23); Bz_param = P(24);
    else
        Xu_m = -4.0; Yv_m = -4.0; Zw_m = -0.1; Bz_param = -0.5;
    end

    p = y(1);
    q = y(2);
    r = y(3);

    % Desempacotar inércias e computar constantes G (corpo rígido)
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
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

    k_T = P(5:8);
    k_Q = P(9:12);

    Dp = P(13);
    Dq = P(14);
    Dr = P(15);

    Bp = P(16);
    Bq = P(17);
    Br = P(18);

    % Interpolar sinais PWM no tempo atual t
    current_pwm = zeros(1,4);
    for i = 1:4
        current_pwm(i) = interp1(pwm_time, pwm_signals(:,i), t, 'linear', 'extrap');
    end
    PWM1 = current_pwm(1);
    PWM2 = current_pwm(2);
    PWM3 = current_pwm(3);
    PWM4 = current_pwm(4);

    % Empuxo individual (N) e torques reativos (N·m)
    Tmr1 = k_T(1) * func_T_ref(PWM1);
    Tmr2 = k_T(2) * func_T_ref(PWM2);
    Tmr3 = k_T(3) * func_T_ref(PWM3);
    Tmr4 = k_T(4) * func_T_ref(PWM4);

    Qmr1 = k_Q(1) * func_Q_ref(PWM1);
    Qmr2 = k_Q(2) * func_Q_ref(PWM2);
    Qmr3 = k_Q(3) * func_Q_ref(PWM3);
    Qmr4 = k_Q(4) * func_Q_ref(PWM4);

    % Braços efetivos com offset do CG
    Lx_r = 0.232 - dy_cg;   % direita (motores 1,4)
    Lx_l = 0.232 + dy_cg;   % esquerda (motores 2,3)
    Ly_f = 0.311185 - dx_cg; % frente (motores 1,3)
    Ly_r = 0.342865 + dx_cg; % traseira (motores 2,4)

    % Momentos Mx, My, Mz
    Mx = -(Lx_r*Tmr1 - Lx_l*Tmr2 - Lx_l*Tmr3 + Lx_r*Tmr4);
    My = (Ly_f*Tmr1 - Ly_r*Tmr2 + Ly_f*Tmr3 - Ly_r*Tmr4);
    Mz = Qmr1 + Qmr2 - Qmr3 - Qmr4;

    % Acelerações angulares
    p_dot = G1*p*q - G2*q*r + G3*Mx + G4*Mz - Dp*p + Bp;
    q_dot = G5*p*r - G6*(p^2 - r^2) + invJy*My - Dq*q + Bq;
    r_dot = G7*p*q - G1*q*r + G4*Mx + G8*Mz - Dr*r + Br;

    if n_states == 3
        dydt = [p_dot; q_dot; r_dot];
        return;
    end

    % =====================================================================
    %  Modo 9 estados: cinemática de atitude + dinâmica translacional
    % =====================================================================
    phi   = y(4);
    theta = y(5);
    psi   = y(6);
    u = y(7);
    v = y(8);
    w = y(9);

    if nargin >= 8 && ~isempty(constants)
        m_body = constants.m;
        g_acc  = constants.g;
    else
        m_body = 1.6011;
        g_acc  = 9.81;
    end

    [phi_dot, theta_dot, psi_dot] = euler_kinematics(p, q, r, phi, theta);
    [u_dot, v_dot, w_dot] = translational_eqs(p, q, r, u, v, w, phi, theta, psi, ...
        Tmr1+Tmr2+Tmr3+Tmr4, m_body, g_acc, Xu_m, Yv_m, Zw_m, Bz_param);

    dydt = [p_dot; q_dot; r_dot; phi_dot; theta_dot; psi_dot; u_dot; v_dot; w_dot];
end

% =========================================================================
%  Equações comuns: cinemática de Euler
% =========================================================================
function [phi_dot, theta_dot, psi_dot] = euler_kinematics(p, q, r, phi, theta)
    cos_theta = cos(theta);
    if abs(cos_theta) < 1e-7
        cos_theta = 1e-7 * sign(cos_theta);
    end
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    tan_theta = sin(theta) / cos_theta;

    phi_dot   = p + (q*sin_phi + r*cos_phi) * tan_theta;
    theta_dot = q*cos_phi - r*sin_phi;
    psi_dot   = (q*sin_phi + r*cos_phi) / cos_theta;
end

% =========================================================================
%  Equações comuns: dinâmica translacional (com arrasto e bias do Simulink)
% =========================================================================
function [u_dot, v_dot, w_dot] = translational_eqs(p, q, r, u, v, w, phi, theta, psi, T_total, m, g, Xu_m, Yv_m, Zw_m, Bz_param)
    sin_phi   = sin(phi);
    cos_phi   = cos(phi);
    cos_theta = cos(theta);

    R_nb = [ cos(theta)*cos(psi), sin_phi*sin(theta)*cos(psi)-cos_phi*sin(psi), cos_phi*sin(theta)*cos(psi)+sin_phi*sin(psi);
             cos(theta)*sin(psi), sin_phi*sin(theta)*sin(psi)+cos_phi*cos(psi), cos_phi*sin(theta)*sin(psi)-sin_phi*cos(psi);
            -sin(theta),          sin_phi*cos_theta,                            cos_phi*cos_theta];
    G_body = R_nb' * [0; 0; m*g];

    Fx = G_body(1);
    Fy = G_body(2);
    Fz = -T_total + G_body(3);

    u_dot = r*v - q*w + Fx/m + Xu_m*u;
    v_dot = p*w - r*u + Fy/m + Yv_m*v;
    w_dot = q*u - p*v + Fz/m + Zw_m*w + Bz_param;
end
