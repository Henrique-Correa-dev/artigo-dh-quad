function x_dot = vtol_dynamics(t, x, t_data, pwm_data, params)
% VTOL_DYNAMICS  Modelo dinamico completo do DH Hybrid-Drone (modo quadrotor)
%
% Equacoes de movimento de corpo rigido com 9 estados, baseadas no
% Beard & McLain "Small Unmanned Aircraft", Cap. 3 (eq. 3.13-3.16).
%
% Estados (9):
%   x = [p; q; r; phi; theta; psi; u; v; w]
%
%   p, q, r       - velocidades angulares no corpo [rad/s]
%   phi, theta, psi - angulos de Euler [rad]
%   u, v, w       - velocidades lineares no corpo [m/s]
%
% Entradas:
%   t        - tempo [s]
%   x        - vetor de estados 9x1
%   t_data   - vetor de tempo dos dados interpolados
%   pwm_data - matriz Nx4 de PWM [us] (colunas: C1, C2, C3, C4)
%   params   - struct com parametros do modelo
%
% Saida:
%   x_dot    - derivada dos estados 9x1
%
% =========================================================================
%
%  DIAGRAMA DO MODELO:
%
%  PWM(t)──► T_ref(pwm) ──► Mixing ──► [T, Mx, My, Mz]
%  [C1..C4]  Q_ref(pwm)     Matrix       │
%                                          │
%       ┌──────────────────────────────────┘
%       │
%       ▼
%  ┌─────────────────────────────────────────────────────────┐
%  │  DINAMICA ROTACIONAL (eq. 3.11)                         │
%  │                                                         │
%  │  ω_dot = J⁻¹ · ( -ω × (J·ω) + m )                     │
%  │                                                         │
%  │  onde ω = [p;q;r],  m = [Mx;My;Mz]                     │
%  │  J = tensor de inercia completo (com Jxy, Jxz, Jyz)    │
%  │                                                         │
%  │  Saida: p_dot, q_dot, r_dot                             │
%  └─────────────────────────────────────────────────────────┘
%       │
%       ▼
%  ┌─────────────────────────────────────────────────────────┐
%  │  CINEMATICA EULER (eq. 3.3)                             │
%  │                                                         │
%  │  phi_dot   = p + sin(phi)*tan(theta)*q                  │
%  │                + cos(phi)*tan(theta)*r                   │
%  │  theta_dot = cos(phi)*q - sin(phi)*r                    │
%  │  psi_dot   = sin(phi)/cos(theta)*q                      │
%  │                + cos(phi)/cos(theta)*r                   │
%  │                                                         │
%  │  Saida: phi_dot, theta_dot, psi_dot                     │
%  └─────────────────────────────────────────────────────────┘
%       │
%       ▼
%  ┌─────────────────────────────────────────────────────────┐
%  │  DINAMICA TRANSLACIONAL (eq. 3.7)                       │
%  │                                                         │
%  │  u_dot = r*v - q*w + (1/m) * fx                        │
%  │  v_dot = p*w - r*u + (1/m) * fy                        │
%  │  w_dot = q*u - p*v + (1/m) * fz                        │
%  │                                                         │
%  │  fx,fy,fz = gravidade + empuxo + arrasto                │
%  │                                                         │
%  │  Saida: u_dot, v_dot, w_dot                             │
%  └─────────────────────────────────────────────────────────┘
%
% =========================================================================

    %% Extrair estados
    p     = x(1);
    q     = x(2);
    r     = x(3);
    phi   = x(4);
    theta = x(5);
    % psi = x(6);   % nao aparece nas equacoes dinamicas
    u     = x(7);
    v     = x(8);
    w     = x(9);

    %% Extrair parametros
    m_kg  = params.mass;        % massa [kg]
    g     = params.g;           % gravidade [m/s^2]
    J     = params.J;           % tensor de inercia 3x3 [kg.m^2]
    J_inv = params.J_inv;       % inversa do tensor de inercia

    % Geometria
    x_m   = params.x_m;        % posicao longitudinal dos motores [m] (4x1)
    y_m   = params.y_m;        % posicao lateral dos motores [m] (4x1)
    d_m   = params.d_m;        % direcao yaw: +1=CW, -1=CCW (4x1)

    % Fatores de escala individuais
    k_T   = params.k_T;        % fator empuxo por motor (4x1)
    k_Q   = params.k_Q;        % fator torque por motor (4x1)

    % Amortecimento rotacional
    Dp    = params.Dp;          % amortecimento roll [1/s]
    Dq    = params.Dq;          % amortecimento pitch [1/s]
    Dr    = params.Dr;          % amortecimento yaw [1/s]

    % Bias rotacional
    Bp    = params.Bp;          % bias roll [rad/s^2]
    Bq    = params.Bq;          % bias pitch [rad/s^2]
    Br    = params.Br;          % bias yaw [rad/s^2]

    % Arrasto translacional (derivadas de estabilidade)
    Xu    = params.Xu;          % arrasto em x [1/s]
    Yv    = params.Yv;          % arrasto em y [1/s]
    Zw    = params.Zw;          % arrasto em z [1/s]
    Bz    = params.Bz;          % bias vertical [m/s^2]

    % Modelos de bancada (function handles)
    T_ref = params.T_ref;       % PWM -> empuxo [N]
    Q_ref = params.Q_ref;       % PWM -> torque [Nm]

    %% Interpolar PWM no instante t
    pwm = zeros(4,1);
    for i = 1:4
        pwm(i) = interp1(t_data, pwm_data(:,i), t, 'linear', 'extrap');
    end

    %% Calcular forcas e momentos dos motores
    % Empuxo e torque de cada motor (com fator de escala individual)
    Ti = zeros(4,1);
    Qi = zeros(4,1);
    for i = 1:4
        Ti(i) = k_T(i) * T_ref(pwm(i));
        Qi(i) = k_Q(i) * Q_ref(pwm(i));
    end

    % Empuxo total
    T_total = sum(Ti);

    % Momentos via mixing matrix (tau = r x F + torque reativo)
    %   Mx (roll)  = sum( -yi * Ti )
    %   My (pitch) = sum( +xi * Ti )
    %   Mz (yaw)   = sum( di * Qi )
    Mx = sum(-y_m .* Ti);
    My = sum(+x_m .* Ti);
    Mz = sum( d_m .* Qi);

    % Adicionar amortecimento e bias
    Mx = Mx - Dp*p + Bp;
    My = My - Dq*q + Bq;
    Mz = Mz - Dr*r + Br;

    moment = [Mx; My; Mz];

    %% ====================================================================
    %  DINAMICA ROTACIONAL (eq. 3.11)
    %  omega_dot = J_inv * ( -omega x (J*omega) + moment )
    % =====================================================================
    omega = [p; q; r];

    omega_dot = J_inv * (-cross(omega, J * omega) + moment);

    p_dot = omega_dot(1);
    q_dot = omega_dot(2);
    r_dot = omega_dot(3);

    %% ====================================================================
    %  CINEMATICA DE EULER (eq. 3.3)
    %  Relaciona derivadas dos angulos com velocidades angulares
    % =====================================================================
    sp = sin(phi);   cp = cos(phi);
    st = sin(theta); ct = cos(theta); tt = tan(theta);

    phi_dot   = p + sp*tt*q + cp*tt*r;
    theta_dot = cp*q - sp*r;
    psi_dot   = (sp/ct)*q + (cp/ct)*r;

    %% ====================================================================
    %  DINAMICA TRANSLACIONAL (eq. 3.7)
    %  Forcas no frame do corpo: gravidade + empuxo + arrasto
    % =====================================================================

    % Gravidade projetada no corpo (eq. na secao 4.1)
    fg_x = -m_kg * g * st;
    fg_y =  m_kg * g * sp * ct;
    fg_z =  m_kg * g * cp * ct;

    % Empuxo no corpo: F = [0; 0; -T] (aponta para cima, -k_b)
    ft_x = 0;
    ft_y = 0;
    ft_z = -T_total;

    % Arrasto aerodinamico translacional (linear em velocidade)
    fd_x = Xu * u * m_kg;
    fd_y = Yv * v * m_kg;
    fd_z = (Zw * w + Bz) * m_kg;

    % Forca total no corpo
    fx = fg_x + ft_x + fd_x;
    fy = fg_y + ft_y + fd_y;
    fz = fg_z + ft_z + fd_z;

    % Equacoes de Newton (eq. 3.7)
    u_dot = r*v - q*w + fx/m_kg;
    v_dot = p*w - r*u + fy/m_kg;
    w_dot = q*u - p*v + fz/m_kg;

    %% Montar vetor de derivadas
    x_dot = [p_dot; q_dot; r_dot; ...
             phi_dot; theta_dot; psi_dot; ...
             u_dot; v_dot; w_dot];
end
