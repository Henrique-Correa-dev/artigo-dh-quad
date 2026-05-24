function e = eem_cost(P, params, p, q, r, p_dot, q_dot, r_dot, T_ref_mat, Q_ref_mat, weights)
% EEM_COST  Equation Error Method cost for rotational identification.
%
% Identifica inercias + fatores de escala + amortecimento + bias:
%   P(1:6)   = [Jx, Jy, Jz, Jxy, Jxz, Jyz]  (momentos de inercia)
%   P(7:10)  = k_T1..k_T4   (fatores de escala de empuxo)
%   P(11:14) = k_Q1..k_Q4   (fatores de escala de torque)
%   P(15:17) = Dp, Dq, Dr   (amortecimento rotacional)
%   P(18:20) = Bp, Bq, Br   (bias rotacional)
%
% Residuo: omega_dot_medido - omega_dot_modelo
%   onde omega_dot_modelo = J_inv * (-omega x (J*omega) + M)
%
% Entradas:
%   P          - vetor 20x1 de parametros
%   params     - struct com x_m, y_m, d_m (geometria)
%   p,q,r      - velocidades angulares medidas (Nx1)
%   p_dot,...   - derivadas numericas (Nx1)
%   T_ref_mat  - empuxo de referencia Nx4
%   Q_ref_mat  - torque de referencia Nx4
%   weights    - pesos [w_p; w_q; w_r] (3x1)

    N = length(p);

    % Inercias do vetor de parametros
    Jx = P(1); Jy = P(2); Jz = P(3);
    Jxy = P(4); Jxz = P(5); Jyz = P(6);

    J = [Jx,  Jxy, Jxz;
         Jxy, Jy,  Jyz;
         Jxz, Jyz, Jz ];
    J_inv = inv(J);

    k_T = P(7:10);
    k_Q = P(11:14);
    Dp = P(15); Dq = P(16); Dr = P(17);
    Bp = P(18); Bq = P(19); Br = P(20);

    x_m = params.x_m;
    y_m = params.y_m;
    d_m = params.d_m;

    % Empuxo e torque escalados (Nx4)
    Ti = T_ref_mat .* k_T';
    Qi = Q_ref_mat .* k_Q';

    % Momentos via mixing matrix
    Mx = sum(-y_m' .* Ti, 2);   % roll
    My = sum(+x_m' .* Ti, 2);   % pitch
    Mz = sum( d_m' .* Qi, 2);   % yaw

    % Adicionar amortecimento e bias
    Mx = Mx - Dp*p + Bp;
    My = My - Dq*q + Bq;
    Mz = Mz - Dr*r + Br;

    % Modelo: omega_dot = J_inv * (-omega x (J*omega) + M)
    p_dot_model = zeros(N,1);
    q_dot_model = zeros(N,1);
    r_dot_model = zeros(N,1);

    for k = 1:N
        omega_k = [p(k); q(k); r(k)];
        M_k     = [Mx(k); My(k); Mz(k)];
        gyro    = cross(omega_k, J * omega_k);
        omegad  = J_inv * (-gyro + M_k);

        p_dot_model(k) = omegad(1);
        q_dot_model(k) = omegad(2);
        r_dot_model(k) = omegad(3);
    end

    sqrt_w = sqrt(weights(:));
    e = [sqrt_w(1) * (p_dot - p_dot_model); ...
         sqrt_w(2) * (q_dot - q_dot_model); ...
         sqrt_w(3) * (r_dot - r_dot_model)];
end
