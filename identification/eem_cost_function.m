function e = eem_cost_function(P_scaled, P0_scale, weights, p, q, r, ...
    p_dot, q_dot, r_dot, T_ref, Q_ref)
% EEM_COST_FUNCTION  Equation Error Method cost for VTOL identification.
%
%   No ODE integration. Uses measured states and numerically computed
%   derivatives directly. The residual is the mismatch between measured
%   derivatives and model-predicted derivatives at each time step.
%
%   P(1:4)   = [Jx, Jy, Jz, Jxz]    (momentos de inércia)
%   P(5:8)   = k_T1..k_T4
%   P(9:12)  = k_Q1..k_Q4
%   P(13:15) = Dp, Dq, Dr
%   P(16:18) = Bp, Bq, Br
%   P(19:20) = dx_cg, dy_cg  (offset do CG vs CAD)

    P = P_scaled(:) .* P0_scale(:);

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
    Dp = P(13); Dq = P(14); Dr = P(15);
    Bp = P(16); Bq = P(17); Br = P(18);

    % CG offsets (P pode ter 18 ou 20 elementos)
    if length(P) >= 20
        dx_cg = P(19); dy_cg = P(20);
    else
        dx_cg = 0; dy_cg = 0;
    end

    % Braços efetivos com offset do CG
    Lx_r = 0.232 - dy_cg;   % direita (motores 1,4)
    Lx_l = 0.232 + dy_cg;   % esquerda (motores 2,3)
    Ly_f = 0.311185 - dx_cg; % frente (motores 1,3)
    Ly_r = 0.342865 + dx_cg; % traseira (motores 2,4)

    Tmr = T_ref .* k_T';   % N x 4
    Qmr = Q_ref .* k_Q';   % N x 4

    Mx = -(Lx_r*Tmr(:,1) - Lx_l*Tmr(:,2) - Lx_l*Tmr(:,3) + Lx_r*Tmr(:,4));
    My = Ly_f*Tmr(:,1) - Ly_r*Tmr(:,2) + Ly_f*Tmr(:,3) - Ly_r*Tmr(:,4);
    Mz = Qmr(:,1) + Qmr(:,2) - Qmr(:,3) - Qmr(:,4);

    p_dot_model = G1*p.*q - G2*q.*r + G3*Mx + G4*Mz - Dp*p + Bp;
    q_dot_model = G5*p.*r - G6*(p.^2 - r.^2) + invJy*My - Dq*q + Bq;
    r_dot_model = G7*p.*q - G1*q.*r + G4*Mx + G8*Mz - Dr*r + Br;

    sqrt_w = sqrt(weights(:));
    e = [sqrt_w(1)*(p_dot - p_dot_model); ...
         sqrt_w(2)*(q_dot - q_dot_model); ...
         sqrt_w(3)*(r_dot - r_dot_model)];
end
