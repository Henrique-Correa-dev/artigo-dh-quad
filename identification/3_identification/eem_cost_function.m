function e = eem_cost_function(P_scaled, P0_scale, weights, p, q, r, ...
    p_dot, q_dot, r_dot, T_ref, Q_ref, reg_lambda)
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
%   (Bz removido — modelo é rotacional puro agora)
%   (Bp, Bq, Br REMOVIDOS — CG offset capturado via Lx/Ly assimétricos)
%
%   reg_lambda (opcional, struct):
%     .kt_pair → penaliza |k_T1-k_T2| e |k_T3-k_T4| (default 0 = sem reg)
%     .kq_pair → penaliza |k_Q1-k_Q2| e |k_Q3-k_Q4| (default 0 = sem reg)
%     Maior λ → optimizer força motores do mesmo par CW/CCW ficarem similares.
%     Resolve identifiability quando r_dot único não distingue k_Q individuais.

    if nargin < 12, reg_lambda = struct('kt_pair', 0, 'kq_pair', 0); end
    if ~isfield(reg_lambda, 'kt_pair'), reg_lambda.kt_pair = 0; end
    if ~isfield(reg_lambda, 'kq_pair'), reg_lambda.kq_pair = 0; end

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

    % Braços lidos de parameters() — fonte única de verdade
    proj_p = parameters();
    Lx_r = proj_p.arms.Lx_r;
    Lx_l = proj_p.arms.Lx_l;
    Ly_f = proj_p.arms.Ly_f;
    Ly_r = proj_p.arms.Ly_r;

    Tmr = T_ref .* k_T';   % N x 4
    Qmr = Q_ref .* k_Q';   % N x 4

    Mx = -(Lx_r*Tmr(:,1) - Lx_l*Tmr(:,2) - Lx_l*Tmr(:,3) + Lx_r*Tmr(:,4));
    % My — ArduPilot QuadX (M1=FR, M2=RL, M3=FL, M4=RR — front/rear correto)
    My = Ly_f*Tmr(:,1) - Ly_r*Tmr(:,2) + Ly_f*Tmr(:,3) - Ly_r*Tmr(:,4);
    % Mz: pares DIAGONAIS padrão ArduPilot (M1+M2 CCW, M3+M4 CW)
    Mz = Qmr(:,1) + Qmr(:,2) - Qmr(:,3) - Qmr(:,4);

    p_dot_model = G1*p.*q - G2*q.*r + G3*Mx + G4*Mz - Dp*p;
    q_dot_model = G5*p.*r - G6*(p.^2 - r.^2) + invJy*My - Dq*q;
    r_dot_model = G7*p.*q - G1*q.*r + G4*Mx + G8*Mz - Dr*r;

    sqrt_w = sqrt(weights(:));
    e_dyn = [sqrt_w(1)*(p_dot - p_dot_model); ...
             sqrt_w(2)*(q_dot - q_dot_model); ...
             sqrt_w(3)*(r_dot - r_dot_model)];

    % Regularização: penaliza diferenças dentro do mesmo par CW/CCW
    %   M1+M2 pair (CCW):  k_*1 ≈ k_*2
    %   M3+M4 pair (CW) :  k_*3 ≈ k_*4
    e_reg = [sqrt(reg_lambda.kt_pair) * (k_T(1) - k_T(2)); ...
             sqrt(reg_lambda.kt_pair) * (k_T(3) - k_T(4)); ...
             sqrt(reg_lambda.kq_pair) * (k_Q(1) - k_Q(2)); ...
             sqrt(reg_lambda.kq_pair) * (k_Q(3) - k_Q(4))];

    e = [e_dyn; e_reg];
end
