function P_estimated = P_J_to_simulink(P_J)
%P_J_TO_SIMULINK  Converte vetor P (formulacao J) para P_estimated (formulacao G)
%
% Entrada (P_J), formato usado por identification/vtol_dynamics.m:
%   P_J(1:4)   = [Jx, Jy, Jz, Jxz]
%   P_J(5:8)   = k_T1..k_T4
%   P_J(9:12)  = k_Q1..k_Q4
%   P_J(13:15) = Dp, Dq, Dr
%   P_J(16:18) = Bp, Bq, Br
%   P_J(19:22) = Xu, Yv, Zw, Bz        (opcional, defaults -4,-4,-0.1,-0.5)
%
% NOTA: dx_cg, dy_cg REMOVIDOS do vetor (CG oficial = onde dx=dy=0).
% v4.slx ainda lê dx_cg/dy_cg do workspace — esta função força ambos a 0.
%
% Saida (P_estimated), formato usado por quad_model_v5.slx:
%   P_estimated(1:8)   = G1..G8        (computados a partir de Jx,Jy,Jz,Jxz)
%   P_estimated(9)     = invJy = 1/Jy
%   P_estimated(10:13) = k_T1..k_T4
%   P_estimated(14:17) = k_Q1..k_Q4
%   P_estimated(18:20) = Dp, Dq, Dr
%   P_estimated(21:23) = Bp, Bq, Br    (reservado, nao usado pelo SLX atual)
%
% Os seguintes valores tambem sao assignados no base workspace (o SLX
% le diretamente dali, fora do P_estimated):
%   p_bias = Bp, q_bias = Bq, r_bias = Br
%   dx_cg, dy_cg
%   Xu_param, Yv_param, Zw_param, Bz_param
%   phi0=0, theta0=0, psi0=0  (somente se ainda nao existirem)
%
% Uso tipico:
%   load('resultado_identification.mat', 'P_final');   % do new_identification.m
%   P_estimated = P_J_to_simulink(P_final);
%   % Agora pode-se simular o quad_model_v5.slx.

    P_J = P_J(:);  % garantir coluna

    %% Inercias -> G-constants
    Jx  = P_J(1);
    Jy  = P_J(2);
    Jz  = P_J(3);
    Jxz = P_J(4);
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

    %% Montar P_estimated (23 entradas)
    P_estimated = zeros(23,1);
    P_estimated(1:8)   = [G1; G2; G3; G4; G5; G6; G7; G8];
    P_estimated(9)     = invJy;
    P_estimated(10:13) = P_J(5:8);     % k_T1..k_T4
    P_estimated(14:17) = P_J(9:12);    % k_Q1..k_Q4
    P_estimated(18:20) = P_J(13:15);   % Dp, Dq, Dr
    P_estimated(21:23) = P_J(16:18);   % Bp, Bq, Br (slot reservado)

    %% Variaveis avulsas que o SLX le do workspace
    % Biases rotacionais
    assignin('base', 'p_bias', P_J(16));
    assignin('base', 'q_bias', P_J(17));
    assignin('base', 'r_bias', P_J(18));

    % CG: forçado a 0 (oficial do CAD). v4.slx ainda lê esses do workspace.
    assignin('base', 'dx_cg', 0);
    assignin('base', 'dy_cg', 0);

    % Parametros translacionais (opcional, com defaults)
    if numel(P_J) >= 22
        assignin('base', 'Xu_param', P_J(19));
        assignin('base', 'Yv_param', P_J(20));
        assignin('base', 'Zw_param', P_J(21));
        assignin('base', 'Bz_param', P_J(22));
    else
        assignin('base', 'Xu_param', -4.0);
        assignin('base', 'Yv_param', -4.0);
        assignin('base', 'Zw_param', -0.1);
        assignin('base', 'Bz_param', -0.5);
    end

    % Condicoes iniciais de atitude (somente se nao definidas)
    for nm = {'phi0','theta0','psi0'}
        if ~evalin('base', sprintf('exist(''%s'',''var'')', nm{1}))
            assignin('base', nm{1}, 0);
        end
    end

    fprintf('P_J_to_simulink: P_estimated montado (%dx1).\n', numel(P_estimated));
    fprintf('   Workspace populado: p_bias, q_bias, r_bias, dx_cg, dy_cg,\n');
    fprintf('                       Xu_param, Yv_param, Zw_param, Bz_param.\n');
end
