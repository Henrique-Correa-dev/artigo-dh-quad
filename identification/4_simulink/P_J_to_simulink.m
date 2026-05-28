function P_estimated = P_J_to_simulink(P_J)
%P_J_TO_SIMULINK  Converte P_J (15 elementos, J-formulation) para P_estimated (G-formulation, 20 elementos)
%
% Entrada (P_J), formato usado por identification/vtol_dynamics.m:
%   P_J(1:4)   = [Jx, Jy, Jz, Jxz]
%   P_J(5:8)   = k_T1..k_T4
%   P_J(9:12)  = k_Q1..k_Q4
%   P_J(13:15) = Dp, Dq, Dr
%
% Saída (P_estimated), formato usado por quad_model_v4.slx:
%   P_estimated(1:8)   = G1..G8        (computados a partir de Jx,Jy,Jz,Jxz)
%   P_estimated(9)     = invJy = 1/Jy
%   P_estimated(10:13) = k_T1..k_T4
%   P_estimated(14:17) = k_Q1..k_Q4
%   P_estimated(18:20) = Dp, Dq, Dr
%
% NOTA: v4.slx atualizado (sem biases, sem drag, sem Bz). P_estimated agora
% tem 20 elementos (era 23 — slots 21:23 de Bp/Bq/Br removidos).
%
% Variáveis avulsas (lidas pelo SLX fora do P_estimated):
%   Lx_r, Lx_l, Ly_f, Ly_r  (braços de momento, fonte: parameters().arms)
%   r_imu                     (offset CG→IMU, fonte: parameters().imu_offset)
%   bias_acc                  (bias DC do acelerômetro, hardcoded em accelerometer_model.m)
%   phi0, theta0, psi0        (ICs de atitude — preserva se já existir)
%
% Uso típico:
%   load('P_identified.mat', 'P_final');
%   P_estimated = P_J_to_simulink(P_final);
%   sim('quad_model_v4', ...);

    P_J = P_J(:);  % garantir coluna
    if numel(P_J) ~= 15
        error('P_J_to_simulink: P_J deve ter 15 elementos (recebeu %d).', numel(P_J));
    end

    %% Inércias → G-constants
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

    %% Montar P_estimated (20 entradas)
    P_estimated = zeros(20,1);
    P_estimated(1:8)   = [G1; G2; G3; G4; G5; G6; G7; G8];
    P_estimated(9)     = invJy;
    P_estimated(10:13) = P_J(5:8);     % k_T1..k_T4
    P_estimated(14:17) = P_J(9:12);    % k_Q1..k_Q4
    P_estimated(18:20) = P_J(13:15);   % Dp, Dq, Dr

    %% Variáveis avulsas que o SLX lê do workspace
    proj_p = parameters();

    % Braços de momento (substituem os antigos 0.232 ± dy_cg, 0.311 ± dx_cg)
    assignin('base', 'Lx_r', proj_p.arms.Lx_r);
    assignin('base', 'Lx_l', proj_p.arms.Lx_l);
    assignin('base', 'Ly_f', proj_p.arms.Ly_f);
    assignin('base', 'Ly_r', proj_p.arms.Ly_r);

    % IMU offset (CG → sensor, body frame) — pro accelerometer subsystem
    assignin('base', 'r_imu', proj_p.imu_offset);

    % Bias DC do acelerômetro — pegar do accelerometer_model.m
    % (lê o arquivo e procura a linha "bias = [...];")
    bias_acc = read_bias_from_acc_model();
    assignin('base', 'bias_acc', bias_acc);

    % Condições iniciais de atitude (preserva se já definidas pelo setup)
    for nm = {'phi0','theta0','psi0'}
        if ~evalin('base', sprintf('exist(''%s'',''var'')', nm{1}))
            assignin('base', nm{1}, 0);
        end
    end

    fprintf('P_J_to_simulink: P_estimated montado (%dx1).\n', numel(P_estimated));
    fprintf('   Workspace populado:\n');
    fprintf('     Lx_r=%.4f  Lx_l=%.4f  Ly_f=%.4f  Ly_r=%.4f\n', ...
            proj_p.arms.Lx_r, proj_p.arms.Lx_l, proj_p.arms.Ly_f, proj_p.arms.Ly_r);
    fprintf('     r_imu=[%.3f; %.3f; %.3f]\n', proj_p.imu_offset);
    fprintf('     bias_acc=[%.3f; %.3f; %.3f]\n', bias_acc);
end


function bias = read_bias_from_acc_model()
%READ_BIAS_FROM_ACC_MODEL  Lê o bias hardcoded em accelerometer_model.m
% (assim o Simulink fica em sync com o .m sem duplicar valores)
    paths = setup_paths();
    fid = fopen(fullfile(paths.model, 'accelerometer_model.m'), 'r');
    if fid < 0
        warning('Não consegui ler accelerometer_model.m. Usando bias=[0;0;0].');
        bias = [0; 0; 0];
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    bias = [0; 0; 0];
    while ~feof(fid)
        ln = fgetl(fid);
        if ischar(ln)
            % procurar "bias = [...];"
            tok = regexp(strtrim(ln), '^bias\s*=\s*\[([^\]]+)\]\s*;', 'tokens', 'once');
            if ~isempty(tok)
                vals = sscanf(tok{1}, '%f;');
                if numel(vals) == 3
                    bias = vals(:);
                    return;
                end
            end
        end
    end
end
