function [fx, fy, fz] = accelerometer_model(p, q, r, u, v, w, T_m, p_dot, q_dot, r_dot, r_imu)
%ACCELEROMETER_MODEL  Modelo de sensor de acelerômetro do drone.
%
% Retorna specific force (m/s²) lida pelo sensor montado em ponto offset
% r_imu do CG. Drag NÃO modelado (Xu=Yv=Zw=0 removidos por design).
% Bias DC do sensor: hardcoded como constante neste arquivo (ajustar à mão).
%
%   f_imu = f_cg + α × r_imu + ω × (ω × r_imu)
%             ↑          ↑              ↑
%        proof mass   Euler       centrípeta
%
% Específico do drone:
%   f_cg = [0; 0; -T_m]    (só thrust contribui — sem drag, sem gravidade)
%
% A gravidade NÃO aparece porque o acelerômetro mede SPECIFIC FORCE
% (não-gravitacional). Derivação: f = a_inertial - g_body. Substituindo
% u̇ da EOM, gravidade e Coriolis cancelam.
%
% INPUTS (todos podem ser escalar OU vetor Nx1):
%   p, q, r              : velocidade angular [rad/s]
%   u, v, w              : velocidade body [m/s]      (sem efeito sem drag)
%   T_m                  : T_total/m [m/s²]
%   p_dot, q_dot, r_dot  : aceleração angular [rad/s²] (opcional, default 0)
%   r_imu                : [rx; ry; rz] vetor CG→IMU body frame [m] (opcional, default [0;0;0])
%
% OUTPUTS:
%   fx, fy, fz           : specific force no IMU [m/s²]

    % Defaults para argumentos opcionais
    if nargin < 11 || isempty(r_imu)
        r_imu = [0; 0; 0];
    end
    if nargin < 10 || isempty(p_dot)
        p_dot = zeros(size(p)); q_dot = zeros(size(p)); r_dot = zeros(size(p));
    end

    %% ====================================================================
    %  BIAS DC do sensor — HARDCODED (ajustar à mão conforme análise visual)
    %  Em hover ideal: f = [0; 0; -g]. Bias = mean(acc_IMU_hover) - [0;0;-g].
    %  Não entra em P, não é identificado — é propriedade do hardware.
    %  ====================================================================
    bias = [-0.3; -0.2; +0.4];   % [bx; by; bz] em m/s²

    rx = r_imu(1); ry = r_imu(2); rz = r_imu(3);

    %% Specific force no CG
    % Sem drag (Xu=Yv=Zw=0), sem bias (Bz vai pra sensor model).
    % Só thrust em Z body.
    fx_cg = 0 .* u;     % drag X removido → escala com u só pra preservar dim Nx1
    fy_cg = 0 .* v;
    fz_cg = -T_m;

    %% Correção por IMU offset: a_imu = a_cg + α × r + ω × (ω × r)
    % Termo Euler (α × r)
    eul_x = q_dot.*rz - r_dot.*ry;
    eul_y = r_dot.*rx - p_dot.*rz;
    eul_z = p_dot.*ry - q_dot.*rx;

    % Termo centrípeto (ω × (ω × r))
    cen_x = p.*q.*ry + p.*r.*rz - (q.^2 + r.^2).*rx;
    cen_y = p.*q.*rx + q.*r.*rz - (p.^2 + r.^2).*ry;
    cen_z = p.*r.*rx + q.*r.*ry - (p.^2 + q.^2).*rz;

    %% Saída (com bias hardcoded somado)
    fx = fx_cg + eul_x + cen_x + bias(1);
    fy = fy_cg + eul_y + cen_y + bias(2);
    fz = fz_cg + eul_z + cen_z + bias(3);
end
