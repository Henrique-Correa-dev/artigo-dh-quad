function [fx, fy, fz] = accelerometer_model(p, q, r, u, v, w, T_m, p_dot, q_dot, r_dot, r_imu, kdm_in)
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
%   f_cg = R(ε)·[0; 0; -T_m]  (só thrust contribui — sem drag, sem gravidade),
%   com R(ε) o DESALINHAMENTO ANGULAR entre o plano dos rotores (eixo do
%   empuxo) e os eixos do sensor: parameters().imu_tilt = [eps_x; eps_y] (rad).
%   Um empuxo inclinado de eps_y em torno de y aparece no acelerômetro como
%   f_x = −T_m·sin(eps_y). Achado no drag_probe.m: o resíduo de a_x é
%   proporcional a T/m (R² 0,46), com eps_y ≈ 1,8°; no pairado isso vale
%   −0,31 m/s², que era o "bias" de x escrito à mão. Bias constante acerta o
%   nível médio mas erra sempre que T muda; o termo em T corrige isso.
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
    %  Com o desalinhamento eps (abaixo) o bias de x deixa de conter a parcela
    %  do empuxo inclinado no pairado (−T_m·sin(eps_y) ≈ −0,31 m/s²): o que
    %  sobra em x é ~0. Valores em parameters().imu_bias (fonte única).
    bias = pp_bias_default();

    rx = r_imu(1); ry = r_imu(2); rz = r_imu(3);

    %% Specific force no CG
    % Arrasto translacional (Eq. 6.31): f=-k_drag·v em x,y → força específica
    % −(k_drag/m)·u, −(k_drag/m)·v. Só thrust em Z. k_drag=0 → modelo antigo.
    pp = parameters();
    if nargin >= 12 && ~isempty(kdm_in)
        kdm = kdm_in;                         % = g·C_d, vindo de P(16)
    else
        kdm = pp.k_drag / pp.m;
        cd_ = getappdata(0,'diag_drag');
        if ~isempty(cd_), kdm = pp.g * cd_; end   % arrasto induzido do rotor (C_d, Beard 14.4.2)
    end
    % Desalinhamento empuxo↔sensor: empuxo T ao longo de −z do plano dos rotores,
    % visto pelo sensor girado de eps_x (em torno de x) e eps_y (em torno de y):
    %   f_x = −T_m·sin(eps_y),  f_y = +T_m·sin(eps_x)·cos(eps_y),  f_z = −T_m·cos(eps_x)·cos(eps_y)
    eps = [0;0]; if isfield(pp,'imu_tilt'), eps = pp.imu_tilt(:); end
    % kdm escalar (só arrasto de rotor) ou vetor de 3 (com as forças da estrutura)
    if size(kdm,2) >= 3        % Nx3 (um por amostra) ou 1x3
        kx = kdm(:,1); ky = kdm(:,2); kz = kdm(:,3);
    elseif numel(kdm) >= 3     % vetor coluna de 3
        kx = kdm(1);   ky = kdm(2);   kz = kdm(3);
    else                        % escalar: só arrasto de rotor, nada em z
        kx = kdm;      ky = kdm;      kz = 0;
    end
    fx_cg = -kx .* u - T_m .* sin(eps(2));                      % arrasto + inclinação
    fy_cg = -ky .* v + T_m .* sin(eps(1)) .* cos(eps(2));
    fz_cg = -kz .* w - T_m .* cos(eps(1)) .* cos(eps(2));

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

function b = pp_bias_default()
%PP_BIAS_DEFAULT  Bias DC do acelerômetro (m/s²) — lido de parameters().imu_bias.
    pp = parameters();
    if isfield(pp,'imu_bias'), b = pp.imu_bias(:);
    else,                      b = [0.0; -0.2; +0.4];   % fallback (x sem a parcela do empuxo)
    end
end
