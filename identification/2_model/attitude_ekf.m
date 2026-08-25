function est = attitude_ekf(time, gyro, accel, opts)
%ATTITUDE_EKF  EKF contínuo-discreto de 2 estados [phi, theta] (roll, pitch).
%
% Estima atitude fundindo gyro (predição) e acelerômetro (correção pela
% gravidade). Cura a deriva da integração open-loop da cinemática de Euler.
% Formulação Beard & McLain, cap. 8 (attitude estimation).
%
% NÃO estima yaw — gravidade não observa yaw (precisa de magnetômetro).
% Yaw é tratado num bloco separado (ver attitude_yaw / mag).
%
% ── MODELO ──────────────────────────────────────────────────────────────
% Estado:    x = [phi; theta]                       (rad)
% Entrada:   u = [p; q; r]  (gyro, rad/s)            — entra na predição
% Medida:    y = [ax; ay; az]  (accel, m/s²)         — entra na correção
%
% Predição (cinemática de Euler, gyro como entrada):
%   phi_dot   = p + (q·sinφ + r·cosφ)·tanθ
%   theta_dot = q·cosφ − r·sinφ
%
% Medição (força específica em quase-equilíbrio, Va≈0):
%   h(x) = [ +g·sinθ;  −g·sinφ·cosθ;  −g·cosφ·cosθ ]
%   (em hover nivelado → [0; 0; −g], bate com o accel medido)
%
% ── GATE DE CONFIANÇA (a "zona morta" feita certo) ──────────────────────
% O accel só é referência válida de gravidade perto do equilíbrio. Em
% manobra agressiva |a| afasta de g e a medida fica corrompida. Então R
% (ruído de medição) é INFLADO suavemente conforme ||a|−g| cresce:
%   R_eff = R · (1 + (max(0, ||a|−g| − gate_lo)/gate_scale)²)
% Perto do equilíbrio (||a|−g| < gate_lo): R normal, accel confiável.
% Em manobra: R cresce, EKF confia mais no gyro. Sem descontinuidade.
%
% ── USO ─────────────────────────────────────────────────────────────────
%   est = attitude_ekf(time, gyro, accel);                    % defaults
%   est = attitude_ekf(time, gyro, accel, struct('Q',Q,'R',R,'bias_gyro',bg));
%
% ENTRADAS:
%   time   [Nx1]  tempo (s)
%   gyro   [Nx3]  [p q r] rad/s
%   accel  [Nx3]  [ax ay az] m/s²
%   opts   struct opcional:
%     .g          (9.81)            gravidade
%     .Q          (2x2)             ruído de processo (default diag([1e-4 1e-4]))
%     .R          (3x3 ou escalar)  ruído de medição accel (default (0.5)²·I)
%     .bias_gyro  [3x1] (0)         bias de gyro a subtrair (de estimate_bias)
%     .bias_accel [3x1] (0)         bias de accel a subtrair
%     .phi0,.theta0 (auto)          IC; default = ângulos do 1º accel
%     .gate_lo    (0.5)             ||a|-g| abaixo disto: accel 100% confiável [m/s²]
%     .gate_scale (1.0)             escala de inflação de R [m/s²]
%     .n_sub      (5)               sub-passos de predição por amostra
%
% SAÍDA (struct est):
%   .phi, .theta        [Nx1] estimativa (rad)
%   .phi_deg,.theta_deg [Nx1] estimativa (graus)
%   .P_trace            [Nx1] traço da covariância (diagnóstico de convergência)
%   .innov_norm         [Nx1] norma da inovação (accel medido − previsto)
%   .R_inflation        [Nx1] fator de inflação de R aplicado (1=sem inflação)
%   .opts_used          struct com opções efetivas

    %% Defaults
    d = struct('g', 9.81, ...
               'Q', diag([1e-4, 1e-4]), ...
               'R', 1.0 * eye(3), ...   % σ_acc≈1 m/s² — tuned: RMSE φ,θ ≈1.4° (vs 4° open-loop)
               'bias_gyro',  [0;0;0], ...
               'bias_accel', [0;0;0], ...
               'phi0', [], 'theta0', [], ...
               'gate_lo', 0.5, 'gate_scale', 1.0, ...
               'n_sub', 5);
    if nargin < 4, opts = struct(); end
    opts = merge_opts_(d, opts);
    if isscalar(opts.R), opts.R = opts.R * eye(3); end

    g = opts.g;
    N = numel(time);

    % Subtrair bias dos sensores
    gyro  = gyro  - opts.bias_gyro(:)';
    accel = accel - opts.bias_accel(:)';

    %% Condição inicial: ângulos do 1º accel (se não fornecido)
    a0 = accel(1,:)';
    if isempty(opts.phi0)
        phi = atan2(-a0(2), -a0(3));
    else
        phi = opts.phi0;
    end
    if isempty(opts.theta0)
        theta = atan2(a0(1), hypot(a0(2), a0(3)));
    else
        theta = opts.theta0;
    end

    x = [phi; theta];
    P = diag([0.1, 0.1]);   % incerteza inicial moderada

    %% Saídas
    est.phi = zeros(N,1);  est.theta = zeros(N,1);
    est.P_trace = zeros(N,1);  est.innov_norm = zeros(N,1);
    est.R_inflation = ones(N,1);
    est.phi(1) = x(1);  est.theta(1) = x(2);  est.P_trace(1) = trace(P);

    I2 = eye(2);

    %% Loop EKF
    for k = 1:N-1
        dt = time(k+1) - time(k);
        if dt <= 0 || ~isfinite(dt), dt = 1e-3; end

        % ---- PREDIÇÃO (sub-stepping) ----
        u = gyro(k,:)';   % [p;q;r] no instante k (ZOH)
        h_sub = dt / opts.n_sub;
        for s = 1:opts.n_sub
            [f, A] = euler_kin_(x, u);
            x = x + h_sub * f;
            P = P + h_sub * (A*P + P*A' + opts.Q);
        end

        % ---- CORREÇÃO (accel em k+1) ----
        y = accel(k+1,:)';
        a_norm = norm(y);

        % Gate: inflar R conforme afastamento de g
        excess = max(0, abs(a_norm - g) - opts.gate_lo);
        infl   = 1 + (excess / opts.gate_scale)^2;
        R_eff  = opts.R * infl;

        [hx, C] = accel_meas_(x, g);
        innov = y - hx;
        S = C*P*C' + R_eff;
        K = (P*C') / S;
        x = x + K*innov;
        P = (I2 - K*C)*P;
        P = 0.5*(P + P');   % manter simétrica

        % Armazenar
        est.phi(k+1)   = x(1);
        est.theta(k+1) = x(2);
        est.P_trace(k+1)    = trace(P);
        est.innov_norm(k+1) = norm(innov);
        est.R_inflation(k+1)= infl;
    end

    est.phi_deg   = rad2deg(est.phi);
    est.theta_deg = rad2deg(est.theta);
    est.opts_used = opts;
end


%% ========================================================================
%  SUBFUNÇÕES
%  ========================================================================
function [f, A] = euler_kin_(x, u)
%EULER_KIN_  Cinemática de Euler f(x,u) e Jacobiano A = ∂f/∂x.
    phi = x(1); theta = x(2);
    p = u(1); q = u(2); r = u(3);

    sp = sin(phi); cp = cos(phi);
    tt = tan(theta); ct = cos(theta);
    sec2 = 1/ct^2;     % sec(theta)^2

    f = [ p + (q*sp + r*cp)*tt;
          q*cp - r*sp ];

    A = [ (q*cp - r*sp)*tt,   (q*sp + r*cp)*sec2;
          -q*sp - r*cp,       0 ];
end

function [h, C] = accel_meas_(x, g)
%ACCEL_MEAS_  Modelo de medição do accel h(x) e Jacobiano C = ∂h/∂x.
%   Força específica em quase-equilíbrio (gravidade-reação).
    phi = x(1); theta = x(2);
    sp = sin(phi); cp = cos(phi);
    st = sin(theta); ct = cos(theta);

    h = [  g*st;
          -g*sp*ct;
          -g*cp*ct ];

    C = [ 0,            g*ct;
         -g*cp*ct,      g*sp*st;
          g*sp*ct,      g*cp*st ];
end

function opts = merge_opts_(d, u)
    opts = d;
    if isstruct(u)
        fn = fieldnames(u);
        for k = 1:numel(fn), opts.(fn{k}) = u.(fn{k}); end
    end
end
