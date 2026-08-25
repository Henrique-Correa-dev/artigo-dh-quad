function [rpm, fT, fQ, info] = motor_chain(time, pwm, varargin)
%MOTOR_CHAIN  Cadeia de atuação ÚNICA: PWM (log) → atraso → RPM regime → lag → RPM.
%
%   [rpm, fT, fQ, info] = motor_chain(time, pwm)
%   [rpm, fT, fQ, info] = motor_chain(time, pwm, 'delay_s',0.10, 'tau_m',0.05)
%
%   É a única definição do modelo de motor usada em identificação, validação e
%   figuras (identify_plant, sim_window, results_id_figures, compare_*). Antes,
%   cada script tinha a sua versão (spline × física, lag discreto × RPM como
%   estado, atraso no chamador × no script), e os R² não eram comparáveis.
%
%   Cadeia (Quan 2017; constantes de parameters().motor, ajustadas na bancada):
%     [0] atraso puro de transporte  PWM(t) → PWM(t − delay_s)      (ESC + log)
%     [1] estático linear            RPM_ss = CR·PWM + Ω_b
%     [2] lag de 1ª ordem            τ_m·dRPM/dt = RPM_ss − RPM     (spin-up)
%     [3] estático quadrático        T = kT_rpm·RPM², Q = kQ_rpm·RPM²  (fT, fQ)
%
%   O lag [2] é resolvido em forma FECHADA para entrada linear por partes
%   (o PWM do log é interpolado linearmente entre amostras), com condição
%   inicial em regime. Isso é exatamente o que vtol_dynamics_motor faz com a RPM
%   como estado, mas sem os 4 estados extras: os consumidores continuam
%   integrando 9 estados e recebem `rpm` no lugar de `pwm`, com fT/fQ agindo
%   sobre RPM. Como vtol_dynamics interpola o sinal linearmente entre amostras,
%   T(t) fica contínuo (não mais constante por 0,1 s).
%
%   Entradas
%     time   [N×1] grade regular (s)
%     pwm    [N×4] PWM dos 4 motores (μs)
%   Saídas
%     rpm    [N×4] rotação após atraso + lag
%     fT,fQ  handles: T = fT(rpm) [N], Q = fQ(rpm) [N·m]  (escalas k_T,k_Q ficam em P)
%     info   struct com delay_s, delay_n, tau_m, rpm_ss, e os coeficientes A,B,E
%
%   Escala por motor (P(5:8), P(9:12)) NÃO entra aqui: é aplicada por quem chama.

    mo = parameters().motor;
    % Modo utilitário: só os handles estáticos [3], p/ funções locais que já
    % recebem a RPM pronta:  [fT, fQ] = motor_chain('handles')
    if ischar(time) && strcmpi(time,'handles')
        rpm = @(w) mo.kT_rpm .* max(w,0).^2;    % 1ª saída = fT
        fT  = @(w) mo.kQ_rpm .* max(w,0).^2;    % 2ª saída = fQ
        fQ = []; info = [];
        return;
    end
    opt = struct('delay_s', mo.delay_s, 'tau_m', mo.tau_m);
    for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

    time = time(:);  N = numel(time);
    if size(pwm,1) ~= N && size(pwm,2) == N, pwm = pwm'; end
    dt = median(diff(time));

    % [0] atraso puro (shift inteiro de amostras, segurando o 1º valor)
    n_d = round(opt.delay_s / dt);
    if n_d > 0
        pwm = [repmat(pwm(1,:), n_d, 1); pwm(1:end-n_d, :)];
    end

    % [1] RPM de regime (saturação física: PWM em [1000,2000], RPM ≥ 0)
    rpm_ss = max(0, mo.CR .* min(max(pwm,1000),2000) + mo.Omega_b);

    % [2] lag de 1ª ordem, solução exata para entrada linear por partes:
    %     y[k+1] = E·y[k] + A·u[k] + B·(u[k+1] − u[k]),
    %     E = e^{−h/τ},  A = 1 − E,  B = A·(1 − τ/h) + E
    tau = opt.tau_m;
    if tau > 0
        E = exp(-dt/tau);  A = 1 - E;  B = A*(1 - tau/dt) + E;
        rpm = zeros(N,4);
        rpm(1,:) = rpm_ss(1,:);                       % IC em regime
        for k = 1:N-1
            rpm(k+1,:) = E*rpm(k,:) + A*rpm_ss(k,:) + B*(rpm_ss(k+1,:) - rpm_ss(k,:));
        end
    else
        E = 0; A = 1; B = 0;
        rpm = rpm_ss;
    end

    % [3] estático quadrático
    fT = @(w) mo.kT_rpm .* max(w,0).^2;
    fQ = @(w) mo.kQ_rpm .* max(w,0).^2;

    info = struct('delay_s',opt.delay_s,'delay_n',n_d,'tau_m',tau,'dt',dt, ...
                  'rpm_ss',rpm_ss,'E',E,'A',A,'B',B,'kT_rpm',mo.kT_rpm,'kQ_rpm',mo.kQ_rpm);
end
