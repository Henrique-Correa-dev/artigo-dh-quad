function torque_model_function = create_torque_model(pwm_experimental, torque_Nm_experimental, varargin)
%CREATE_TORQUE_MODEL  Modelo de contra-torque PWM->Q via spline cubica da bancada.
%
% Inputs:
%   pwm_experimental:        vetor PWM dos pontos de bancada
%   torque_Nm_experimental:  vetor de torque em N.m
%   (3o arg, opcional, ignorado: mantido por compat com chamadas antigas)
%
% Output:
%   torque_model_function: handle Q(pwm) [N.m]
%     - PWM saturado em [1000, 2000] us
%     - Interpolacao SPLINE CUBICA passando exato pelos pontos de bancada
%     - Saida >= 0 (max(0, .))

    pwm_min = 1000;
    pwm_max = 2000;

    pwm_bp = pwm_experimental(:);
    Q_bp   = torque_Nm_experimental(:);
    torque_model_function = @(pwm_input) ...
        max(0, interp1(pwm_bp, Q_bp, ...
                       min(max(pwm_input, pwm_min), pwm_max), ...
                       'makima', 'extrap'));

    fprintf('Modelo de Torque: Akima (makima) sobre %d pontos\n', numel(pwm_bp));
    fprintf('  PWM saturado em [%d, %d] us; torque >= 0\n', pwm_min, pwm_max);
    fprintf('  Bench: PWM=[%s]  Q=[%s] Nm\n', ...
        strjoin(arrayfun(@(x) sprintf('%g',x), pwm_bp, 'UniformOutput', false), ' '), ...
        strjoin(arrayfun(@(x) sprintf('%.4f',x), Q_bp, 'UniformOutput', false), ' '));
end
