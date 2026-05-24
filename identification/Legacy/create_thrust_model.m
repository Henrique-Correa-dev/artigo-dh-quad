function thrust_model_function = create_thrust_model(pwm_experimental, thrust_grams_experimental, varargin)
%CREATE_THRUST_MODEL  Modelo de empuxo PWM->T via spline cubica da tabela de bancada.
%
% Inputs:
%   pwm_experimental:          vetor PWM dos pontos de bancada
%   thrust_grams_experimental: vetor de empuxo em GRAMAS
%   (3o arg, opcional, ignorado: mantido por compat com chamadas antigas)
%
% Output:
%   thrust_model_function: handle T(pwm) [N]
%     - PWM saturado em [1000, 2000] us
%     - Interpolacao SPLINE CUBICA passando exato pelos pontos de bancada
%     - Saida >= 0 (max(0, .))
%
% NOTA: o terceiro argumento (anteriormente polynomial_degree) e mantido por
% compatibilidade com codigo existente mas IGNORADO. O modelo nao usa mais
% polyfit/polyval -- usa interp1 'spline'.

    % Converter empuxo de gramas para Newtons
    thrust_N_experimental = thrust_grams_experimental * 9.80665 / 1000;

    pwm_min = 1000;
    pwm_max = 2000;

    % Handle do modelo: spline cubica, saturando PWM e clampando T>=0
    pwm_bp = pwm_experimental(:);
    T_bp   = thrust_N_experimental(:);
    thrust_model_function = @(pwm_input) ...
        max(0, interp1(pwm_bp, T_bp, ...
                       min(max(pwm_input, pwm_min), pwm_max), ...
                       'makima', 'extrap'));

    fprintf('Modelo de Empuxo: Akima (makima) sobre %d pontos\n', numel(pwm_bp));
    fprintf('  PWM saturado em [%d, %d] us; empuxo >= 0\n', pwm_min, pwm_max);
    fprintf('  Bench: PWM=[%s]  T=[%s] N\n', ...
        strjoin(arrayfun(@(x) sprintf('%g',x), pwm_bp, 'UniformOutput', false), ' '), ...
        strjoin(arrayfun(@(x) sprintf('%.3f',x), T_bp, 'UniformOutput', false), ' '));
end
