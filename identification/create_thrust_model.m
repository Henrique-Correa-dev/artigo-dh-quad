function thrust_model_function = create_thrust_model(pwm_experimental, thrust_grams_experimental, polynomial_degree)
    % CREATE_THRUST_MODEL Ajusta um modelo polinomial aos dados de empuxo.
    %
    % Inputs:
    %   pwm_experimental: Vetor coluna de valores PWM experimentais.
    %   thrust_grams_experimental: Vetor coluna de valores de empuxo experimentais [gramas].
    %   polynomial_degree: Grau do polinômio a ser ajustado (ex: 2).
    %
    % Output:
    %   thrust_model_function: Um handle de função que aceita PWM e retorna empuxo [N].

    % Converter empuxo de gramas para Newtons
    thrust_N_experimental = thrust_grams_experimental * 9.80665 / 1000;

    % Ajustar o polinômio: T(pwm) = c(1)*pwm^N + ... + c(N+1)
    coeffs_thrust = polyfit(pwm_experimental, thrust_N_experimental, polynomial_degree);

    % Encontrar o PWM mínimo para o qual o empuxo é > 0 nos dados experimentais
    % Isso define uma "dead zone" para o modelo.
    idx_first_active_thrust = find(thrust_N_experimental > 1e-9, 1, 'first'); % Usar uma pequena tolerância
    if isempty(idx_first_active_thrust)
        min_pwm_for_active_thrust = pwm_experimental(1); % Caso padrão (improvável com seus dados)
    else
        min_pwm_for_active_thrust = pwm_experimental(idx_first_active_thrust);
    end
    
    % Criar o handle da função do modelo de empuxo
    % A função garante que o empuxo seja >= 0 e respeite a "dead zone"
    thrust_model_function = @(pwm_input) (pwm_input >= min_pwm_for_active_thrust) .* max(0, polyval(coeffs_thrust, pwm_input));
    
    fprintf('Coeficientes do Modelo de Empuxo (grau %d) [N]:\n', polynomial_degree);
    disp(coeffs_thrust);
    fprintf('Empuxo será zero para PWM < %.0f com base nos dados fornecidos.\n', min_pwm_for_active_thrust);
    
    % Opcional: Plotar para verificar o ajuste
    % figure;
    % plot(pwm_experimental, thrust_N_experimental, 'o', 'DisplayName', 'Dados Experimentais (N)');
    % hold on;
    % pwm_range_plot = linspace(min(pwm_experimental), max(pwm_experimental), 200);
    % plot(pwm_range_plot, thrust_model_function(pwm_range_plot), '-', 'DisplayName', sprintf('Modelo Polinomial (Grau %d)', polynomial_degree));
    % xlabel('PWM');
    % ylabel('Empuxo (N)');
    % legend show;
    % title('Ajuste do Modelo de Empuxo');
    % grid on;
end