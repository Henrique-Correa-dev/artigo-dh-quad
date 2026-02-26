function torque_model_function = create_torque_model(pwm_experimental, torque_Nm_experimental, polynomial_degree)
    % CREATE_TORQUE_MODEL Ajusta um modelo polinomial aos dados de contra-torque.
    %
    % Inputs:
    %   pwm_experimental: Vetor coluna de valores PWM experimentais.
    %   torque_Nm_experimental: Vetor coluna de valores de contra-torque experimentais [Nm].
    %   polynomial_degree: Grau do polinômio a ser ajustado (ex: 2).
    %
    % Output:
    %   torque_model_function: Um handle de função que aceita PWM e retorna contra-torque [Nm].

    % Ajustar o polinômio: Q(pwm) = c(1)*pwm^N + ... + c(N+1)
    coeffs_torque = polyfit(pwm_experimental, torque_Nm_experimental, polynomial_degree);
    fprintf("%.15f \n",coeffs_torque)

    % Encontrar o PWM mínimo para o qual o torque é > 0 nos dados experimentais
    idx_first_active_torque = find(torque_Nm_experimental > 1e-9, 1, 'first'); % Usar uma pequena tolerância
    if isempty(idx_first_active_torque)
        min_pwm_for_active_torque = pwm_experimental(1); % Caso padrão
    else
        min_pwm_for_active_torque = pwm_experimental(idx_first_active_torque);
    end

    % Criar o handle da função do modelo de torque
    % A função garante que o torque seja >= 0 e respeite a "dead zone"
    torque_model_function = @(pwm_input) (pwm_input >= min_pwm_for_active_torque) .* max(0, polyval(coeffs_torque, pwm_input));

    fprintf('Coeficientes do Modelo de Torque (grau %d) [Nm]:\n', polynomial_degree);
    disp(coeffs_torque);
    fprintf('Torque será zero para PWM < %.0f com base nos dados fornecidos.\n', min_pwm_for_active_torque);

    % Opcional: Plotar para verificar o ajuste
    % figure;
    % plot(pwm_experimental, torque_Nm_experimental, 'o', 'DisplayName', 'Dados Experimentais (Nm)');
    % hold on;
    % pwm_range_plot = linspace(min(pwm_experimental), max(pwm_experimental), 200);
    % plot(pwm_range_plot, torque_model_function(pwm_range_plot), '-', 'DisplayName', sprintf('Modelo Polinomial (Grau %d)', polynomial_degree));
    % xlabel('PWM');
    % ylabel('Contra-Torque (Nm)');
    % legend show;
    % title('Ajuste do Modelo de Torque');
    % grid on;
end