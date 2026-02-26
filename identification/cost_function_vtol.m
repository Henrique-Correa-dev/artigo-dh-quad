function error_vector = cost_function_vtol(P, time_vec_exp, pqr_exp, initial_conditions, ...
                                           pwm_time, pwm_signals, func_T_ref, func_Q_ref, time_sim_output)
    % Função de custo para estimação de parâmetros do VTOL
    % P: vetor de parâmetros
    % time_vec_exp: vetor de tempo dos dados experimentais
    % pqr_exp: matriz de dados experimentais [p, q, r] (N_exp x 3)
    % initial_conditions: [p0; q0; r0] (condições iniciais para simulação)
    % pwm_time: vetor de tempo para os sinais PWM
    % pwm_signals: matriz dos sinais PWM [PWM1, PWM2, PWM3, PWM4]
    % func_T_ref: handle para a função de referência Torque = f(PWM)
    % func_Q_ref: handle para a função de referência TorqueReativo = f(PWM)
    % time_sim_output: pontos de tempo para avaliar a saída simulada (deve coincidir com time_vec_exp)

    % Define o handle da função ODE com os parâmetros P atuais
    ode_func = @(t,y) vtol_dynamics(t, y, P, pwm_time, pwm_signals, func_T_ref, func_Q_ref);

    % Opções do solver ODE (opcional, pode ser ajustado)
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

    % Simula o sistema
    % Para obter saída nos pontos de tempo específicos 'time_sim_output':
    [~, y_sim_full] = ode45(ode_func, time_sim_output, initial_conditions, options);
    
    p_sim = y_sim_full(:,1);
    q_sim = y_sim_full(:,2);
    r_sim = y_sim_full(:,3);

    % Dados experimentais
    p_exp = pqr_exp(:,1);
    q_exp = pqr_exp(:,2);
    r_exp = pqr_exp(:,3);

    % Vetor de erro (diferença entre simulado e experimental)
    % Pesos podem ser adicionados aqui se necessário (ex: se os sensores tiverem diferentes níveis de ruído)
    error_p = p_sim - p_exp;
    error_q = q_sim - q_exp;
    error_r = r_sim - r_exp;

    % Concatena em um único vetor de erro coluna
    error_vector = [error_p; error_q; error_r];
end