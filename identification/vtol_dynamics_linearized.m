function dydt = vtol_dynamics_linearized(t, y, A, B, f0, x0, u0, pwm_time, pwm_signals)
    % Modelo linearizado do VTOL em torno do ponto de trim (hover)
    %
    % Equacao: dx/dt = f(x0,u0) + A*(x - x0) + B*(u - u0)
    %
    % onde:
    %   A = df/dx |_(x0,u0)   — Jacobiano dos estados (9x9)
    %   B = df/du |_(x0,u0)   — Jacobiano das entradas (9x4)
    %   f0 = f(x0, u0)        — residuo no trim (biases, etc.)
    %   x0 = estado de trim   (9x1)
    %   u0 = entrada de trim  (4x1, PWM hover)
    %
    % Estados: y = [p; q; r; phi; theta; psi; u; v; w]
    % Entradas: u = [PWM1; PWM2; PWM3; PWM4]
    %
    % Compativel com ode45:
    %   ode_lin = @(t,y) vtol_dynamics_linearized(t, y, A, B, f0, x0, u0, time, pwm);
    %   [t_s, y_s] = ode45(ode_lin, tspan, y0);

    % Interpolar sinais PWM no tempo atual t
    current_pwm = zeros(4, 1);
    for i = 1:4
        current_pwm(i) = interp1(pwm_time, pwm_signals(:,i), t, 'linear', 'extrap');
    end

    % Perturbacoes em relacao ao trim
    dx = y - x0;
    du = current_pwm - u0;

    % Dinamica linearizada
    dydt = f0 + A * dx + B * du;
end
