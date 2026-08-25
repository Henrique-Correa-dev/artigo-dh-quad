function dydt = vtol_dynamics_motor(t, y, P, pwm_time, pwm_signals, constants)
%VTOL_DYNAMICS_MOTOR  Dinâmica VTOL com MODELO DE MOTOR FÍSICO (RPM como estado).
%
%  Estende vtol_dynamics adicionando 4 estados de RPM (1 por motor):
%      estado y = [p q r  phi theta psi  u v w  RPM1 RPM2 RPM3 RPM4]  (13x1)
%
%  Cadeia do motor (Quan 2017 / rotorcraft), constantes de parameters().motor:
%    [1] RPM_ss_i = CR·PWM_i + Omega_b            (estático linear, bancada)
%    [2] d(RPM_i)/dt = (RPM_ss_i − RPM_i)/tau_m   (lag 1ª ordem)
%    [3] T_i = k_T_i·kT_rpm·RPM_i² , Q_i = k_Q_i·kQ_rpm·RPM_i²  (estático quadr.)
%  k_T_i = P(5:8), k_Q_i = P(9:12) entram como ESCALA POR MOTOR.
%  O atraso de transporte (tau_d ≈ 1 amostra) é aplicado FORA, no PWM (DELAY_PWM).
%
%  Reusa as subfunções de vtol_dynamics('get_handles') → corpo rígido idêntico.

    mp = parameters();
    CR = mp.motor.CR;  Ob = mp.motor.Omega_b;
    kT = mp.motor.kT_rpm;  kQ = mp.motor.kQ_rpm;  tau = mp.motor.tau_m;
    Lx_r = mp.arms.Lx_r;  Lx_l = mp.arms.Lx_l;  Ly_f = mp.arms.Ly_f;  Ly_r = mp.arms.Ly_r;
    kt = P(5:8);  kq = P(9:12);

    dyn = vtol_dynamics('get_handles');

    % PWM interpolado em t
    pwm = zeros(1,4);
    for i = 1:4
        pwm(i) = interp1(pwm_time, pwm_signals(:,i), t, 'linear', 'extrap');
    end

    % --- Estados de RPM ---
    w   = y(10:13)';                       % 1x4  (RPM atual)
    w_ss = max(0, CR*pwm + Ob);            % RPM de regime (sat ≥ 0)
    w_dot = (w_ss - w) / tau;              % lag 1ª ordem

    % --- Estágio [3]: RPM → empuxo/torque (escala por motor) ---
    Tmr = kt(:)' .* (kT * w.^2);           % 1x4  [N]
    Qmr = kq(:)' .* (kQ * w.^2);           % 1x4  [N·m]

    [Mx, My, Mz] = dyn.moments(Tmr, Qmr, Lx_r, Lx_l, Ly_f, Ly_r);
    T_total = sum(Tmr, 2);

    % --- Corpo rígido (9 estados) via subfunção (fonte única) ---
    drb = dyn.rigid_body(y(1:9), T_total, Mx, My, Mz, P, constants);

    dydt = [drb; w_dot(:)];
end
