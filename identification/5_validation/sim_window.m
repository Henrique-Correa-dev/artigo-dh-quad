function res = sim_window(mode, P, time, pwm, pqr_meas, att_meas, constants)
%SIM_WINDOW  Simula uma janela do voo em um dos 3 modos: full / hybrid / semi.
%
% USO:
%   res = sim_window('full',   P, time, pwm, pqr_meas, att_meas, constants);
%   res = sim_window('hybrid', P, time, pwm, pqr_meas, att_meas, constants);
%   res = sim_window('semi',   P, time, pwm, pqr_meas, att_meas, constants);
%
% MODOS:
%   'full'   → integra TUDO (p,q,r,φ,θ,ψ,u,v,w). Modelo puro.
%   'hybrid' → integra p,q,r e u,v,w. Usa att_meas pra gravidade. Sem drift de att.
%   'semi'   → passthrough de pqr_meas. Integra só u,v,w com att_meas.
%
% ARQUITETURA:
%   moments e trans_dot vêm de vtol_dynamics('get_handles') — fonte única.
%   accelerometer_model (sensor) é arquivo separado em 2_model/.
%
% ENTRADAS:
%   P         [15x1] vetor de parâmetros (P_J)
%   time      [Nx1] grade de tempo (s)
%   pwm       [Nx4] PWMs dos 4 motores (µs)
%   pqr_meas  [Nx3] gyro medido [rad/s]
%   att_meas  [Nx3] atitude medida do EKF [graus]
%   constants struct: .m .g .tau_motor
%
% SAÍDA (struct res):
%   .p, .q, .r          [Nx1] velocidades angulares simuladas
%   .phi, .theta, .psi  [Nx1] atitude simulada (graus)
%   .u, .v, .w          [Nx1] velocidades body
%   .accX, .accY, .accZ [Nx1] força específica IMU (do accelerometer_model)
%   .T_total            [Nx1] empuxo total
%   .Mx, .My, .Mz       [Nx1] momentos

    %% Setup
    N  = length(time);
    dt = time(2) - time(1);
    g  = constants.g;
    m  = constants.m;

    % Modelo de motor
    [func_T, func_Q] = motor_models();

    % Parâmetros do P
    k_T = P(5:8);
    k_Q = P(9:12);

    % Braços e IMU offset (de parameters.m)
    proj_p_sw = parameters();
    Lx_r = proj_p_sw.arms.Lx_r;  Lx_l = proj_p_sw.arms.Lx_l;
    Ly_f = proj_p_sw.arms.Ly_f;  Ly_r = proj_p_sw.arms.Ly_r;
    r_imu = proj_p_sw.imu_offset;

    % Handles centralizados do vtol_dynamics (fonte única)
    dyn_h = vtol_dynamics('get_handles');
    moments_fn   = dyn_h.moments;
    trans_dot_fn = dyn_h.trans_dot;

    %% T_ref, Q_ref e momentos (vetorizado sobre toda a janela)
    T_ref = zeros(N,4);  Q_ref = zeros(N,4);
    for j = 1:4
        T_ref(:,j) = func_T(pwm(:,j));
        Q_ref(:,j) = func_Q(pwm(:,j));
    end
    Tmr = T_ref .* k_T';
    Qmr = Q_ref .* k_Q';
    T_total = sum(Tmr, 2);

    % Momentos via subfunção (forma vetorizada — moments_local detecta Nx4)
    [Mx, My, Mz] = moments_fn(Tmr, Qmr, Lx_r, Lx_l, Ly_f, Ly_r);

    res.Mx = Mx;  res.My = My;  res.Mz = Mz;  res.T_total = T_total;

    %% ============= ROTACIONAL =============
    switch lower(mode)
        case 'full'
            % Integra p,q,r,φ,θ,ψ,u,v,w via vtol_dynamics (9 estados)
            att_rad0 = deg2rad(att_meas(1,:));
            y0 = [pqr_meas(1,:)'; att_rad0(:); 0; 0; 0];
            ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
            ode_func = @(t,y) vtol_dynamics(t, y, P, time, pwm, func_T, func_Q, constants);
            [t_s, y_s] = ode45(ode_func, time, y0, ode_opts);
            y_out = interp1(t_s, y_s, time, 'linear', 'extrap');

            res.p   = y_out(:,1);  res.q     = y_out(:,2);  res.r   = y_out(:,3);
            res.phi = rad2deg(y_out(:,4));
            res.theta = rad2deg(y_out(:,5));
            res.psi = rad2deg(y_out(:,6));
            res.u = y_out(:,7);  res.v = y_out(:,8);  res.w = y_out(:,9);

        case 'hybrid'
            % Integra só p,q,r via vtol_dynamics em modo 3-estados.
            % Atitude vem da medida.
            y0 = pqr_meas(1,:)';
            ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
            ode_func = @(t,y) vtol_dynamics(t, y, P, time, pwm, func_T, func_Q, constants);
            [t_s, y_s] = ode45(ode_func, time, y0, ode_opts);
            y_out = interp1(t_s, y_s, time, 'linear', 'extrap');
            res.p = y_out(:,1); res.q = y_out(:,2); res.r = y_out(:,3);
            res.phi   = att_meas(:,1);
            res.theta = att_meas(:,2);
            res.psi   = att_meas(:,3);

        case 'semi'
            % p,q,r passthrough do medido (não integra rotacional)
            res.p = pqr_meas(:,1);  res.q = pqr_meas(:,2);  res.r = pqr_meas(:,3);
            res.phi   = att_meas(:,1);
            res.theta = att_meas(:,2);
            res.psi   = att_meas(:,3);

        otherwise
            error('sim_window: modo desconhecido "%s". Use full|hybrid|semi.', mode);
    end

    %% ============= TRANSLACIONAL (u, v, w) =============
    % Gravidade no body frame
    if strcmp(mode, 'full')
        phi_use   = deg2rad(res.phi);
        theta_use = deg2rad(res.theta);
    else
        phi_use   = deg2rad(att_meas(:,1));
        theta_use = deg2rad(att_meas(:,2));
    end
    gx = -g * sin(theta_use);
    gy =  g * cos(theta_use) .* sin(phi_use);
    gz =  g * cos(theta_use) .* cos(phi_use);

    % Pra full, u,v,w já vieram do ode45. Pra hybrid/semi, integra aqui via subfunção.
    if ~strcmp(mode, 'full')
        n_sub  = 5;
        dt_sub = dt / n_sub;
        h2     = dt_sub / 2;
        u = zeros(N,1); v = zeros(N,1); w = zeros(N,1);
        for k = 1:N-1
            pk = res.p(k); qk = res.q(k); rk = res.r(k);
            gxk = gx(k); gyk = gy(k); gzk = gz(k);
            Tk_m = T_total(k) / m;
            us = u(k); vs = v(k); ws = w(k);

            % RK4 sub-stepping (chama trans_dot_fn — fonte única)
            for si = 1:n_sub
                [ud1,vd1,wd1] = trans_dot_fn(pk,qk,rk, us,vs,ws, gxk,gyk,gzk, Tk_m);
                u2=us+h2*ud1; v2=vs+h2*vd1; w2=ws+h2*wd1;

                [ud2,vd2,wd2] = trans_dot_fn(pk,qk,rk, u2,v2,w2, gxk,gyk,gzk, Tk_m);
                u3=us+h2*ud2; v3=vs+h2*vd2; w3=ws+h2*wd2;

                [ud3,vd3,wd3] = trans_dot_fn(pk,qk,rk, u3,v3,w3, gxk,gyk,gzk, Tk_m);
                u4=us+dt_sub*ud3; v4=vs+dt_sub*vd3; w4=ws+dt_sub*wd3;

                [ud4,vd4,wd4] = trans_dot_fn(pk,qk,rk, u4,v4,w4, gxk,gyk,gzk, Tk_m);

                us = us + dt_sub/6 * (ud1 + 2*ud2 + 2*ud3 + ud4);
                vs = vs + dt_sub/6 * (vd1 + 2*vd2 + 2*vd3 + vd4);
                ws = ws + dt_sub/6 * (wd1 + 2*wd2 + 2*wd3 + wd4);
            end
            u(k+1) = us; v(k+1) = vs; w(k+1) = ws;
        end
        res.u = u; res.v = v; res.w = w;
    end

    %% ============= Acelerômetro (modelo de sensor) =============
    % α calculado por derivada numérica de p, q, r (gradient).
    p_dot_sig = gradient(res.p, dt);
    q_dot_sig = gradient(res.q, dt);
    r_dot_sig = gradient(res.r, dt);

    [res.accX, res.accY, res.accZ] = accelerometer_model( ...
        res.p, res.q, res.r, ...
        res.u, res.v, res.w, ...
        T_total/m, ...
        p_dot_sig, q_dot_sig, r_dot_sig, ...
        r_imu);
end
