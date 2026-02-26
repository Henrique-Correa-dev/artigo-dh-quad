function error_vector = oem_cost_function(P_scaled, P0_scale, weights, time_exp, pqr_exp, ...
    pwm_time, pwm_data, func_T_ref, func_Q_ref, seg_len)
% OEM_COST_FUNCTION  Multiple-shooting OEM cost for VTOL identification.
%
%   Uses short integration windows (multiple shooting) so the ODE never
%   diverges far from the data. Each window starts from the experimental
%   state at that time, providing well-conditioned gradient information.
%
%   seg_len: number of time steps per shooting segment (e.g., 10 => 1 s at dt=0.1)

    N = length(time_exp);
    P_real = P_scaled(:) .* P0_scale(:);
    sqrt_w = sqrt(weights(:));

    error_p = zeros(N, 1);
    error_q = zeros(N, 1);
    error_r = zeros(N, 1);

    ode_opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
    ode_func = @(t, y) vtol_dynamics(t, y, P_real, pwm_time, pwm_data, func_T_ref, func_Q_ref);

    seg_starts = 1:seg_len:N;
    for k = 1:length(seg_starts)
        i1 = seg_starts(k);
        i2 = min(i1 + seg_len, N);
        idx = i1:i2;

        ic_seg = pqr_exp(i1, :)';
        t_seg  = time_exp(idx);

        try
            [t_out, y_out] = ode45(ode_func, t_seg(:), ic_seg, ode_opts);
        catch
            error_p(idx) = 1e4;
            error_q(idx) = 1e4;
            error_r(idx) = 1e4;
            continue;
        end

        if length(t_out) < 2 || any(~isfinite(y_out(:)))
            error_p(idx) = 1e4;
            error_q(idx) = 1e4;
            error_r(idx) = 1e4;
            continue;
        end

        p_s = interp1(t_out, y_out(:,1), t_seg(:), 'linear', 'extrap');
        q_s = interp1(t_out, y_out(:,2), t_seg(:), 'linear', 'extrap');
        r_s = interp1(t_out, y_out(:,3), t_seg(:), 'linear', 'extrap');

        p_s(~isfinite(p_s)) = pqr_exp(idx(~isfinite(p_s)), 1);
        q_s(~isfinite(q_s)) = pqr_exp(idx(~isfinite(q_s)), 2);
        r_s(~isfinite(r_s)) = pqr_exp(idx(~isfinite(r_s)), 3);

        error_p(idx) = p_s(:) - pqr_exp(idx, 1);
        error_q(idx) = q_s(:) - pqr_exp(idx, 2);
        error_r(idx) = r_s(:) - pqr_exp(idx, 3);
    end

    error_vector = [sqrt_w(1)*error_p; sqrt_w(2)*error_q; sqrt_w(3)*error_r];
end
