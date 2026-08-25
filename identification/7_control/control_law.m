function [u_cmd, cs, dbg] = control_law(x, h, sp, K, cs, dt)
%CONTROL_LAW  Lei de controle em CASCATA (PID) — saída u = [T, Mx, My, Mz].
%
% Estrutura (forward via pitch, roll regulado a zero):
%   ── Altitude h        → T   (PD + I + feedforward de gravidade m·g)
%   ── Vel. forward u    → θc  (PI externo) → θ → My  (PD interno de atitude)
%   ── Roll φ = 0        → Mx  (PD: segura roll nivelado)
%   ── Heading ψ         → Mz  (PD com wrap de ±180°)
%
% É a MESMA lei usada na simulação LINEAR e NÃO-LINEAR. No linear u entra
% direto na planta; no não-linear u passa pela ponte forces_to_pwm (alocação).
%
% ENTRADAS
%   x   [9x1]  estado [p q r φ θ ψ u v w]  (rad, rad/s, m/s)
%   h   escalar  altitude atual (m)        (estado aumentado, ḣ = -w)
%   sp  struct  setpoints: .h, .u (vel forward), .psi (heading, rad)
%   K   struct  ganhos (de design_control.m / control_gains.mat)
%   cs  struct  estado dos integradores: .int_h, .int_u  (passa e retorna)
%   dt  escalar  passo de tempo (s)
%
% SAÍDAS
%   u_cmd [4x1]  [T; Mx; My; Mz]
%   cs    struct integradores atualizados
%   dbg   struct sinais internos (θc, e_h, e_u, e_psi) p/ diagnóstico

    p = x(1); q = x(2); r = x(3);
    phi = x(4); theta = x(5); psi = x(6);
    u = x(7); w = x(9);

    g = K.g;  m = K.m;

    %% ---- Malha de ALTITUDE (h → T) ----
    e_h = sp.h - h;
    cs.int_h = clamp(cs.int_h + e_h*dt, -K.int_h_max, K.int_h_max);   % anti-windup
    % ḣ = -w  ⇒  termo derivativo -Kd·ḣ = +Kd·w
    dT = K.Kp_h*e_h + K.Ki_h*cs.int_h + K.Kd_h*w;
    T  = m*g + dT;
    T  = max(T, 0);                       % empuxo não-negativo

    %% ---- Malha de VEL. FORWARD (u → θc) ----  (externa, lenta)
    e_u = sp.u - u;
    int_u_try = clamp(cs.int_u + e_u*dt, -K.int_u_max, K.int_u_max);
    % u̇ = -g·θ  ⇒  para acelerar forward (e_u>0) precisa de θ<0 (nariz baixo)
    theta_c_raw = -(K.Kp_u*e_u + K.Ki_u*int_u_try) / g;
    theta_c = clamp(theta_c_raw, -K.theta_max, K.theta_max);   % limite de tilt
    % anti-windup: só acumula se θc NÃO está saturado (integração condicional)
    if theta_c == theta_c_raw
        cs.int_u = int_u_try;
    end   % senão congela o integrador (evita windup durante saturação)

    %% ---- Malha de PITCH (θc → My) ----  (interna, rápida)
    My = K.Kp_th*(theta_c - theta) - K.Kd_th*q;

    %% ---- Malha de ROLL (φ=0 → Mx) ----  (regula roll nivelado)
    phi_c = 0;
    Mx = K.Kp_ph*(phi_c - phi) - K.Kd_ph*p;

    %% ---- Malha de HEADING (ψ → Mz) ----
    e_psi = wrap_pi(sp.psi - psi);
    Mz = K.Kp_ps*e_psi - K.Kd_ps*r;

    u_cmd = [T; Mx; My; Mz];
    dbg = struct('theta_c', theta_c, 'e_h', e_h, 'e_u', e_u, 'e_psi', e_psi, 'dT', dT);
end

function y = clamp(x, lo, hi)
    y = min(max(x, lo), hi);
end

function a = wrap_pi(a)
    a = mod(a + pi, 2*pi) - pi;
end
