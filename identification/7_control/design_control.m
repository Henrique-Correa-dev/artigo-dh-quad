% design_control.m — Projeto do controle CASCATA (PID) sobre o modelo LINEAR
%
% Pega o modelo linear [T,Mx,My,Mz]→[p q r φ θ ψ u v w] (linearize.m),
% aumenta com altitude h (ḣ = -w), e projeta cada malha SISO por ALOCAÇÃO
% DE POLOS (ωn, ζ escolhidos por banda). Estrutura:
%
%   altitude h  → T   (duplo integrador (1/m)/s²)      : PID
%   vel.forward → θc  (u̇ = -g·θ) → θ → My             : PI externo + PD interno
%   roll φ=0    → Mx  (φ̈ = -Dp·φ̇ + bφ·Mx)             : PD
%   heading ψ   → Mz  (ψ̈ = -Dr·ψ̇ + bψ·Mz)             : PD
%
% Separação de banda: atitude (interna) RÁPIDA; velocidade/altitude (externa)
% LENTA (~4×). Salva control_gains.mat.  Uso: >> design_control

addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ===================== CARREGAR MODELO LINEAR =====================
lm = load(fullfile(paths.outputs, 'linear_model.mat'));
A = lm.A;  B = lm.B;
proj = parameters();
m = proj.m;  g = proj.g;

% Coeficientes físicos extraídos do modelo (SISO desacoplado perto do hover)
a_p = -A(1,1);   b_ph = B(1,2);    % roll:  φ̈ = -a_p·φ̇ + b_ph·Mx
a_q = -A(2,2);   b_th = B(2,3);    % pitch: θ̈ = -a_q·θ̇ + b_th·My
a_r = -A(3,3);   b_ps = B(3,4);    % yaw:   ψ̈ = -a_r·ψ̇ + b_ps·Mz
b_w = B(9,1);                      % ẇ = b_w·T  (b_w = -1/m)

fprintf('==========================================================\n');
fprintf('  PROJETO DE CONTROLE CASCATA (modelo linear)\n');
fprintf('==========================================================\n');
fprintf('  Coef físicos:  roll b=%.3f a=%.3f | pitch b=%.3f a=%.3f | yaw b=%.3f a=%.3f\n', ...
    b_ph, a_p, b_th, a_q, b_ps, a_r);
fprintf('  Vertical: ẇ = %.4f·T  (esperado -1/m = %.4f)\n', b_w, -1/m);

%% ===================== ESPECIFICAÇÃO (ωn, ζ por malha) =====================
% Interna (atitude) rápida; externa (vel/altitude) ~4× mais lenta.
wn_th = 6.0;  z_th = 0.8;    % pitch  (interno)
wn_ph = 6.0;  z_ph = 0.8;    % roll   (interno)
wn_ps = 2.5;  z_ps = 0.8;    % yaw
wn_u  = 1.5;  z_u  = 0.9;    % vel forward (externo)
wn_h  = 2.0;  z_h  = 0.9;    % altitude

%% ===================== ALOCAÇÃO DE POLOS =====================
% Atitude (2ª ordem ang/torque = b/(s(s+a))): M = Kp(ang_c-ang) - Kd·rate
%   malha fechada: s² + (a + b·Kd)s + b·Kp  →  Kp=ωn²/b, Kd=(2ζωn - a)/b
K.Kp_th = wn_th^2 / b_th;      K.Kd_th = (2*z_th*wn_th - a_q) / b_th;
K.Kp_ph = wn_ph^2 / b_ph;      K.Kd_ph = (2*z_ph*wn_ph - a_p) / b_ph;
K.Kp_ps = wn_ps^2 / b_ps;      K.Kd_ps = (2*z_ps*wn_ps - a_r) / b_ps;

% Vel forward (externo): com θ interno rápido, u̇ ≈ -g·θ. Lei θc=-(Kp·e+Ki·∫)/g
%   ⇒ u̇ = Kp·e_u + Ki·∫e_u  →  s² + Kp·s + Ki  →  Kp=2ζωn, Ki=ωn²
K.Kp_u = 2*z_u*wn_u;           K.Ki_u = wn_u^2;

% Altitude (duplo integrador ḧ = (1/m)·δT): PID
%   ḧ + (Kd/m)ḣ + (Kp/m)h = (Kp/m)h_sp  →  Kp=m·ωn², Kd=m·2ζωn
K.Kp_h = m*wn_h^2;             K.Kd_h = m*2*z_h*wn_h;
K.Ki_h = 0.5*m*wn_h^3*0;       % I de altitude (FF de gravidade já zera bias) — 0 por ora

% Limites / anti-windup
K.theta_max = deg2rad(20);     % tilt máximo (vel forward)
K.int_u_max = 5.0;             % anti-windup vel
K.int_h_max = 5.0;             % anti-windup altitude
K.g = g;  K.m = m;

%% ===================== POLOS DE MALHA FECHADA (checagem) =====================
fprintf('\n  --- Ganhos projetados ---\n');
fprintf('    Pitch   : Kp=%.3f Kd=%.3f   (ωn=%.1f ζ=%.2f)\n', K.Kp_th, K.Kd_th, wn_th, z_th);
fprintf('    Roll    : Kp=%.3f Kd=%.3f   (ωn=%.1f ζ=%.2f)\n', K.Kp_ph, K.Kd_ph, wn_ph, z_ph);
fprintf('    Yaw     : Kp=%.3f Kd=%.3f   (ωn=%.1f ζ=%.2f)\n', K.Kp_ps, K.Kd_ps, wn_ps, z_ps);
fprintf('    Vel fwd : Kp=%.3f Ki=%.3f   (ωn=%.1f ζ=%.2f)\n', K.Kp_u, K.Ki_u, wn_u, z_u);
fprintf('    Altitude: Kp=%.3f Kd=%.3f   (ωn=%.1f ζ=%.2f)\n', K.Kp_h, K.Kd_h, wn_h, z_h);

show_poles('Pitch   ', [1, a_q + b_th*K.Kd_th, b_th*K.Kp_th]);
show_poles('Roll    ', [1, a_p + b_ph*K.Kd_ph, b_ph*K.Kp_ph]);
show_poles('Yaw     ', [1, a_r + b_ps*K.Kd_ps, b_ps*K.Kp_ps]);
show_poles('Vel fwd ', [1, K.Kp_u, K.Ki_u]);
show_poles('Altitude', [1, K.Kd_h/m, K.Kp_h/m]);

%% ===================== SALVAR =====================
gains_path = fullfile(paths.control, 'control_gains.mat');
save(gains_path, 'K');
fprintf('\n  Ganhos salvos: %s\n', gains_path);
fprintf('==========================================================\n');


%% ===================== HELPER =====================
function show_poles(name, den)
    pr = roots(den);
    s = sprintf('%+.3f%+.3fj  ', [real(pr), imag(pr)]');
    if all(real(pr) < 0)
        tag = 'estavel';
    else
        tag = 'INSTAVEL';
    end
    fprintf('    polos %s: %s (%s)\n', name, s, tag);
end
