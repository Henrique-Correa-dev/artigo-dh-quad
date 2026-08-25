% prior_damping_moment.m — a priori de L_p, M_q, N_r direto da física
% =========================================================================
% Na forma de MOMENTO o valor a priori sai mais direto e mais forte do que na
% forma de taxa, por um motivo estrutural: a derivação física produz um MOMENTO
% em N·m·s, e a inércia só entrava depois, na divisão que gerava c_p. Tirando
% essa divisão, o a priori deixa de depender do tensor de inércia, que é o
% número mais incerto que temos (o CAD tem 1,072 kg contra 1,993 kg pesados).
%
% Consequência prática: este a priori vale para QUALQUER campanha que use os
% mesmos rotores e os mesmos braços, mesmo com massa, CG e inércia diferentes.
%
% MECANISMOS (todos em N·m·s, nenhum divide por J)
%   [A] influxo      L_p = 4·k_v·l_y²        M_q = 4·k_v·l_x²
%       k_v é a queda de empuxo por m/s de velocidade axial. MEDIDO no próprio
%       voo em 2_model/measure_kv.m: 0,226 ± 0,013 N/(m/s), contra 0,239 da
%       teoria de elemento de pá. Não é mais hipótese, é medida.
%   [B] momento de cubo   4·N_b·ρ·a·c_pá·Ω·R⁴/16      (rolagem e arfagem)
%   [C] força H      N_r = 4·k_h·(l_x² + l_y²)        (guinada)
%   [D] esteira      (2·D/w)·l²   com D o download medido pelo excesso de empuxo
%   [E] estrutura    ¼ρSb²·|C_lp|·V  etc., com os coeficientes MEDIDOS no voo
%                    de asa fixa, avaliados no p95 do envelope de treino
%
% Uso:  >> prior_damping_moment
% =========================================================================
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
rho = p.rho;  g = p.g;  m = p.m;

%% rotor (hélice 1045) — mesmas hipóteses do prior_damping
R = 0.127;  A = pi*R^2;  Nb = 2;  c_blade = 0.025;  a_lift = 5.0;  Cd0 = 0.03;
pitch_in = 4.5;
sigma = Nb*c_blade/(pi*R);
theta75 = atan(pitch_in/(0.75*2*R/0.0254*pi));
T0 = m*g/4;
rpm_h = interp1(p.bench.T_grams*1e-3*g, p.bench.RPM, T0, 'linear');
Om = rpm_h*2*pi/60;  OmR = Om*R;
vh = sqrt(T0/(2*rho*A));  lam = vh/OmR;
kap_teo = a_lift*sigma/(16*lam + a_lift*sigma);
kv_teo  = (T0/vh)*kap_teo;
k_h     = rho*A*OmR*(sigma/4)*(Cd0 + a_lift*theta75*lam);

% k_v MEDIDO no voo, se disponível
kv = kv_teo;  kv_se = NaN;  fonte = 'teoria (elemento de pá)';
fkv = fullfile(paths.outputs,'measure_kv.mat');
if exist(fkv,'file')
    K = load(fkv);  kv = K.kv_hat;  kv_se = K.se(2)/4;  fonte = 'MEDIDO no voo';
end

ly = mean([p.arms.Lx_r, p.arms.Lx_l]);      % braço lateral (rolagem)
lx = mean([p.arms.Ly_f, p.arms.Ly_r]);      % braço longitudinal (arfagem)
d2 = lx^2 + ly^2;

%% download e esteira
DL_FRAC = 0.078;  D_dl = DL_FRAC*m*g;  w_wake = 1.5*vh;  kD = 2*D_dl/w_wake;

%% aerodinâmica da estrutura, coeficientes medidos em asa fixa
S = p.wing.S;  b = p.wing.b;  cbar = p.wing.c;  V_ENV = 1.29;
Clp = -0.0769; Cmq = -5.540; Cnr = -0.2037;    % 11_fixed_wing (médias de 4 doublets)
fl = fullfile(paths.outputs,'fw_lateral.mat');  fo = fullfile(paths.outputs,'fw_longitudinal.mat');
if exist(fl,'file') && exist(fo,'file')
    T_ = load(fl);  L_ = load(fo);
    gi = @(S_,nm) mean(S_.Pid(strcmp(S_.pn,nm),:));
    Clp = gi(T_,'Clp');  Cmq = gi(L_,'Cmq');  Cnr = gi(T_,'Cnr');
end

%% balanço em N·m·s
B = zeros(4,3);
B(1,:) = [4*kv*ly^2,            4*kv*lx^2,             0];                     % influxo
B(2,:) = [4*Nb*rho*a_lift*c_blade*Om*R^4/16, 4*Nb*rho*a_lift*c_blade*Om*R^4/16, 0];  % cubo
B(3,:) = [kD*ly^2,              kD*lx^2,               0];                     % esteira
B(4,:) = [0.25*rho*S*b^2*abs(Clp), 0.25*rho*S*cbar^2*abs(Cmq), ...
          0.25*rho*S*b^2*abs(Cnr)] * V_ENV;                                    % estrutura
B(3,3) = 4*k_h*d2;                                                             % força H (guinada)

mech = {'[A] influxo (4·k_v·l²)','[B] momento de cubo','[C+D] esteira / força H','[E] estrutura ∝ V'};
ID  = [0.262265, 0.363162, 0.093578];    % identificado (moment_Jxz_2026 / moment_2026)

fprintf('\n  ROTOR: %.0f rpm | v_i = %.2f m/s | λ_i = %.4f\n', rpm_h, vh, lam);
fprintf('  k_v = %.4f N/(m/s)  [%s', kv, fonte);
if isfinite(kv_se), fprintf(', ± %.4f', kv_se); end
fprintf(']   k_h = %.4f\n', k_h);
fprintf('  braços: l_y = %.4f m (rolagem)   l_x = %.4f m (arfagem)\n', ly, lx);
fprintf('  C_lp = %.4f | C_mq = %.2f | C_nr = %.4f  (voo de asa fixa, em V = %.2f m/s)\n\n', Clp, Cmq, Cnr, V_ENV);
fprintf('  %-26s %12s %12s %12s\n', 'mecanismo [N·m·s]', 'L_p', 'M_q', 'N_r');
for i = 1:4, fprintf('  %-26s %12.5f %12.5f %12.5f\n', mech{i}, B(i,:)); end
fprintf('  %-26s %12.5f %12.5f %12.5f\n', 'SOMA (a priori)', sum(B,1));
fprintf('  %-26s %12.5f %12.5f %12.5f\n', 'SÓ influxo', B(1,:));
fprintf('  %-26s %12.5f %12.5f %12.5f\n', 'identificado (OEM)', ID);
fprintf('  %-26s %11.0f%% %11.0f%% %11.0f%%\n', 'soma / identificado', 100*sum(B,1)./ID);
fprintf('  %-26s %11.0f%% %11.0f%% %11.0f%%\n', 'influxo / identificado', 100*B(1,:)./ID);

prior_m = struct('L_p_influxo',B(1,1), 'M_q_influxo',B(1,2), 'N_r_H',B(3,3), ...
                 'L_p_total',sum(B(:,1)), 'M_q_total',sum(B(:,2)), 'N_r_total',sum(B(:,3)), ...
                 'k_v',kv, 'k_h',k_h, 'B',B, 'mech',{mech});
save(fullfile(paths.outputs,'prior_damping_moment.mat'), 'prior_m');
fprintf('\n  Salvo em outputs/prior_damping_moment.mat\n');
fprintf('  P0 OFICIAL (só influxo, k_v medido; força H em guinada):\n');
fprintf('     p.P0_J(13:15) = [%.4f; %.4f; %.4f]\n', B(1,1), B(1,2), B(3,3));
fprintf('  (asa é termo separado; cubo e esteira ficam só como estimativa no balanço)\n');
