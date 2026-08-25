% damping_vs_speed.m — como o amortecimento angular muda do pairado ao cruzeiro
% =========================================================================
% Curva L_p(V), M_q(V), N_r(V) com as parcelas empilhadas, do pairado até a
% velocidade de estol, mostrando o que decai e o que cresce.
%
% PARCELAS E LEIS DE ESCALA
%   rotor (influxo + cubo)   ∝ T/mg         k_v ∝ √T·(T/v_h) e o cubo ∝ Ω,
%                                           no conjunto ≈ linear em T; cai conforme a asa assume
%   efetivo (não explicado)  ∝ T/mg         atribuído à dinâmica de motor; existe enquanto o
%                                           amortecimento for de rotor; vai junto com T
%   asa                      ∝ V            ¼ρSb²·|C_lp|·V, C_lp MEDIDO nos doublets
%
% A lei de transferência de sustentação é o balanço vertical a incidência
% constante,  T/mg = 1 − (V/V_stall)²,  que Hu et al. (2024) mediram em voo
% para um VTOL composto de 5,5 kg entre 0 e 9 m/s.
%
% O que está MEDIDO: os dois extremos. Pairado (V ≈ 1,3 m/s, L_p identificado)
% e cruzeiro (V = 17 m/s, C_lp dos doublets). O meio, de 3 m/s ao estol, é
% interpolação física com as leis acima. Sem dado nessa faixa.
%
% Uso:  >> damping_vs_speed
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
rho = p.rho; S = p.wing.S; b = p.wing.b; c = p.wing.c; g = p.g; m = p.m;
V_STALL = 9.28;   V_ENV = 1.29;  V_VAL = 2.21;  V_FW = 17.0;

% parcelas em N·m·s (prior_damping_moment sem a esteira, e o identificado oficial)
PR = load(fullfile(paths.outputs,'prior_damping_moment.mat')).prior_m;
ID = load(fullfile(paths.outputs,'runs','oficial_2026_final','P_identified.mat')).P_final(:);
Lp_id = ID(13); Mq_id = ID(14); Nr_id = ID(15);
Clp = abs(ID(17)); Cmq = abs(ID(19)); Cnr = abs(ID(21));      % medidos em asa fixa
rot = [PR.B(1,1), PR.B(1,2), PR.B(3,3)];   % só influxo (+ força H em guinada)  % influxo+cubo (N_r: força H)
kasa = [0.25*rho*S*b^2*Clp, 0.25*rho*S*c^2*Cmq, 0.25*rho*S*b^2*Cnr];   % por m/s
asa_env = kasa*V_ENV;
% No modelo oficial_fw a asa entra como termo SEPARADO (¼ρSb²·C_lp·V_a·p), então
% o L_p identificado NÃO a contém. O efetivo é o identificado menos só o rotor.
efet = [Lp_id Mq_id Nr_id] - rot;                             % o que a identificação acrescenta

V = linspace(0, V_STALL, 300)';
Tfrac = max(1 - (V/V_STALL).^2, 0);
nomes = {'L_p (rolagem)','M_q (arfagem)','N_r (guinada)'};
cor = [0.20 0.40 0.65; 0.55 0.55 0.55; 0.45 0.70 0.35];

fprintf('\n  Parcelas em V = %.2f m/s [N·m·s]      rotor    efetivo   = L_p id.    asa (separada)   total\n', V_ENV);
for k = 1:3
    fprintf('  %-16s %8.4f %8.4f %10.4f %12.4f %12.4f\n', nomes{k}, rot(k), efet(k), rot(k)+efet(k), asa_env(k), rot(k)+efet(k)+asa_env(k));
end
fprintf('\n  L_p(V) nos marcos:\n');
for Vk = [V_ENV V_VAL 3 6.27 V_STALL]
    tf = max(1-(Vk/V_STALL)^2,0);
    fprintf('    V = %4.2f m/s:  rotor %.4f  efetivo %.4f  asa %.4f  TOTAL %.4f   (T/mg %.2f)\n', ...
        Vk, rot(1)*tf, efet(1)*tf, kasa(1)*Vk, (rot(1)+efet(1))*tf + kasa(1)*Vk, tf);
end
fprintf('    V = %4.1f m/s (asa fixa, medido): asa %.4f\n', V_FW, kasa(1)*V_FW);

f = figure('Position',[40 40 1350 470],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile; hold on; grid on;
    y_rot  = rot(k)*Tfrac;
    y_efet = efet(k)*Tfrac;
    y_asa  = kasa(k)*V;
    area(V, [y_rot, y_efet, y_asa], 'LineStyle','none');
    colororder(gca, cor);
    % marcos
    xregion(0, V_VAL, 'FaceColor',[0 0.45 0.7], 'FaceAlpha',0.10);
    xregion(3, V_STALL, 'FaceColor',[0.9 0.9 0.9], 'FaceAlpha',0.5);
    % ponto: L_p identificado MAIS a asa em V_ENV (= amortecimento total no pairado)
    plot(V_ENV, (Lp_id*(k==1)+Mq_id*(k==2)+Nr_id*(k==3)) + kasa(k)*V_ENV, 'o', 'MarkerSize',9, ...
        'MarkerFaceColor',[0.85 0.37 0.01], 'MarkerEdgeColor','k');
    xline(6.27, 'k:', 'V = v_i', 'LabelVerticalAlignment','bottom', 'FontSize',8);
    xline(V_STALL, 'k--', 'V_{stall}', 'LabelVerticalAlignment','bottom', 'FontSize',8);
    xlabel('V [m/s]'); ylabel([nomes{k}(1:3) ' [N·m·s]']); xlim([0 V_STALL*1.02]);
    title(nomes{k});
    if k == 1
        legend({'rotor (influxo, k_v medido)  \propto T', 'efetivo (motor)  \propto T', ...
                'asa, C medido em voo  \propto V', 'envelope ensaiado', 'sem dado (interpolação)', ...
                'L_p identificado + asa, em pairado'}, 'Location','northeast','FontSize',7.5);
    end
end
title(tl, 'Amortecimento angular do pairado ao estol: o rotor decai com o empuxo, a asa cresce com a velocidade', ...
    'FontWeight','bold');
fn = fullfile(paths.images,'damping_vs_speed.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',150);
copyfile(fn, fullfile(getenv('HOME'),'Desktop','DH_modelo_oficial','damping_vs_speed.png'));
fprintf('\n  Figura: %s  (cópia na Mesa)\n', fn);
