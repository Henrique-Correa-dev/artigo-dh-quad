% damping_terms_validation.m — L_p·p contra o termo da asa, no trecho de validação
% =========================================================================
% Duas figuras:
%  (1) L_p, M_q, N_r total contra V até 3 m/s (envelope real), com rotor,
%      efetivo e asa empilhados. Nada de extrapolação até o estol.
%  (2) Série temporal no trecho de validação 610–630 s: o momento de
%      amortecimento de rotor L_p·p e o momento da asa ¼ρSb²·C_lp·V_a·p,
%      calculados com p e V_a MEDIDOS/ESTIMADOS no trecho, lado a lado. É o que
%      o modelo está de fato somando em cada instante.
% Uso:  >> damping_terms_validation
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
rho = p.rho; S = p.wing.S; b = p.wing.b; c = p.wing.c;
T_WIN = [610 630];  V_ENV = 1.29;  V_VAL = 2.21;

PR = load(fullfile(paths.outputs,'prior_damping_moment.mat')).prior_m;
ID = load(fullfile(paths.outputs,'runs','oficial_2026_final','P_identified.mat')).P_final(:);
Lid = ID(13:15)';
Cw  = abs([ID(17) ID(19) ID(21)]);                                % C_lp, C_mq, C_nr medidos
rot = [PR.B(1,1), PR.B(1,2), PR.B(3,3)];   % só influxo (+ força H em guinada)
efe = Lid - rot;
kasa = [0.25*rho*S*b^2*Cw(1), 0.25*rho*S*c^2*Cw(2), 0.25*rho*S*b^2*Cw(3)];  % N·m·s por m/s
nomes = {'rolagem','arfagem','guinada'};  sym = {'L_p','M_q','N_r'};  tx = {'p','q','r'};

%% ---------------- figura 1: até 3 m/s
V = linspace(0, 3, 200)';
f1 = figure('Position',[40 40 1350 430],'Color','w'); try, f1.Theme='light'; catch, end
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile; hold on; grid on;
    area(V, [rot(k)*ones(size(V)), efe(k)*ones(size(V)), kasa(k)*V], 'LineStyle','none');
    colororder(gca, [0.20 0.40 0.65; 0.55 0.55 0.55; 0.45 0.70 0.35]);
    xline(V_ENV, '--', sprintf('p95 treino %.2f', V_ENV), 'Color',[0 0.45 0.7], 'LabelVerticalAlignment','bottom','FontSize',8);
    xline(V_VAL, '--', sprintf('p95 validação %.2f', V_VAL), 'Color',[0.85 0.37 0.01], 'LabelVerticalAlignment','bottom','FontSize',8);
    fr = @(Vk) 100*kasa(k)*Vk/(Lid(k)+kasa(k)*Vk);
    text(0.03, 0.10, sprintf('asa / total:  %.1f%% em %.2f m/s,  %.1f%% em %.2f m/s', fr(V_ENV), V_ENV, fr(V_VAL), V_VAL), ...
        'Units','normalized','FontSize',8.5,'BackgroundColor',[1 1 1 0.8]);
    xlabel('V [m/s]'); ylabel([sym{k} ' total [N·m·s]']); xlim([0 3]); title(nomes{k});
    if k == 1
        legend({'rotor (influxo, k_v medido)','efetivo (identificado − rotor)','asa, C medido em voo, \propto V'}, ...
            'Location','northwest','FontSize',8);
    end
end
title(tl, 'Amortecimento angular no envelope ensaiado (até 3 m/s): parcela da asa contra o L_p identificado', 'FontWeight','bold');
fn1 = fullfile(paths.images,'damping_terms_envelope.png');
exportgraphics(f1, fn1, 'BackgroundColor','white','Resolution',150);

%% ---------------- figura 2: séries no trecho de validação
L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU) min(L.time_ATT) min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU) max(L.time_ATT) max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';  VE = estimate_velocity(L, tg);
idx = tg>=T_WIN(1) & tg<=T_WIN(2);  t = tg(idx);
Va = VE.V(idx); Va(~isfinite(Va)) = 0;
ip = @(tt,xx) interp1(tt, xx, t, 'linear');
om = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
% momentos de controle dos rotores, pela mesma cadeia de motor e k_T, k_Q identificados
W4 = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
[rpm, fT, fQ] = motor_chain(t, W4);
Tmr = (ID(5:8)') .* fT(rpm);   Qmr = (ID(9:12)') .* fQ(rpm);      % N×4
Lx_r = p.arms.Lx_r; Lx_l = p.arms.Lx_l; Ly_f = p.arms.Ly_f; Ly_r = p.arms.Ly_r;
Mx_ctl = -(Lx_r*Tmr(:,1) - Lx_l*Tmr(:,2) - Lx_l*Tmr(:,3) + Lx_r*Tmr(:,4));
My_ctl =   Ly_f*Tmr(:,1) - Ly_r*Tmr(:,2) + Ly_f*Tmr(:,3) - Ly_r*Tmr(:,4);
Mz_ctl =   Qmr(:,1) + Qmr(:,2) - Qmr(:,3) - Qmr(:,4);
M_ctl = [Mx_ctl - mean(Mx_ctl), My_ctl - mean(My_ctl), Mz_ctl - mean(Mz_ctl)];   % sem o trim

f2 = figure('Position',[40 40 1250 900],'Color','w'); try, f2.Theme='light'; catch, end
tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
fprintf('\n  Trecho %d–%d s: RMS dos momentos de amortecimento [N·m]\n', T_WIN);
fprintf('  %-9s %10s %12s %12s %10s %10s\n', 'eixo', 'M_rot', 'L·taxa (id.)', 'asa', 'L/M_rot', 'asa/M_rot');
for k = 1:3
    M_id  = Lid(k)*om(:,k);                 % o que o L_p identificado produz
    M_asa = kasa(k)*Va.*om(:,k);            % o que a asa produz
    M_rot = rot(k)*om(:,k);                 % a parcela física de rotor dentro do L_p
    nexttile; hold on; grid on;
    plot(t, M_ctl(:,k), '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.2);
    plot(t, M_id,  '-', 'Color',[0.30 0.30 0.30], 'LineWidth',1.6);
    plot(t, M_asa, '-', 'Color',[0.45 0.70 0.35], 'LineWidth',1.8);
    ylabel(sprintf('%s [N·m]', nomes{k})); xlim(T_WIN);
    r = rms(M_asa)/max(rms(M_id),eps);  rc = rms(M_ctl(:,k));
    text(0.006, 0.94, sprintf('RMS:  M_rot = %.4f    %s·%s = %.4f (%.0f%% de M_rot)    asa = %.4f (%.1f%% de %s·%s, %.1f%% de M_rot)', ...
        rc, sym{k}, tx{k}, rms(M_id), 100*rms(M_id)/rc, rms(M_asa), 100*r, sym{k}, tx{k}, 100*rms(M_asa)/rc), ...
        'Units','normalized','VerticalAlignment','top','FontSize',8.5,'FontName','Menlo','BackgroundColor',[1 1 1 0.8]);
    if k == 1
        legend({'M_{rot}  (controle dos rotores)', sprintf('%s·p   (amortecimento identificado)', sym{1}), ...
                '¼\rhoSb²·C_{lp}·V_a·p   (asa, C_{lp} medido em voo)'}, 'Location','southeast','FontSize',8.5);
    end
    fprintf('  %-9s %10.4f %12.4f %12.4f %9.1f%% %9.1f%%\n', nomes{k}, rc, rms(M_id), rms(M_asa), 100*rms(M_id)/rc, 100*rms(M_asa)/rc);
end
nexttile; hold on; grid on;
plot(t, Va, 'k-', 'LineWidth',1.3); ylabel('V_a estimada [m/s]'); xlabel('t [s]'); xlim(T_WIN);
title(tl, sprintf('Trecho de validação %d–%d s: controle dos rotores, amortecimento identificado e asa, com p, q, r medidos e V_a estimada', T_WIN), 'FontWeight','bold');
fn2 = fullfile(paths.images,'damping_terms_validation.png');
exportgraphics(f2, fn2, 'BackgroundColor','white','Resolution',140);
dd = fullfile(getenv('HOME'),'Desktop','DH_modelo_oficial');
copyfile(fn1, fullfile(dd,'damping_terms_envelope.png'));  copyfile(fn2, fullfile(dd,'damping_terms_validation.png'));
fprintf('\n  Figuras: %s\n           %s  (cópias na Mesa)\n', fn1, fn2);
