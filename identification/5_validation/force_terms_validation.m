% force_terms_validation.m — as parcelas que compõem F_x, F_y, F_z
% =========================================================================
% Análogo do damping_terms_validation para as forças. O modelo tem:
%   F = F_rot + F_grav + F_asa
%   F_rot  = [0; 0; −T]                          empuxo dos rotores (k_T identificado)
%   F_asa  = ½ρS·[ −C_D0·V² ;  C_Yβ·V·v ;  −(C_L0·V² + C_Lα·V·w) ]
%            C_L0 = 0,34, C_Lα = 2,94, C_Yβ = −0,196 (asa fixa), C_D0 = 0,05
% A gravidade é estado (atitude), não parâmetro, e o acelerômetro não a vê.
% Então a comparação de escala é F_asa contra F_rot, que é o que o sensor mede.
%
% Figura 1: F_asa por eixo contra V até 3 m/s, com o empuxo de pairado (mg)
%           como referência, e a fração em cada marco.
% Figura 2: trecho de validação 610–630 s, F_rot e F_asa instante a instante.
% Uso:  >> force_terms_validation
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
rho = p.rho; S = p.wing.S; m = p.m; g = p.g;
T_WIN = [610 630];  V_ENV = 1.29;  V_VAL = 2.21;
CD0 = 0.05; CYb = 0.196; CL0 = 0.34; CLa = 2.94;     % módulos
qS = 0.5*rho*S;
ID = load(fullfile(paths.outputs,'runs','oficial_2026_final','P_identified.mat')).P_final(:);

%% ---------------- figura 1: até 3 m/s, parcelas da força da asa
V = linspace(0, 3, 200)';
% em w e v típicos do envelope (mediana do módulo): w ~ 0,47, v ~ 0,27 (janela 610-630)
w_typ = 0.47; v_typ = 0.27;
Fx = qS*CD0*V.^2;
Fy = qS*CYb*V*v_typ;
Fz0 = qS*CL0*V.^2;  Fza = qS*CLa*V*w_typ;
f1 = figure('Position',[40 40 1350 430],'Color','w'); try, f1.Theme='light'; catch, end
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on;
area(V, Fx, 'FaceColor',[0.45 0.70 0.35], 'LineStyle','none');
xline(V_ENV,'--','Color',[0 0.45 0.7]); xline(V_VAL,'--','Color',[0.85 0.37 0.01]);
xlabel('V [m/s]'); ylabel('|F_x| asa [N]'); title('longitudinal');
text(0.03,0.92, sprintf('½\\rhoS·C_{D0}·V²\n%.1f%% de mg em %.2f m/s, %.1f%% em %.2f', ...
    100*qS*CD0*V_ENV^2/(m*g), V_ENV, 100*qS*CD0*V_VAL^2/(m*g), V_VAL), 'Units','normalized','VerticalAlignment','top','FontSize',8.5,'BackgroundColor',[1 1 1 .8]);
nexttile; hold on; grid on;
area(V, Fy, 'FaceColor',[0.45 0.70 0.35], 'LineStyle','none');
xline(V_ENV,'--','Color',[0 0.45 0.7]); xline(V_VAL,'--','Color',[0.85 0.37 0.01]);
xlabel('V [m/s]'); ylabel('|F_y| asa [N]'); title(sprintf('lateral (v = %.2f m/s típico)', v_typ));
text(0.03,0.92, sprintf('½\\rhoS·|C_{Y\\beta}|·V·v\n%.1f%% de mg em %.2f m/s, %.1f%% em %.2f', ...
    100*qS*CYb*V_ENV*v_typ/(m*g), V_ENV, 100*qS*CYb*V_VAL*v_typ/(m*g), V_VAL), 'Units','normalized','VerticalAlignment','top','FontSize',8.5,'BackgroundColor',[1 1 1 .8]);
nexttile; hold on; grid on;
area(V, [Fz0, Fza], 'LineStyle','none'); colororder(gca, [0.30 0.55 0.25; 0.60 0.80 0.45]);
xline(V_ENV,'--','Color',[0 0.45 0.7]); xline(V_VAL,'--','Color',[0.85 0.37 0.01]);
xlabel('V [m/s]'); ylabel('|F_z| asa [N]'); title(sprintf('vertical (w = %.2f m/s típico)', w_typ));
legend({'½\rhoS·C_{L0}·V²','½\rhoS·C_{L\alpha}·V·w'},'Location','northwest','FontSize',8);
text(0.03,0.72, sprintf('total %.1f%% de mg em %.2f m/s, %.1f%% em %.2f\n(mg = %.1f N)', ...
    100*(qS*CL0*V_ENV^2+qS*CLa*V_ENV*w_typ)/(m*g), V_ENV, 100*(qS*CL0*V_VAL^2+qS*CLa*V_VAL*w_typ)/(m*g), V_VAL, m*g), ...
    'Units','normalized','VerticalAlignment','top','FontSize',8.5,'BackgroundColor',[1 1 1 .8]);
title(tl, 'Força aerodinâmica da asa no envelope ensaiado (até 3 m/s), coeficientes medidos em voo, contra o peso', 'FontWeight','bold');
fn1 = fullfile(paths.images,'force_terms_envelope.png');
exportgraphics(f1, fn1, 'BackgroundColor','white','Resolution',150);

%% ---------------- figura 2: trecho de validação
L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU) min(L.time_ATT) min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU) max(L.time_ATT) max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';  VE = estimate_velocity(L, tg);
idx = tg>=T_WIN(1) & tg<=T_WIN(2);  t = tg(idx);
Va = VE.V(idx); va = VE.v(idx); wa = VE.w(idx);
Va(~isfinite(Va))=0; va(~isfinite(va))=0; wa(~isfinite(wa))=0;
ip = @(tt,xx) interp1(tt, xx, t, 'linear');
W4 = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
[rpm, fT] = motor_chain(t, W4);
T = sum((ID(5:8)') .* fT(rpm), 2);
Fasa = [ -qS*CD0*Va.^2,  -qS*CYb*Va.*va,  -qS*(CL0*Va.^2 + CLa*Va.*wa) ];
Frot = [ zeros(size(t)), zeros(size(t)), -T ];

f2 = figure('Position',[40 40 1250 900],'Color','w'); try, f2.Theme='light'; catch, end
tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
nomes = {'F_x','F_y','F_z'};
fprintf('\n  Trecho %d–%d s [N]:   RMS F_rot (sem trim)   RMS F_asa    asa/rot    asa/mg\n', T_WIN);
for k = 1:3
    nexttile; hold on; grid on;
    if k == 3
        plot(t, Frot(:,k) - mean(Frot(:,k)), '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.2);
        lg = {'−T + T̄  (empuxo sem o trim)', 'F_z asa'};
        rc = rms(Frot(:,k) - mean(Frot(:,k)));
    else
        lg = {'F asa'}; rc = NaN;
    end
    plot(t, Fasa(:,k), '-', 'Color',[0.45 0.70 0.35], 'LineWidth',1.8);
    ylabel([nomes{k} ' [N]']); xlim(T_WIN);
    if k == 3
        legend(lg, 'Location','southeast','FontSize',8.5);
        text(0.006,0.94, sprintf('RMS:  empuxo (sem trim) = %.3f N    asa = %.3f N (%.1f%% do empuxo variável, %.2f%% de mg)', ...
            rc, rms(Fasa(:,k)), 100*rms(Fasa(:,k))/rc, 100*rms(Fasa(:,k))/(m*g)), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8.5,'FontName','Menlo','BackgroundColor',[1 1 1 .8]);
    else
        text(0.006,0.94, sprintf('RMS:  asa = %.3f N  (%.2f%% de mg)   |   não há força de rotor neste eixo (só acoplamento e gravidade)', ...
            rms(Fasa(:,k)), 100*rms(Fasa(:,k))/(m*g)), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8.5,'FontName','Menlo','BackgroundColor',[1 1 1 .8]);
    end
    fprintf('  %-4s %18.3f %12.3f %9s %9.2f%%\n', nomes{k}, rc, rms(Fasa(:,k)), ...
        ternary(isfinite(rc), sprintf('%.1f%%',100*rms(Fasa(:,k))/rc), '-'), 100*rms(Fasa(:,k))/(m*g));
end
nexttile; hold on; grid on;
plot(t, Va, 'k-', 'LineWidth',1.3); plot(t, va, '-', 'Color',[0 0.45 0.7]); plot(t, wa, '-', 'Color',[0.85 0.37 0.01]);
legend({'V_a','v','w'},'Location','northeast','Orientation','horizontal','FontSize',8);
ylabel('velocidade estimada [m/s]'); xlabel('t [s]'); xlim(T_WIN);
title(tl, sprintf('Trecho de validação %d–%d s: força da asa (coeficientes medidos) contra o empuxo dos rotores', T_WIN), 'FontWeight','bold');
fn2 = fullfile(paths.images,'force_terms_validation.png');
exportgraphics(f2, fn2, 'BackgroundColor','white','Resolution',140);
dd = fullfile(getenv('HOME'),'Desktop','DH_modelo_oficial');
copyfile(fn1, fullfile(dd,'force_terms_envelope.png'));  copyfile(fn2, fullfile(dd,'force_terms_validation.png'));
fprintf('\n  Figuras: %s\n           %s  (cópias na Mesa)\n', fn1, fn2);

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
