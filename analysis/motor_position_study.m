%MOTOR_POSITION_STUDY  Verifica se o mapeamento motor→posição no modelo
%                       é consistente com o que aparece nos PWMs do voo real.
%
% HIPÓTESE DO USUÁRIO:
%   Olhando o log, PWMs aparecem em DOIS PATAMARES:
%       - Par A (C1, C2): patamar ALTO  (~1760 us)
%       - Par B (C3, C4): patamar BAIXO (~1690 us)
%   Diferença ~70 us, sistemática.
%
% INTERPRETAÇÃO FÍSICA:
%   Geometria do drone: Ly_f = 0.311 m  <  Ly_r = 0.343 m
%   Em hover steady (My = 0):
%       Ly_f · T_front = Ly_r · T_rear   →   T_front/T_rear = Ly_r/Ly_f = 1.10
%   Logo: FRONT precisa de ~10% MAIS empuxo (~50-70 us mais PWM).
%
% LOGO O PAR DE "PATAMAR ALTO" SÃO OS MOTORES FRONT (não os diagonais
% como a convenção ArduPilot QuadX sugere).
%
% Este script testa 2 hipóteses de mapeamento e diz qual minimiza My residual:
%
%   H1 (ArduPilot QuadX padrão, USADO no código atual):
%        M1=FR, M2=RL, M3=FL, M4=RR  →  front = {M1, M3}, rear = {M2, M4}
%
%   H2 (Sugerida pelo usuário, pares frente/trás):
%        M1, M2 = front;  M3, M4 = rear
%
%   Hipótese vencedora = a que dá |My_hover| ≈ 0 (drone equilibrado).

clear; clc; close all;

addpath(fileparts(fileparts(mfilename('fullpath'))));   % raiz do projeto
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'identification'));
setup_paths();

%% 1. Carregar log novo (motores novos, primeira identificação)
LOG_FILE = '4 25-05-2026 09-31-48.log-132954.mat';
L = load_log_data(fullfile(setup_paths().data, LOG_FILE));

%% 2. Detectar janelas de "quase hover" (gyro pequeno + |acc| ≈ g)
[func_T_ref, ~] = motor_models();
p = parameters();
m  = p.m;
g  = p.g;
Ly_f = p.arms.Ly_f;
Ly_r = p.arms.Ly_r;
Lx_base = p.arms.Lx_base;

% Interpolar tudo na grade do IMU
t   = L.time_IMU;
gyr = [L.gyrX_raw, L.gyrY_raw, L.gyrZ_raw];
acc = [L.accX_raw, L.accY_raw, L.accZ_raw];
pwm = [interp1(L.time_RCOU, L.pwm1_raw, t, 'previous','extrap'), ...
       interp1(L.time_RCOU, L.pwm2_raw, t, 'previous','extrap'), ...
       interp1(L.time_RCOU, L.pwm3_raw, t, 'previous','extrap'), ...
       interp1(L.time_RCOU, L.pwm4_raw, t, 'previous','extrap')];

phi   = interp1(L.time_ATT, deg2rad(L.roll_deg),  t, 'linear','extrap');
theta = interp1(L.time_ATT, deg2rad(L.pitch_deg), t, 'linear','extrap');

gyr_norm = sqrt(sum(gyr.^2, 2));
acc_norm = sqrt(sum(acc.^2, 2));

is_hover = (gyr_norm < 0.10) & (abs(acc_norm - g) < 0.5) & ...
           (max(pwm,[],2) > 1100) & (abs(phi) < deg2rad(5)) & (abs(theta) < deg2rad(5));

n_hover = sum(is_hover);
if n_hover < 20
    error('Poucas amostras de hover steady (%d). Relaxe critérios.', n_hover);
end
fprintf('Janelas de hover steady: %d amostras (%.1f s total)\n', n_hover, n_hover*(t(2)-t(1)));

% Janelas pra plotar (segmentos contínuos)
d = diff([0; is_hover; 0]);
seg_start_idx = find(d == +1);
seg_end_idx   = find(d == -1) - 1;
fprintf('  Segmentos contínuos: %d\n', numel(seg_start_idx));

%% 3. PWM e empuxo médios em hover (sob cada motor)
pwm_hover = pwm(is_hover, :);
fprintf('\n=== PWMs medianos em HOVER STEADY ===\n');
fprintf('  C1 = %.0f us\n', median(pwm_hover(:,1)));
fprintf('  C2 = %.0f us\n', median(pwm_hover(:,2)));
fprintf('  C3 = %.0f us\n', median(pwm_hover(:,3)));
fprintf('  C4 = %.0f us\n', median(pwm_hover(:,4)));

% Empuxo via bench data (k_T = 1.0)
T = zeros(size(pwm_hover));
for i = 1:4
    T(:,i) = func_T_ref(pwm_hover(:,i));
end
T_med = median(T, 1);
fprintf('\n=== Empuxos medianos (k_T=1, motor referência) ===\n');
fprintf('  T1 = %.2f N\n', T_med(1));
fprintf('  T2 = %.2f N\n', T_med(2));
fprintf('  T3 = %.2f N\n', T_med(3));
fprintf('  T4 = %.2f N\n', T_med(4));
fprintf('  Sum = %.2f N   (peso mg = %.2f N)\n', sum(T_med), m*g);

%% 4. Testar 2 hipóteses de mapeamento
%
%   H1 (atual no código): M1=FR(front), M3=FL(front), M2=RL(rear), M4=RR(rear)
%       My = Ly_f*T1 - Ly_r*T2 + Ly_f*T3 - Ly_r*T4
%
%   H2 (sugerida user): M1, M2 = front;  M3, M4 = rear
%       My = Ly_f*T1 + Ly_f*T2 - Ly_r*T3 - Ly_r*T4
%
%   (Variantes: trocar quais ficam à esquerda/direita afetam Mx, não My)

% Para o estudo de My usaremos médias (em hover My deveria ser ≈ 0)
My_H1 = zeros(n_hover,1);
My_H2 = zeros(n_hover,1);

for k = 1:n_hover
    Tk = T(k,:);
    My_H1(k) =  Ly_f*Tk(1) - Ly_r*Tk(2) + Ly_f*Tk(3) - Ly_r*Tk(4);   % atual (ArduPilot QuadX)
    My_H2(k) =  Ly_f*Tk(1) + Ly_f*Tk(2) - Ly_r*Tk(3) - Ly_r*Tk(4);   % pares frente/trás
end

fprintf('\n=== MOMENTO DE PITCH RESIDUAL (My) em hover ===\n');
fprintf('  Esperado: ~0 N·m (drone equilibrado)\n');
fprintf('  H1 (atual, ArduPilot QuadX, diagonais): mediana = %+.3f N·m  |  RMS = %.3f\n', ...
        median(My_H1), sqrt(mean(My_H1.^2)));
fprintf('  H2 (sugerida, pares frente/trás):       mediana = %+.3f N·m  |  RMS = %.3f\n', ...
        median(My_H2), sqrt(mean(My_H2.^2)));

if abs(median(My_H2)) < abs(median(My_H1))
    fprintf('\n  >>> H2 VENCE: pares frente/trás explicam melhor o equilíbrio.\n');
    fprintf('  >>> Mapeamento atual no código está provavelmente TROCADO.\n');
else
    fprintf('\n  >>> H1 vence: mapeamento atual está correto.\n');
end

%% 5. Análise de quanto m_y_wing seria necessário pra explicar o residual em H1
fprintf('\n=== ALTERNATIVA: M_y_wing pra "salvar" H1 ===\n');
fprintf('  Se H1 está certo, então em hover My_motor ≠ 0 indica momento externo.\n');
fprintf('  M_y_wing necessário = -median(My_H1) = %+.3f N·m\n', -median(My_H1));
fprintf('  Isso é compatível com massa da asa × distância CG? (depende do drone)\n');

%% 6. Plot diagnóstico
fig = figure('Position', [80 50 1500 950], 'Color','w');

% PWMs em hover (ponto por janela)
subplot(3,2,1); hold on; grid on;
boxplot(pwm_hover, {'C1','C2','C3','C4'});
ylabel('PWM [\mus]'); title('Distribuição PWM em hover steady');
yline(median(pwm_hover(:,1)), 'b:', 'C1', 'LineWidth', 0.5);
yline(median(pwm_hover(:,3)), 'r:', 'C3', 'LineWidth', 0.5);

% Pares de pwm vs tempo
subplot(3,2,2); hold on; grid on;
plot(t, pwm(:,1), 'r', 'DisplayName','C1');
plot(t, pwm(:,2), 'b', 'DisplayName','C2');
plot(t, pwm(:,3), 'g', 'DisplayName','C3');
plot(t, pwm(:,4), 'm', 'DisplayName','C4');
ylim([1500 2000]);
% Sombrear janelas de hover
for s = 1:numel(seg_start_idx)
    x_s = t(seg_start_idx(s));  x_e = t(seg_end_idx(s));
    patch([x_s x_e x_e x_s], [1500 1500 2000 2000], [0.6 0.6 0.6], ...
          'FaceAlpha', 0.15, 'EdgeColor','none', 'HandleVisibility','off');
end
xlabel('t [s]'); ylabel('PWM [\mus]');
title('PWMs (cinza = hover steady detectado)');
legend('Location','best');

% Histograma My sob cada hipótese
subplot(3,2,3); hold on; grid on;
histogram(My_H1, 30, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName','H1 (ArduPilot QuadX)');
histogram(My_H2, 30, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'DisplayName','H2 (pares frente/trás)');
xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('My residual [N·m]'); ylabel('contagem');
title('Distribuição de My em hover (esperado ~0)');
legend('Location','best');

% My vs tempo
subplot(3,2,4); hold on; grid on;
t_hov = t(is_hover);
plot(t_hov, My_H1, 'r.', 'MarkerSize', 8, 'DisplayName','H1');
plot(t_hov, My_H2, 'g.', 'MarkerSize', 8, 'DisplayName','H2');
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
xlabel('t [s]'); ylabel('My [N·m]');
title('My em cada amostra de hover');
legend('Location','best');

% Comparação numérica
subplot(3,2,5); axis off;
txt = {
    sprintf('\\bfDiagnóstico em %d amostras de hover steady\\rm', n_hover);
    sprintf('PWM medianos: C1=%.0f  C2=%.0f  C3=%.0f  C4=%.0f us', median(pwm_hover));
    sprintf('Diff (C1+C2)/2 vs (C3+C4)/2 = %.0f us', (median(pwm_hover(:,1))+median(pwm_hover(:,2)))/2 - (median(pwm_hover(:,3))+median(pwm_hover(:,4)))/2);
    '';
    '\bfHipótese 1 — ArduPilot QuadX (atual no código)\rm';
    sprintf('  My mediana = %+.3f N·m,  RMS = %.3f', median(My_H1), sqrt(mean(My_H1.^2)));
    '';
    '\bfHipótese 2 — Pares frente/trás (sugestão usuário)\rm';
    sprintf('  My mediana = %+.3f N·m,  RMS = %.3f', median(My_H2), sqrt(mean(My_H2.^2)));
    '';
    '\bfConclusão:\rm';
    sprintf('  Melhor (|My| menor): %s', ternary(abs(median(My_H1)) < abs(median(My_H2)), 'H1', 'H2'));
};
text(0.05, 0.5, txt, 'Interpreter','tex', 'FontSize', 11, 'VerticalAlignment','middle');

% Esquema dos 2 mappings
subplot(3,2,6); hold on; axis equal; axis off;
title('Esquema dos mapeamentos');
% Drone visto de cima — desenhar 2x: H1 e H2
draw_mapping(0, 0, 'H1: ArduPilot QuadX', {'M1\rightarrowFR','M3\rightarrowFL','M2\rightarrowRL','M4\rightarrowRR'});
draw_mapping(1.2, 0, 'H2: pares F/T', {'M1\rightarrowF','M2\rightarrowF','M3\rightarrowR','M4\rightarrowR'});

sgtitle(sprintf('Estudo de mapeamento motor\\rightarrowposição  |  log: %s', strrep(LOG_FILE,'_','\_')));

img_dir = fileparts(mfilename('fullpath'));
out_path = fullfile(img_dir, 'images', 'motor_position_study.png');
saveas(fig, out_path);
fprintf('\nFigura salva: %s\n', out_path);


%% ========== HELPERS ==========
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function draw_mapping(x0, y0, label, motor_labels)
% Drone simples visto de cima, marcando posições dos 4 motores
sz = 0.4;
% Body
rectangle('Position',[x0-sz/2, y0-sz/2, sz, sz], 'EdgeColor',[0.6 0.6 0.6]);
% Motors (FR, FL, RR, RL) — convenção física
positions = [+sz/2 +sz/2;  % FR
             -sz/2 +sz/2;  % FL
             +sz/2 -sz/2;  % RR
             -sz/2 -sz/2]; % RL
labels_phys = {'FR','FL','RR','RL'};
for i = 1:4
    plot(x0 + positions(i,1), y0 + positions(i,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor','w');
    text(x0 + positions(i,1)*1.3, y0 + positions(i,2)*1.3, labels_phys{i}, ...
        'HorizontalAlignment','center','FontSize',8);
end
text(x0, y0+sz, label, 'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
text(x0, y0-sz-0.05, strjoin(motor_labels, ', '), ...
    'HorizontalAlignment','center','FontSize', 8, 'Interpreter','tex');
end
