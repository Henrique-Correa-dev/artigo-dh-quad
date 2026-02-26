% motor_mapping_check.m
% Diagnóstico: identifica qual canal PWM (C1-C4) corresponde a qual
% posição física do motor, usando correlação cruzada entre PWM e p_dot, q_dot, r_dot.
%
% Lógica:
%   - Se PWM_i correlaciona positivamente com p_dot → motor do lado ESQUERDO
%   - Se PWM_i correlaciona negativamente com p_dot → motor do lado DIREITO
%   - Se PWM_i correlaciona positivamente com q_dot → motor FRONTAL
%   - Se PWM_i correlaciona negativamente com q_dot → motor TRASEIRO
%   - Correlação com r_dot indica sentido CW/CCW

%% 1. Carregar e interpolar dados (mesmo que new_identification)
load("log_data.mat")

ATT.TimeS  = double(ATT.TimeUS) / 1e6;
IMU.TimeS  = double(IMU.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;

idx = IMU.I == 0;
gyrX_raw = IMU.GyrX(idx); gyrY_raw = IMU.GyrY(idx); gyrZ_raw = IMU.GyrZ(idx);
time_IMU = IMU.TimeS(idx);
time_RCOU = RCOU.TimeS;

pwm1_raw = double(RCOU.C1); pwm2_raw = double(RCOU.C2);
pwm3_raw = double(RCOU.C3); pwm4_raw = double(RCOU.C4);

t_start = max([min(time_IMU), min(time_RCOU)]);
t_end   = min([max(time_IMU), max(time_RCOU)]);
dt = 0.1;
t_common = t_start:dt:t_end;

p = interp1(time_IMU, gyrX_raw, t_common, 'linear')';
q = interp1(time_IMU, gyrY_raw, t_common, 'linear')';
r = interp1(time_IMU, gyrZ_raw, t_common, 'linear')';

pwm1 = interp1(time_RCOU, pwm1_raw, t_common, 'linear')';
pwm2 = interp1(time_RCOU, pwm2_raw, t_common, 'linear')';
pwm3 = interp1(time_RCOU, pwm3_raw, t_common, 'linear')';
pwm4 = interp1(time_RCOU, pwm4_raw, t_common, 'linear')';

%% 2. Usar janela de treino com mais excitação (147-187s)
t_range = [147, 187];
idx_r = (t_common >= t_range(1)) & (t_common <= t_range(2));
t = t_common(idx_r)';
p_r = p(idx_r); q_r = q(idx_r); r_r = r(idx_r);
pwm_r = [pwm1(idx_r), pwm2(idx_r), pwm3(idx_r), pwm4(idx_r)];

%% 3. Derivadas numéricas suavizadas
smooth_win = 5;
p_dot = gradient(movmean(p_r, smooth_win), dt);
q_dot = gradient(movmean(q_r, smooth_win), dt);
r_dot = gradient(movmean(r_r, smooth_win), dt);

% Detrend PWM (remover média — foco na variação)
pwm_det = pwm_r - mean(pwm_r);

%% 4. Correlação cruzada (normalizada)
fprintf('\n==========================================================\n');
fprintf('  DIAGNÓSTICO DE MAPEAMENTO DE MOTORES\n');
fprintf('  Janela: %d-%ds\n', t_range(1), t_range(2));
fprintf('==========================================================\n\n');

corr_matrix = zeros(4, 3);
labels_axis = {'p_dot (Roll)', 'q_dot (Pitch)', 'r_dot (Yaw)'};
labels_motor = {'C1', 'C2', 'C3', 'C4'};

for i = 1:4
    corr_matrix(i, 1) = corr(pwm_det(:,i), p_dot);
    corr_matrix(i, 2) = corr(pwm_det(:,i), q_dot);
    corr_matrix(i, 3) = corr(pwm_det(:,i), r_dot);
end

fprintf('  Correlação PWM vs derivadas angulares:\n\n');
fprintf('  %6s  |  %12s  %12s  %12s\n', '', labels_axis{1}, labels_axis{2}, labels_axis{3});
fprintf('  %s\n', repmat('-', 1, 58));
for i = 1:4
    fprintf('  %6s  |  %+12.4f  %+12.4f  %+12.4f\n', ...
        labels_motor{i}, corr_matrix(i,1), corr_matrix(i,2), corr_matrix(i,3));
end

%% 5. Interpretar resultados
fprintf('\n  --- Interpretação ---\n');
fprintf('  (corr > 0 com p_dot → motor ESQUERDO, < 0 → DIREITO)\n');
fprintf('  (corr > 0 com q_dot → motor FRONTAL,  < 0 → TRASEIRO)\n');
fprintf('  (corr com r_dot → indica sentido CW/CCW)\n\n');

pos_labels = {'??', '??', '??', '??'};
for i = 1:4
    if corr_matrix(i,1) > 0 && corr_matrix(i,2) > 0
        pos_labels{i} = 'FRENTE-ESQ';
    elseif corr_matrix(i,1) > 0 && corr_matrix(i,2) < 0
        pos_labels{i} = 'TRÁS-ESQ';
    elseif corr_matrix(i,1) < 0 && corr_matrix(i,2) > 0
        pos_labels{i} = 'FRENTE-DIR';
    elseif corr_matrix(i,1) < 0 && corr_matrix(i,2) < 0
        pos_labels{i} = 'TRÁS-DIR';
    end

    if corr_matrix(i,3) > 0
        yaw_dir = 'CW (visto de cima)';
    else
        yaw_dir = 'CCW (visto de cima)';
    end

    fprintf('  %s → %s  |  Rotação: %s\n', labels_motor{i}, pos_labels{i}, yaw_dir);
end

%% 6. Comparar com o modelo atual
fprintf('\n  --- Modelo atual assume ---\n');
fprintf('  Mx = -(0.232*T1 - 0.232*T2 - 0.232*T3 + 0.232*T4)\n');
fprintf('     → T1,T4 roll negativo (DIREITO)  |  T2,T3 roll positivo (ESQUERDO)\n');
fprintf('  My = 0.311*T1 - 0.343*T2 + 0.311*T3 - 0.343*T4\n');
fprintf('     → T1,T3 pitch positivo (FRONTAL)  |  T2,T4 pitch negativo (TRASEIRO)\n');
fprintf('  Mz = Q1 + Q2 - Q3 - Q4\n');
fprintf('     → M1,M2 yaw positivo (CW)  |  M3,M4 yaw negativo (CCW)\n\n');

fprintf('  Modelo:  C1=FRENTE-DIR | C2=TRÁS-ESQ | C3=FRENTE-ESQ | C4=TRÁS-DIR\n');
fprintf('  Dados:   C1=%-11s | C2=%-11s | C3=%-11s | C4=%-11s\n', ...
    pos_labels{1}, pos_labels{2}, pos_labels{3}, pos_labels{4});

% Verificar se há mismatch
model_pos = {'FRENTE-DIR', 'TRÁS-ESQ', 'FRENTE-ESQ', 'TRÁS-DIR'};
n_match = 0;
for i = 1:4
    if strcmp(pos_labels{i}, model_pos{i})
        n_match = n_match + 1;
    end
end
fprintf('\n  >>> Match: %d/4\n', n_match);
if n_match < 4
    fprintf('  >>> ATENÇÃO: Mapeamento de motores inconsistente!\n');
    fprintf('  >>> Corrija Mx, My, Mz ou reordene os canais PWM.\n');
end

%% 7. Plots
fig = figure('Position', [50 50 1100 900], 'Visible', 'off');

subplot(4,1,1);
plot(t, pwm_r(:,1), t, pwm_r(:,2), t, pwm_r(:,3), t, pwm_r(:,4), 'LineWidth', 1.2);
legend('C1','C2','C3','C4'); ylabel('PWM'); grid on;
title('Sinais PWM individuais');

subplot(4,1,2);
yyaxis left;
plot(t, p_dot, 'b-', 'LineWidth', 1.3); ylabel('p\_dot (rad/s²)');
yyaxis right;
plot(t, pwm_det(:,1), '--', t, pwm_det(:,2), '--', t, pwm_det(:,3), '--', t, pwm_det(:,4), '--');
ylabel('\DeltaPWM'); legend('p\_dot','C1','C2','C3','C4'); grid on;
title('Roll: p\_dot vs \DeltaPWM');

subplot(4,1,3);
yyaxis left;
plot(t, q_dot, 'b-', 'LineWidth', 1.3); ylabel('q\_dot (rad/s²)');
yyaxis right;
plot(t, pwm_det(:,1), '--', t, pwm_det(:,2), '--', t, pwm_det(:,3), '--', t, pwm_det(:,4), '--');
ylabel('\DeltaPWM'); legend('q\_dot','C1','C2','C3','C4'); grid on;
title('Pitch: q\_dot vs \DeltaPWM');

subplot(4,1,4);
yyaxis left;
plot(t, r_dot, 'b-', 'LineWidth', 1.3); ylabel('r\_dot (rad/s²)');
yyaxis right;
plot(t, pwm_det(:,1), '--', t, pwm_det(:,2), '--', t, pwm_det(:,3), '--', t, pwm_det(:,4), '--');
ylabel('\DeltaPWM'); xlabel('Tempo (s)');
legend('r\_dot','C1','C2','C3','C4'); grid on;
title('Yaw: r\_dot vs \DeltaPWM');

sgtitle('Diagnóstico de Mapeamento de Motores');
saveas(fig, fullfile('C:/Users/Henrique/ARTIGO/identification/images', 'motor_mapping.png'));

% Heatmap da correlação
fig2 = figure('Position', [200 200 500 350], 'Visible', 'off');
imagesc(corr_matrix');
colorbar; colormap(redblue_cmap());
caxis([-1 1]);
set(gca, 'XTick', 1:4, 'XTickLabel', labels_motor);
set(gca, 'YTick', 1:3, 'YTickLabel', {'Roll (p\_dot)', 'Pitch (q\_dot)', 'Yaw (r\_dot)'});
title('Correlação: PWM vs Derivadas Angulares');
for i = 1:4
    for j = 1:3
        text(i, j, sprintf('%.3f', corr_matrix(i,j)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    end
end
saveas(fig2, fullfile('C:/Users/Henrique/ARTIGO/identification/images', 'motor_corr_heatmap.png'));

fprintf('\nScript finalizado.\n');

function cmap = redblue_cmap()
    n = 256;
    r = [linspace(0,1,n/2), ones(1,n/2)];
    g = [linspace(0,1,n/2), linspace(1,0,n/2)];
    b = [ones(1,n/2), linspace(1,0,n/2)];
    cmap = [r', g', b'];
end
