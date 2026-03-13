% generate_bench_figure.m — Gera figura de bancada (Thrust e Torque) para o artigo
%
% Saída: images/bench_motor.png (duas curvas lado a lado)

pwm_exp   = [1000; 1200; 1400; 1600; 1800; 2000];
thrust_g  = [0; 143; 328; 532; 784; 843];
torque_Nm = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];

% Converter para Newtons
thrust_N = thrust_g * 9.80665 / 1000;

% Interpolação fina para curva suave
pwm_fine = linspace(1000, 2000, 200);

% Ajuste polinomial grau 3 (mesmo do modelo)
cT = polyfit(pwm_exp, thrust_N, 3);
cQ = polyfit(pwm_exp, torque_Nm, 3);

thrust_fine = max(0, polyval(cT, pwm_fine));
torque_fine = max(0, polyval(cQ, pwm_fine));

fig = figure('Position', [100 100 800 320], 'Color', 'w');

% --- Subplot 1: Empuxo ---
subplot(1,2,1);
plot(pwm_fine, thrust_fine, 'b-', 'LineWidth', 1.5); hold on;
plot(pwm_exp, thrust_N, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
hold off;
xlabel('PWM (\mus)', 'FontSize', 10);
ylabel('Empuxo (N)', 'FontSize', 10);
title('Empuxo vs PWM', 'FontSize', 11);
legend('Polinômio 3° grau', 'Dados experimentais', 'Location', 'northwest', 'FontSize', 8);
grid on;
xlim([950 2050]);

% --- Subplot 2: Torque ---
subplot(1,2,2);
plot(pwm_fine, torque_fine, 'b-', 'LineWidth', 1.5); hold on;
plot(pwm_exp, torque_Nm, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
hold off;
xlabel('PWM (\mus)', 'FontSize', 10);
ylabel('Torque (N\cdotm)', 'FontSize', 10);
title('Torque vs PWM', 'FontSize', 11);
legend('Polinômio 3° grau', 'Dados experimentais', 'Location', 'northwest', 'FontSize', 8);
grid on;
xlim([950 2050]);

% Salvar
img_dir = fullfile(fileparts(mfilename('fullpath')), 'images');
saveas(fig, fullfile(img_dir, 'bench_motor.png'));
fprintf('Salvo: %s\n', fullfile(img_dir, 'bench_motor.png'));
