function B = estimate_bias(L, opts)
%ESTIMATE_BIAS  Estima bias do gyro e do acelerômetro a partir de um log de voo.
%
% Estratégia (sem hardware extra, robusto):
%
%   1. Detecta JANELAS CALMAS automaticamente:
%      (a) Pré-arme:    PWM_max < pwm_idle_thr  (drone no chão, motores em idle)
%      (b) Hover steady: gyro_norm < gyro_thr  AND  |accel_norm - g| < acc_thr
%                       durante pelo menos min_duration segundos
%
%   2. Em cada janela calma, calcula:
%      bias_gyro_xyz  = mediana(gyro)                (drone NÃO está rotacionando)
%      bias_accel_xyz = mediana(accel) - g_body      (drone em equilíbrio com g)
%
%      onde g_body = R(phi,theta)' * [0;0;g] projeta a gravidade no body frame
%      usando a atitude do EKF (que está bem mais limpa que integrar gyro).
%
%   3. Combina TODAS as janelas via mediana ponderada por número de amostras.
%
%   4. Reporta std entre janelas (indicador de qualidade).
%
% USO:
%   L = load_log_data('1_data/4 25-05-2026 ...mat');
%   B = estimate_bias(L);                  % opções default
%   B = estimate_bias(L, struct('verbose', true, 'plot', true));
%
% SAÍDA (struct B):
%   B.gyro          [3x1] vetor de bias do gyro [rad/s]   (subtrair: gyr_clean = gyr_raw - B.gyro)
%   B.accel         [3x1] vetor de bias do accel [m/s²]   (subtrair: acc_clean = acc_raw - B.accel)
%   B.gyro_std      [3x1] desvio padrão entre janelas
%   B.accel_std     [3x1] desvio padrão entre janelas
%   B.n_windows     [1x1] número de janelas calmas detectadas
%   B.n_samples     [1x1] total de amostras usadas
%   B.windows       [Nx2] tabela de [t_start, t_end] de cada janela
%   B.opts_used     struct com opções efetivas

% Opções default
default_opts = struct( ...
    'pwm_idle_thr',  1100,  ...   % PWM <= isto = motor idle (drone no chão)
    'gyro_thr',      0.10,  ...   % rad/s — limiar de "calmo" rotacional
    'acc_thr',       0.50,  ...   % m/s² — |accel_norm - g| máximo aceitável
    'min_duration',  0.3,   ...   % segundos — janela mínima
    'g',             9.81,  ...
    'verbose',       true,  ...
    'plot',          false);

if nargin < 2, opts = struct(); end
opts = merge_opts(default_opts, opts);

%% 1. Grade comum (IMU é o sinal mais rápido)
t = L.time_IMU;
N = numel(t);

% Interpolar PWM e atitude pra grade do IMU
pwm = [interp1(L.time_RCOU, L.pwm1_raw, t, 'previous', 'extrap'), ...
       interp1(L.time_RCOU, L.pwm2_raw, t, 'previous', 'extrap'), ...
       interp1(L.time_RCOU, L.pwm3_raw, t, 'previous', 'extrap'), ...
       interp1(L.time_RCOU, L.pwm4_raw, t, 'previous', 'extrap')];

phi   = interp1(L.time_ATT, deg2rad(L.roll_deg),  t, 'linear', 'extrap');
theta = interp1(L.time_ATT, deg2rad(L.pitch_deg), t, 'linear', 'extrap');

gyr = [L.gyrX_raw, L.gyrY_raw, L.gyrZ_raw];
acc = [L.accX_raw, L.accY_raw, L.accZ_raw];

gyr_norm = sqrt(sum(gyr.^2, 2));
acc_norm = sqrt(sum(acc.^2, 2));

%% 2. Critérios de "janela calma"
%
%   (a) PRÉ-ARME: PWM máximo dos 4 motores <= pwm_idle_thr
%   (b) HOVER STEADY: gyro_norm < thr  AND  |acc_norm - g| < acc_thr
%
is_idle  = max(pwm, [], 2) <= opts.pwm_idle_thr;
is_hover = (gyr_norm < opts.gyro_thr) & (abs(acc_norm - opts.g) < opts.acc_thr) & ~is_idle;
is_calm  = is_idle | is_hover;

%% 3. Encontrar segmentos contínuos com mín duração
[seg_start, seg_end] = find_segments(is_calm, t, opts.min_duration);
n_windows = numel(seg_start);

if n_windows == 0
    warning('estimate_bias:noWindow', ...
        'Nenhuma janela calma detectada com critérios atuais. Relaxe gyro_thr/acc_thr.');
    B = empty_result();
    return;
end

%% 4. Para cada janela, calcular bias parcial (mediana robusta)
gyro_per_win  = zeros(n_windows, 3);
accel_per_win = zeros(n_windows, 3);
n_samp_per_win = zeros(n_windows, 1);
typ_per_win = strings(n_windows, 1);

for w = 1:n_windows
    idx_w = (t >= seg_start(w)) & (t <= seg_end(w));
    n_samp_per_win(w) = sum(idx_w);

    % gyro bias = mediana do gyro (drone parado)
    gyro_per_win(w,:) = median(gyr(idx_w, :), 1);

    % accel bias = mediana(acc) - g_body_projetada(atitude_EKF)
    phi_w = median(phi(idx_w));
    theta_w = median(theta(idx_w));
    % Gravidade no body frame (NED com z pra baixo):
    %   g_body = R_b/i · g_inertial,  g_inertial = [0;0;+g]
    %   g_body_x = -g·sin(theta)
    %   g_body_y =  g·sin(phi)·cos(theta)
    %   g_body_z =  g·cos(phi)·cos(theta)
    %
    % O acelerômetro lê specific_force = (a_linear_body) - g_body
    % Em equilíbrio (a_linear = 0): a_imu = -g_body
    g_body = [-opts.g*sin(theta_w); opts.g*sin(phi_w)*cos(theta_w); opts.g*cos(phi_w)*cos(theta_w)];
    expected_acc = -g_body;
    accel_per_win(w,:) = median(acc(idx_w, :), 1) - expected_acc';

    % Tipo da janela
    if max(max(pwm(idx_w,:))) <= opts.pwm_idle_thr
        typ_per_win(w) = "idle";
    else
        typ_per_win(w) = "hover";
    end
end

%% 5. Bias agregado (mediana ponderada por n_samples)
% Para mediana ponderada simples: repetir cada janela n_samp vezes e tirar mediana
gyro_expanded  = expand_for_median(gyro_per_win,  n_samp_per_win);
accel_expanded = expand_for_median(accel_per_win, n_samp_per_win);

B = struct();
B.gyro       = median(gyro_expanded, 1)';
B.accel      = median(accel_expanded, 1)';
B.gyro_std   = std(gyro_per_win, 0, 1)';
B.accel_std  = std(accel_per_win, 0, 1)';
B.n_windows  = n_windows;
B.n_samples  = sum(n_samp_per_win);
B.windows    = [seg_start(:), seg_end(:), n_samp_per_win];   % col 3 = n_samples
B.window_types = typ_per_win;
B.opts_used  = opts;

%% 6. Reportar
if opts.verbose
    fprintf('\n=== estimate_bias: %d janelas calmas (%d amostras totais) ===\n', ...
        n_windows, B.n_samples);
    n_idle  = sum(typ_per_win == "idle");
    n_hover = sum(typ_per_win == "hover");
    fprintf('  Tipos: %d idle (pre-arme), %d hover steady\n', n_idle, n_hover);

    fprintf('\n  BIAS GYRO:\n');
    nm = {'X','Y','Z'};
    for i = 1:3
        fprintf('    bias_gyro_%s = %+9.5f rad/s  (%.4f°/s)   std entre janelas = %.5f\n', ...
            nm{i}, B.gyro(i), rad2deg(B.gyro(i)), B.gyro_std(i));
    end

    fprintf('\n  BIAS ACCEL (acc_raw - g_body_esperado):\n');
    for i = 1:3
        fprintf('    bias_accel_%s = %+8.4f m/s²   std entre janelas = %.4f\n', ...
            nm{i}, B.accel(i), B.accel_std(i));
    end

    fprintf('\n  COMO APLICAR:\n');
    fprintf('    gyrX_clean = gyrX_raw - (%+.5f)\n', B.gyro(1));
    fprintf('    gyrY_clean = gyrY_raw - (%+.5f)\n', B.gyro(2));
    fprintf('    gyrZ_clean = gyrZ_raw - (%+.5f)\n', B.gyro(3));
end

%% 7. Plot opcional
if opts.plot
    plot_bias_diagnostic(L, t, pwm, gyr, acc, is_idle, is_hover, seg_start, seg_end, B);
end
end


%% ========================================================================
%  HELPERS
%  ========================================================================

function [seg_start, seg_end] = find_segments(is_calm, t, min_duration)
% Encontra segmentos contínuos onde is_calm == true e duração >= min_duration
d = diff([0; is_calm(:); 0]);
starts = find(d == +1);
ends   = find(d == -1) - 1;
seg_start = []; seg_end = [];
for k = 1:numel(starts)
    if ends(k) > numel(t), ends(k) = numel(t); end
    if starts(k) > numel(t), continue; end
    dur = t(ends(k)) - t(starts(k));
    if dur >= min_duration
        seg_start(end+1,1) = t(starts(k));   %#ok<AGROW>
        seg_end(end+1,1)   = t(ends(k));     %#ok<AGROW>
    end
end
end

function out = expand_for_median(vals, weights)
% Expande matriz Nx3 repetindo cada linha 'weights(i)' vezes para mediana ponderada
out = zeros(sum(weights), size(vals, 2));
idx = 1;
for k = 1:size(vals, 1)
    w = weights(k);
    out(idx:idx+w-1, :) = repmat(vals(k,:), w, 1);
    idx = idx + w;
end
end

function opts = merge_opts(default_opts, user_opts)
opts = default_opts;
if isstruct(user_opts)
    fields_u = fieldnames(user_opts);
    for k = 1:numel(fields_u)
        opts.(fields_u{k}) = user_opts.(fields_u{k});
    end
end
end

function B = empty_result()
B = struct('gyro', [NaN; NaN; NaN], 'accel', [NaN; NaN; NaN], ...
           'gyro_std', [NaN; NaN; NaN], 'accel_std', [NaN; NaN; NaN], ...
           'n_windows', 0, 'n_samples', 0, 'windows', [], 'window_types', strings(0,1));
end

function plot_bias_diagnostic(L, t, pwm, gyr, acc, is_idle, is_hover, seg_s, seg_e, B)
fig = figure('Position', [100 100 1400 800], 'Color','w');

% PWM com sombreamento de janelas
subplot(3,1,1); hold on; grid on;
plot(t, pwm, 'LineWidth', 0.8);
for w = 1:numel(seg_s)
    x = [seg_s(w), seg_e(w), seg_e(w), seg_s(w)];
    y = [1000, 1000, 2000, 2000];
    if any(B.window_types(w) == "idle")
        col = [0.5 0.8 0.5 0.3];
    else
        col = [0.5 0.5 0.8 0.3];
    end
    fill(x, y, col(1:3), 'FaceAlpha', 0.3, 'EdgeColor','none');
end
ylabel('PWM [\mus]');
title(sprintf('Janelas calmas detectadas (%d total: verde=idle, azul=hover)', numel(seg_s)));

% Gyro
subplot(3,1,2); hold on; grid on;
plot(t, gyr(:,1), 'r', 'DisplayName', 'GyrX');
plot(t, gyr(:,2), 'g', 'DisplayName', 'GyrY');
plot(t, gyr(:,3), 'b', 'DisplayName', 'GyrZ');
yline(B.gyro(1), 'r--', sprintf('bias X = %+.4f', B.gyro(1)));
yline(B.gyro(2), 'g--', sprintf('bias Y = %+.4f', B.gyro(2)));
yline(B.gyro(3), 'b--', sprintf('bias Z = %+.4f', B.gyro(3)));
ylabel('Gyro [rad/s]'); legend('Location','best');
title(sprintf('Bias gyro estimado: [%+.4f, %+.4f, %+.4f] rad/s', B.gyro(1), B.gyro(2), B.gyro(3)));

% Accel norm
subplot(3,1,3); hold on; grid on;
plot(t, sqrt(sum(acc.^2, 2)), 'k', 'DisplayName', '|accel|');
yline(B.opts_used.g, 'r--', 'g=9.81');
ylabel('|a| [m/s²]'); xlabel('t [s]'); legend('Location','best');
title('Norm do acelerômetro (deve estar próximo de g em hover steady)');

sgtitle('Diagnóstico estimate_bias');
out_path = '/Users/graest/ita-master/artigo/artigo-dh-quad/identification/outputs/images/bias_estimate.png';
saveas(fig, out_path);
fprintf('  Figura: %s\n', out_path);
end
