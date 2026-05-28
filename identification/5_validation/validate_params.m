%VALIDATE_PARAMS  Simula a planta com P_J fornecido em UMA OU MAIS janelas do
%                 log e gera comparações nos 3 modos (full, hybrid, semi).
%
% NÃO RODA OTIMIZAÇÃO — só simula. Útil pra:
%   - Testar P0 (chute inicial)
%   - Testar P_final identificado
%   - Testar mudanças manuais (ex: Jy do slide, Bz=0, etc)
%   - Comparar parâmetros diferentes lado a lado em VÁRIAS janelas
%
% Edite as configurações abaixo e rode:
%   >> validate_params
%
% Outputs: figuras em outputs/images/ com prefixo "VP_<tag>_w<n>_<modo>.png"
%          + tabela comparativa de R² no console no fim

clear; clc; close all;

addpath(fileparts(fileparts(mfilename('fullpath'))));   % raiz
setup_paths();

% ╔══════════════════════════════════════════════════════════════════╗
% ║  CONFIGURAÇÃO — EDITE AQUI                                       ║
% ╚══════════════════════════════════════════════════════════════════╝

LOG_FILE = 'logs_concat.mat';      % ou um log individual ex: '4 25-05-2026 ...'

% Janelas de validação (em segundos do log) — uma por LINHA da matriz
% Pode ser 1 janela (1x2) ou várias (Nx2).
T_WINDOWS = [ ...
    1,  30;     % janela curta no log 1 (referência)
];

% Modos a rodar (true = roda, false = pula)
RUN_FULL   = false;
RUN_HYBRID = true;
RUN_SEMI   = false;

% Linhas a aparecer nos plots (cada linha tem 3 colunas: X/Y/Z)
%   PLOT_PQR     -> p, q, r            (medido vs sim)
%   PLOT_MOMENTS -> Mx, My, Mz         (só sim — não há sensor de torque)
%   PLOT_ATT     -> phi, theta, psi    (medido vs sim)
%   PLOT_ACC     -> accX, accY, accZ   (medido vs sim)
PLOT_PQR     = true;
PLOT_MOMENTS = true;
PLOT_ATT     = false;
PLOT_ACC     = false;

% Colunas — isola eixos (1=X, 2=Y, 3=Z)
%   Ex: pra focar só em yaw (r, Mz, psi, accZ), deixa só PLOT_Z = true
PLOT_X = true;
PLOT_Y = true;
PLOT_Z = true;

% Fonte de P_J — escolher UMA das opções:
P_SOURCE = 'manual';        % 'P0' | 'P_final' | 'manual'

% Se P_SOURCE='manual', edita aqui (15 elementos — modelo rotacional puro):
%   [Jx Jy Jz Jxz | k_T1..4 | k_Q1..4 | Dp Dq Dr]
%   Removidos: Bp/Bq/Br (CG via Lx/Ly), Xu/Yv/Zw (drag não-identif.), Bz (vai pra sensor depois)
P_MANUAL = [0.05; 0.2; 0.150; 0.003; ...    % Jx Jy Jz Jxz
            0.4; 0.4; 0.4; 0.4; ...          % k_T1..4
            1.0; 1.0; 1.0; 1.0; ...          % k_Q1..4
            2.12; 2.06; 0.4];                % Dp Dq Dr


% Tag pros nomes dos arquivos de saída (curto, pra distinguir testes)
TAG = 'test1';   % gera "VP_test1_w1_full.png", "VP_test1_w2_full.png", etc.


%% ====== Pipeline (não precisa mexer abaixo) ======

% Validação básica
if size(T_WINDOWS, 2) ~= 2
    error('T_WINDOWS deve ser Nx2 (cada linha = [t_start, t_end]). Recebeu %dx%d.', size(T_WINDOWS));
end
n_wins = size(T_WINDOWS, 1);

% 1. Carrega P_J
proj = parameters();
switch lower(P_SOURCE)
    case 'p0'
        P_J = proj.P0_J;
        fprintf('validate_params: usando P0 (chute inicial).\n');
    case 'p_final'
        Pdat = load(fullfile(setup_paths().outputs, 'P_identified.mat'));
        P_J = Pdat.P_final;
        fprintf('validate_params: usando P_final de P_identified.mat.\n');
    case 'manual'
        P_J = P_MANUAL(:);
        fprintf('validate_params: usando P_MANUAL.\n');
    otherwise
        error('P_SOURCE inválido: %s', P_SOURCE);
end
if numel(P_J) ~= 15
    error('P_J deve ter 15 elementos. Recebeu %d.', numel(P_J));
end

% 2. Carrega log (uma vez, com fallback se concat não existir)
log_path = fullfile(setup_paths().data, LOG_FILE);
if ~exist(log_path, 'file')
    fprintf('⚠ Arquivo "%s" não encontrado.\n', LOG_FILE);
    fprintf('  Para usá-lo, rode primeiro: >> test_flight\n');
    fprintf('  Usando log individual como fallback.\n\n');
    LOG_FILE = '4 25-05-2026 09-31-48.log-132954.mat';
    log_path = fullfile(setup_paths().data, LOG_FILE);
end
L = load_log_data(log_path);
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
t_common_full = t_lo:0.1:t_hi;
t_common_full = t_common_full(:);

% 3. Constants
constants = struct('m', proj.m, 'g', proj.g, 'tau_motor', proj.tau_motor);
R2 = @(y_e, y_s) 1 - sum((y_e - y_s).^2) / max(sum((y_e - mean(y_e)).^2), 1e-12);

% Modos
modes = {};
if RUN_FULL,   modes{end+1} = 'full';   end
if RUN_HYBRID, modes{end+1} = 'hybrid'; end
if RUN_SEMI,   modes{end+1} = 'semi';   end

% Storage de R² por janela / modo / sinal (tabela final)
sig_names = {'p','q','r','phi','theta','psi','accX','accY','accZ'};
R2_table = nan(n_wins, numel(modes), numel(sig_names));

img_dir = setup_paths().images;
if ~exist(img_dir, 'dir'), mkdir(img_dir); end

%% ====== Loop sobre janelas ======
for w = 1:n_wins
    t_window = T_WINDOWS(w, :);
    fprintf('\n############### JANELA %d/%d: [%.1f, %.1f] s ###############\n', ...
        w, n_wins, t_window(1), t_window(2));

    if t_window(1) < t_lo || t_window(2) > t_hi
        warning('Janela [%g, %g] fora do log [%.1f, %.1f] — PULANDO', ...
            t_window(1), t_window(2), t_lo, t_hi);
        continue;
    end

    idx = (t_common_full >= t_window(1)) & (t_common_full <= t_window(2));
    time = t_common_full(idx);
    N = numel(time);

    pwm      = [interp1(L.time_RCOU, L.pwm1_raw, time, 'linear'), ...
                interp1(L.time_RCOU, L.pwm2_raw, time, 'linear'), ...
                interp1(L.time_RCOU, L.pwm3_raw, time, 'linear'), ...
                interp1(L.time_RCOU, L.pwm4_raw, time, 'linear')];
    pqr_meas = [interp1(L.time_IMU, L.gyrX_raw, time, 'linear'), ...
                interp1(L.time_IMU, L.gyrY_raw, time, 'linear'), ...
                interp1(L.time_IMU, L.gyrZ_raw, time, 'linear')];
    acc_meas = [interp1(L.time_IMU, L.accX_raw, time, 'linear'), ...
                interp1(L.time_IMU, L.accY_raw, time, 'linear'), ...
                interp1(L.time_IMU, L.accZ_raw, time, 'linear')];
    att_meas = [interp1(L.time_ATT, L.roll_deg,  time, 'linear'), ...
                interp1(L.time_ATT, L.pitch_deg, time, 'linear'), ...
                interp1(L.time_ATT, L.yaw_deg,   time, 'linear')];

    fprintf('  (%d amostras, dt=%.3f s)\n', N, time(2)-time(1));

    % Roda os modos
    results = struct();
    for k = 1:numel(modes)
        mode = modes{k};
        fprintf('\n  ===== modo: %s =====\n', mode);
        tic;
        results.(mode) = sim_window(mode, P_J, time, pwm, pqr_meas, att_meas, constants);
        fprintf('    Tempo sim: %.2f s\n', toc);

        r = results.(mode);
        all_meas = {pqr_meas(:,1), pqr_meas(:,2), pqr_meas(:,3), ...
                    att_meas(:,1), att_meas(:,2), att_meas(:,3), ...
                    acc_meas(:,1), acc_meas(:,2), acc_meas(:,3)};
        all_sim  = {r.p, r.q, r.r, r.phi, r.theta, r.psi, r.accX, r.accY, r.accZ};
        for s = 1:9
            R2_table(w, k, s) = R2(all_meas{s}, all_sim{s});
        end
        % Formato idêntico ao identify_plant.m → print_R2
        fprintf('  [%s] Validação (%g-%gs):\n', mode, t_window(1), t_window(2));
        fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
            R2_table(w,k,1), R2_table(w,k,2), R2_table(w,k,3));
        fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
            R2_table(w,k,7), R2_table(w,k,8), R2_table(w,k,9));
    end

    % Quais linhas plotar (cada uma vira uma linha de 3 colunas)
    plot_rows = {};
    if PLOT_PQR,     plot_rows{end+1} = 'pqr';     end
    if PLOT_MOMENTS, plot_rows{end+1} = 'moments'; end
    if PLOT_ATT,     plot_rows{end+1} = 'att';     end
    if PLOT_ACC,     plot_rows{end+1} = 'acc';     end

    plot_cols = [];   % índices das colunas (1=X, 2=Y, 3=Z)
    if PLOT_X, plot_cols(end+1) = 1; end
    if PLOT_Y, plot_cols(end+1) = 2; end
    if PLOT_Z, plot_cols(end+1) = 3; end

    if isempty(plot_rows) || isempty(plot_cols)
        warning('Nenhuma linha ou coluna selecionada (PLOT_*). Pulando plots.');
    else
        % Plots por modo
        for k = 1:numel(modes)
            mode = modes{k};
            r = results.(mode);
            plot_compare_rows(time, mode, TAG, w, t_window, plot_rows, plot_cols, ...
                pqr_meas, att_meas, acc_meas, r, img_dir);
        end

        % Overlay dos modos
        if numel(modes) > 1
            plot_modes_overlay_rows(time, modes, results, plot_rows, plot_cols, ...
                pqr_meas, att_meas, acc_meas, TAG, w, t_window, img_dir);
        end
    end
end

%% ====== Tabela comparativa final ======
fprintf('\n\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('  TABELA COMPARATIVA DE R² — todas as janelas / modos / sinais\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
for s = 1:numel(sig_names)
    fprintf('\n  %s:\n', sig_names{s});
    fprintf('    %-12s', 'janela');
    for k = 1:numel(modes), fprintf(' %10s', modes{k}); end
    fprintf('\n');
    for w = 1:n_wins
        fprintf('    [%4.0f,%4.0f] ', T_WINDOWS(w,1), T_WINDOWS(w,2));
        for k = 1:numel(modes)
            fprintf(' %+10.3f', R2_table(w, k, s));
        end
        fprintf('\n');
    end
end
fprintf('\n═══════════════════════════════════════════════════════════════════\n');
fprintf('Figuras salvas em: %s\n', img_dir);
fprintf('  VP_%s_w<n>_<mode>.png   (uma por janela × modo)\n', TAG);
if numel(modes) > 1
    fprintf('  VP_%s_w<n>_overlay.png (overlay dos modos por janela)\n', TAG);
end


%% ====== Plot helpers ======

% Define qual conteúdo vai em cada tipo de linha
function [labels, units, meas_cell, sim_cell, has_meas, ylims] = row_signals(row_type, ...
    pqr_m, att_m, acc_m, r)
% has_meas = true se há sinal medido pra comparar (pqr/att/acc).
% Pra 'moments', has_meas = false (sem sensor de torque).
% ylims = {[lo, hi], ...} pra cada coluna; [] = autoscale.
    switch row_type
        case 'pqr'
            labels = {'p','q','r'};
            units  = {'rad/s','rad/s','rad/s'};
            meas_cell = {pqr_m(:,1), pqr_m(:,2), pqr_m(:,3)};
            sim_cell  = {r.p, r.q, r.r};
            has_meas = true;
            ylims = {[-4 4], [-4 4], [-1.2 1.2]};
        case 'moments'
            labels = {'Mx','My','Mz'};
            units  = {'N·m','N·m','N·m'};
            meas_cell = {[], [], []};
            sim_cell  = {r.Mx, r.My, r.Mz};
            has_meas = false;
            ylims = {[-3 3], [-3 3], [-0.3 0.3]};
        case 'att'
            labels = {'phi','theta','psi'};
            units  = {'deg','deg','deg'};
            meas_cell = {att_m(:,1), att_m(:,2), att_m(:,3)};
            sim_cell  = {r.phi, r.theta, r.psi};
            has_meas = true;
            ylims = {[], [], []};
        case 'acc'
            labels = {'accX','accY','accZ'};
            units  = {'m/s²','m/s²','m/s²'};
            meas_cell = {acc_m(:,1), acc_m(:,2), acc_m(:,3)};
            sim_cell  = {r.accX, r.accY, r.accZ};
            has_meas = true;
            ylims = {[], [], []};
        otherwise
            error('row_type desconhecido: %s', row_type);
    end
end


function plot_compare_rows(time, mode, tag, w_idx, t_win, plot_rows, plot_cols, ...
    pqr_m, att_m, acc_m, r, img_dir)
% Plot Nx(M<=3) dinâmico (linhas = plot_rows, colunas = plot_cols).
    n_rows = numel(plot_rows);
    n_cols = numel(plot_cols);
    width = 500*n_cols;
    fig = figure('Position', [80 50 width 300*n_rows], 'Color', 'w', 'Visible', 'off');

    for ri = 1:n_rows
        [labels, units, meas, sim, has_meas, ylims] = row_signals(plot_rows{ri}, ...
            pqr_m, att_m, acc_m, r);
        for ci = 1:n_cols
            c = plot_cols(ci);   % índice real da coluna (1=X, 2=Y, 3=Z)
            subplot(n_rows, n_cols, (ri-1)*n_cols + ci); hold on; grid on;
            if has_meas
                plot(time, meas{c}, 'b-',  'LineWidth', 1.0, 'DisplayName', 'medido');
            end
            plot(time, sim{c},  'r--', 'LineWidth', 1.3, 'DisplayName', 'sim');
            title(labels{c}, 'FontWeight', 'bold');
            xlabel('t [s]'); ylabel(units{c});
            if ~isempty(ylims{c}), ylim(ylims{c}); end
            if ri == 1 && ci == 1, legend('Location','best'); end
        end
    end
    sgtitle(sprintf('validate\\_params — modo %s — w%d [%g, %g]s — tag %s', ...
        mode, w_idx, t_win(1), t_win(2), tag));
    saveas(fig, fullfile(img_dir, sprintf('VP_%s_w%d_%s.png', tag, w_idx, mode)));
    close(fig);
end


function plot_modes_overlay_rows(time, modes, results, plot_rows, plot_cols, ...
    pqr_m, att_m, acc_m, tag, w_idx, t_win, img_dir)
% Overlay dos modos no layout dinâmico.
    n_rows = numel(plot_rows);
    n_cols = numel(plot_cols);
    width = 500*n_cols;
    fig = figure('Position', [80 50 width 300*n_rows], 'Color', 'w', 'Visible', 'off');

    color_map = containers.Map({'full','hybrid','semi'}, {'r--','m-.','b:'});

    % Primeiro modo só pra obter labels/units (são os mesmos pra qualquer modo)
    ref_r = results.(modes{1});

    for ri = 1:n_rows
        [labels, units, meas, ~, has_meas, ylims] = row_signals(plot_rows{ri}, ...
            pqr_m, att_m, acc_m, ref_r);
        for ci = 1:n_cols
            c = plot_cols(ci);
            subplot(n_rows, n_cols, (ri-1)*n_cols + ci); hold on; grid on;
            if has_meas
                plot(time, meas{c}, 'Color', [0 0.7 0.3], 'LineWidth', 1.2, ...
                    'DisplayName', 'medido');
            end
            for m = 1:numel(modes)
                rm = results.(modes{m});
                [~, ~, ~, sim_m, ~, ~] = row_signals(plot_rows{ri}, ...
                    pqr_m, att_m, acc_m, rm);
                plot(time, sim_m{c}, color_map(modes{m}), 'LineWidth', 1.2, ...
                    'DisplayName', modes{m});
            end
            title(labels{c}, 'FontWeight','bold');
            xlabel('t [s]'); ylabel(units{c});
            if ~isempty(ylims{c}), ylim(ylims{c}); end
            if ri == 1 && ci == 1, legend('Location','best'); end
        end
    end
    sgtitle(sprintf('validate\\_params — overlay — w%d [%g, %g]s — tag %s', ...
        w_idx, t_win(1), t_win(2), tag));
    saveas(fig, fullfile(img_dir, sprintf('VP_%s_w%d_overlay.png', tag, w_idx)));
    close(fig);
end
