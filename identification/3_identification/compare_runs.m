% compare_runs.m — tabela comparativa das rodadas de identificação
% =========================================================================
% Lê outputs/runs/*/summary.mat e imprime, lado a lado, os parâmetros que
% interessam e o R² da janela de validação de cada rodada. Também grava
% outputs/runs/compare_runs.csv.
%
% Uso:  >> compare_runs
% =========================================================================
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();
d = dir(fullfile(paths.outputs,'runs','*'));
d = d([d.isdir] & ~startsWith({d.name},'.'));

rows = {};
fprintf('\n%-17s %-11s | %6s %6s %6s | %6s %6s %6s | %6s %6s %6s %6s\n', ...
    'rodada','validação','R² p','R² q','R² r','R² ax','R² ay','R² az','Jx','Jy','c_p','c_q');
fprintf('%s\n', repmat('-',1,110));
for i = 1:numel(d)
    fn = fullfile(d(i).folder, d(i).name, 'summary.mat');
    if ~exist(fn,'file'), continue; end
    S = load(fn).summary;  P = S.P_final;
    fprintf('%-17s %4.0f-%-6.0f | %6.3f %6.3f %6.3f | %6.2f %6.2f %6.2f | %6.4f %6.4f %6.2f %6.2f\n', ...
        S.tag, S.t_val(1), S.t_val(2), S.R2_val, S.R2_acc, P(1), P(2), P(13), P(14));
    rows(end+1,:) = {S.tag, S.t_val(1), S.t_val(2), S.R2_val(1), S.R2_val(2), S.R2_val(3), ...
        S.R2_acc(1), S.R2_acc(2), S.R2_acc(3), P(1), P(2), P(13), P(14), P(15), ...
        mean(P(5:8)), mean(P(9:12))}; %#ok<SAGROW>
end

fprintf('\nCoeficientes aerodinâmicos identificados (só as rodadas com eles livres):\n');
fprintf('%-17s %8s %8s %8s %8s %8s %8s %8s\n', 'rodada','C_d','Cl_p','Cl_b','Cm_q','Cm_a','Cn_r','Cn_b');
for i = 1:numel(d)
    fn = fullfile(d(i).folder, d(i).name, 'summary.mat');
    if ~exist(fn,'file'), continue; end
    S = load(fn).summary;
    if numel(S.P_final) < 22 || ~any(S.P_final(16:22) ~= 0), continue; end
    fprintf('%-17s %8.4f %8.3f %8.3f %8.2f %8.3f %8.3f %8.3f\n', S.tag, S.P_final(16:22));
end

if ~isempty(rows)
    T = cell2table(rows, 'VariableNames', {'rodada','t0','t1','R2_p','R2_q','R2_r', ...
        'R2_ax','R2_ay','R2_az','Jx','Jy','cp','cq','cr','kT_medio','kQ_medio'});
    fn = fullfile(paths.outputs,'runs','compare_runs.csv');
    writetable(T, fn);
    fprintf('\n  CSV: %s\n', fn);
end
