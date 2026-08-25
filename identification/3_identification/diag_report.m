function T = diag_report()
%DIAG_REPORT  Tabela-resumo das rodadas em outputs/diag/*/summary.mat
%   Colunas: tag, tau_m, cp/cq/cr (P_final), CR% de cp/cq/cr, Jx/Jy, kT médio,
%   R² de validação (p,q,r) do P_final e do híbrido, brancura, tempo.
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();
d = dir(fullfile(paths.outputs,'diag','*','summary.mat'));
if isempty(d), fprintf('diag_report: nenhuma rodada em %s\n', fullfile(paths.outputs,'diag')); T = []; return; end
rows = {};
for i = 1:numel(d)
    S = load(fullfile(d(i).folder, d(i).name)).summary;
    P = S.P_final;  cr = S.CR_pct;
    lock = isfield(S.DIAG,'lock_damp');
    dly = 1; if isfield(S.DIAG,'delay_pwm'), dly = S.DIAG.delay_pwm; end
    rows(end+1,:) = {S.tag, S.tau_lag, dly, lock, P(13), P(14), P(15), cr(13), cr(14), cr(15), ...
        P(1), P(2), mean(P(5:8)), mean(P(9:12)), ...
        S.R2_final(1), S.R2_final(2), S.R2_final(3), S.R2_hyb(1), S.R2_hyb(2), S.R2_hyb(3), ...
        S.white_pct(1), S.white_pct(2), S.white_pct(3), S.R2_P0(1), S.R2_P0(2), S.R2_P0(3), S.elapsed_s/60}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'tag','tau_m','delay','lock','cp','cq','cr','CRcp','CRcq','CRcr', ...
    'Jx','Jy','kT_med','kQ_med','R2p','R2q','R2r','R2p_hyb','R2q_hyb','R2r_hyb', ...
    'white_p','white_q','white_r','R2p_P0','R2q_P0','R2r_P0','min'});
fprintf('\n=================== DIAG REPORT (%d rodadas) ===================\n', height(T));
fprintf('%-11s %5s %3s %4s | %6s %6s %6s | %5s %5s %5s | %6s %6s | %5s %5s | %6s %6s %6s | %6s %6s %6s | %4s %4s %4s | %5s\n', ...
    'tag','tau','dly','lock','cp','cq','cr','CR%cp','CR%cq','CR%cr','Jx','Jy','kT','kQ', ...
    'R2p','R2q','R2r','R2p_h','R2q_h','R2r_h','w_p','w_q','w_r','min');
for i = 1:height(T)
    fprintf('%-11s %5.2f %3d %4d | %6.2f %6.2f %6.3f | %5.1f %5.1f %5.1f | %6.4f %6.4f | %5.2f %5.2f | %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f | %4.0f %4.0f %4.0f | %5.1f\n', ...
        T.tag{i}, T.tau_m(i), T.delay(i), T.lock(i), T.cp(i), T.cq(i), T.cr(i), T.CRcp(i), T.CRcq(i), T.CRcr(i), ...
        T.Jx(i), T.Jy(i), T.kT_med(i), T.kQ_med(i), T.R2p(i), T.R2q(i), T.R2r(i), ...
        T.R2p_hyb(i), T.R2q_hyb(i), T.R2r_hyb(i), T.white_p(i), T.white_q(i), T.white_r(i), T.min(i));
end
fprintf('  R² P0 (p,q,r) por rodada: ');
for i = 1:height(T), fprintf('%s=[%.2f %.2f %.2f] ', T.tag{i}, T.R2p_P0(i), T.R2q_P0(i), T.R2r_P0(i)); end
fprintf('\n');
writetable(T, fullfile(paths.outputs,'diag','diag_report.csv'));
end
