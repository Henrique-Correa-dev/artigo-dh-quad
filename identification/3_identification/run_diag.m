function run_diag(which)
%RUN_DIAG  Roda identify_plant_diag.m para um conjunto de casos de diagnóstico.
%   Objetivo: entender por que cp, cq, cr identificados (OEM) ficam bem acima
%   dos valores a priori (prior_damping.m) e o que o dado realmente exige.
%
%   Casos (outputs/diag/<tag>/summary.mat):
%     base        pipeline oficial tal como está (tau_m=0.05, P0 damp=0.5)
%     prior_p0    idem, com P0(cp,cq,cr) = a priori (influxo)  → P_final deve ser igual
%     tau010/015/020  lag do motor maior (P0 = a priori)      → cp cai quando tau sobe?
%     tau002/delay0/nolag  lag e atraso MENORES que o oficial  → cp cai quando o lag cai?
%     lock_prior  cp,cq,cr TRAVADOS no a priori (só rotor)     → quanto R² perde?
%     lock_phys   TRAVADOS em rotor + estrutura a 1,3 m/s      → idem
%
%   Uso:  run_diag                      % todos, na ordem acima
%         run_diag({'base','tau015'})   % subconjunto
%   Relatório: diag_report
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
pr = load(fullfile(paths.outputs,'prior_damping.mat')).prior;
c_prior = [pr.cp0, pr.cq0, pr.cr0];
V_env   = 1.3;                                          % p95 da velocidade GPS na campanha
c_phys  = c_prior + V_env*[pr.Lp_per_V/p.J.Jx, pr.Mq_per_V/p.J.Jy, pr.Nr_per_V/p.J.Jz];

cases = { ...
    struct('tag','base'), ...
    struct('tag','base_cr'), ...   % = base, mas com CR%% também em Θ0 e Θ_EEM (Tabela 5.3 completa)
    struct('tag','base_cr05', 'P0_damp', [0.5 0.5 0.5]), ...
    struct('tag','prior_cr',  'P0_damp', c_prior), ...
    struct('tag','freeJ_cr',  'P0_damp', c_prior, 'lb_J', [0.03243; 0.06330; 0.09464; 0.00000], ...
                                                 'ub_J', [0.05406; 0.10551; 0.15774; 0.01000]), ...   % Jz, Jxz LIVRES + CR%% nos 3 pontos
    ... % ---- planta UNIFICADA (motor_chain) a partir daqui ----
    struct('tag','unified'), ...                              % = identify_plant atual (referência nova)
    struct('tag','lock_kT1', 'lock_kT', 1.0), ...             % cT travado na bancada: az volta? p,q perdem?
    ... % TESTE: amortecimento aerodinâmico ∝ V (Eq. 13 de Salahudden et al.: Cl_p·pb/2V etc.) no lugar de cp,cq,cr
    struct('tag','damp_aero',   'damp_mode','aero'), ...     % só ∝ V: identifica Cl_p, Cm_q, Cn_r (c0 = 0)
    struct('tag','damp_hybrid', 'damp_mode','hybrid'), ...   % rotor (a priori, fixo) + ∝ V identificado
    ... % ARRASTO INDUZIDO DO ROTOR (Beard 14.4.2): varredura de C_d. O drag_probe com a
    ... % velocidade estimada (voos em ginásio, sem vento; GPS com interferência) deu
    ... % C_d ≈ 0,008 (x) e 0,019 (y), com o sinal correto. Aqui se vê o efeito no ajuste.
    struct('tag','cd010', 'drag_cd', 0.010), ...
    struct('tag','cd020', 'drag_cd', 0.020), ...
    struct('tag','cd050', 'drag_cd', 0.050), ...   % Θ0 do amortecimento = a priori (influxo), com CR%% nos 3 pontos   % idem, com o Θ0 EXATO da Tabela 5.3 (cp=cq=cr=0,5)
    struct('tag','prior_p0',   'P0_damp', c_prior), ...
    struct('tag','tau010',     'P0_damp', c_prior, 'tau_m', 0.10), ...
    struct('tag','tau015',     'P0_damp', c_prior, 'tau_m', 0.15), ...
    struct('tag','tau020',     'P0_damp', c_prior, 'tau_m', 0.20), ...
    struct('tag','tau002',     'P0_damp', c_prior, 'tau_m', 0.02), ...
    struct('tag','delay0',     'P0_damp', c_prior, 'delay_pwm', 0), ...
    struct('tag','nolag',      'P0_damp', c_prior, 'delay_pwm', 0, 'tau_m', 0.02), ...
    struct('tag','lock_prior', 'lock_damp', c_prior), ...
    struct('tag','lock_phys',  'lock_damp', c_phys), ...
    ... % Jz e Jxz LIVRES (hoje travados no CAD): tenta reproduzir a Tabela 5.3,
    ... % que mostra Jz = 0,1453 e Jxz = 0,00261 — só possível com limites abertos.
    struct('tag','freeJ', 'lb_J', [0.03243; 0.06330; 0.09464; 0.00000], ...
                          'ub_J', [0.05406; 0.10551; 0.15774; 0.01000]), ...
    ... % Tentativas de REPRODUZIR a Tabela 5.3 (Jy 0,0984 | Jz 0,1453 = Jx+Jy | Jxz 0,00261):
    ... % oldbounds = limites largos do commit b9e7576 (estado intermediário provável)
    struct('tag','oldbounds', 'lb_J', [0.040; 0.050; 0.090; 0.000], ...
                              'ub_J', [0.080; 0.310; 0.300; 0.010]), ...
    ... % Jz ±25 %, Jxz preso acima do CAD (só pode subir)
    struct('tag','JxzUp',     'lb_J', [0.03243; 0.06330; 0.09464; 0.001571], ...
                              'ub_J', [0.05406; 0.10551; 0.15774; 0.005]), ...
    ... % Jz livre até 0,30 e Jxz ±100 % do CAD
    struct('tag','Jz30',      'lb_J', [0.03243; 0.06330; 0.09464; 0.000], ...
                              'ub_J', [0.05406; 0.10551; 0.30000; 0.00314]) };
tags = cellfun(@(c) c.tag, cases, 'UniformOutput', false);
if nargin < 1 || isempty(which), which = tags; end
if ischar(which), which = {which}; end

fprintf('run_diag: a priori c = [%.3f %.3f %.3f] | rotor+estrutura(%.1f m/s) = [%.3f %.3f %.3f]\n', ...
    c_prior, V_env, c_phys);
for k = 1:numel(cases)
    if ~ismember(cases{k}.tag, which), continue; end
    fprintf('\n################ CASO %s ################\n', cases{k}.tag);
    try
        run_one(cases{k});
    catch ME
        fprintf(2, '!!! caso %s falhou: %s\n', cases{k}.tag, ME.message);
        for e = 1:min(3,numel(ME.stack)), fprintf(2, '    em %s (linha %d)\n', ME.stack(e).name, ME.stack(e).line); end
    end
    close all;
end
diag_report;
end

function run_one(D)
% Workspace limpo por caso: o script roda dentro desta função e enxerga DIAG.
DIAG = D; %#ok<NASGU>
identify_plant_diag;
end
