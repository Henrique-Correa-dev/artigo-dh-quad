% make_cases.m — gera os arquivos de configuração (CFG) das rodadas de identificação
% =========================================================================
% Cada caso vira outputs/runs/<tag>/cfg.mat. Para rodar um caso:
%
%   IDP_CFG_FILE=identification/outputs/runs/<tag>/cfg.mat \
%       matlab -batch "cd identification/3_identification; identify_plant"
%
% O identify_plant lê a variável de ambiente, sobrescreve log/janelas/flags e
% grava tudo em outputs/runs/<tag>/ (P_identified.mat, summary.mat) e
% outputs/images/runs/<tag>/. Sem a variável de ambiente ele roda o caso
% oficial de sempre (voo de 2026, modelo estendido, saída na raiz de outputs).
%
% Uso:  >> make_cases
% =========================================================================
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
AVL = proj.P0_J(16:22);                       % [C_d Cl_p Cl_b Cm_q Cm_a Cn_r Cn_b]
OFF = struct('C_d',false,'Cl_p',false,'Cl_b',false,'Cm_q',false,'Cm_a',false,'Cn_r',false,'Cn_b',false);
ON  = struct('C_d',true, 'Cl_p',true, 'Cl_b',true, 'Cm_q',true, 'Cm_a',true, 'Cn_r',false,'Cn_b',true);

C = {};
% ---- voo de 2026 (dissertação): três tratamentos da aerodinâmica ----------
C{end+1} = struct('tag','aero_off_2026',  'AERO_FREE',OFF, 'P0_aero',zeros(7,1), ...
    'nota','15 parâmetros: aerodinâmica desligada (referência da dissertação)');
C{end+1} = struct('tag','aero_avl_2026',  'AERO_FREE',OFF, 'P0_aero',AVL, ...
    'nota','aerodinâmica FIXA nos valores da AVL (não identificada)');
C{end+1} = struct('tag','aero_free_2026', 'AERO_FREE',ON,  'P0_aero',AVL, ...
    'nota','aerodinâmica identificada (Cn_r travado: não observável em guinada)');
AVL_nocd = AVL;  AVL_nocd(1) = 0;    % só os momentos, sem o arrasto induzido
C{end+1} = struct('tag','aero_avlM_2026', 'AERO_FREE',OFF, 'P0_aero',AVL_nocd, ...
    'nota','momentos aerodinâmicos FIXOS na AVL, C_d = 0 (isola o arrasto)');

% ---- campanha de fev/2024 (motores antigos, voo único de 73 s) -----------
FEV = fullfile(fileparts(paths.root), 'babyshark_vtol_model', 'examples', ...
               'example_inputs', 'Logs do Voo Vertical', 'VooVert_001.mat');
tr_fev  = {[141,151];[151,161];[161,171];[171,181];[181,190]};
val_fev = [191, 211];
C{end+1} = struct('tag','fev2024_noaero', 'LOG_FILE',FEV, 't_trains',{tr_fev}, 't_val',val_fev, ...
    'AERO_FREE',OFF, 'P0_aero',zeros(7,1), 'nota','fev/2024, 15 parâmetros, CG do CAD (sem campanha)');
C{end+1} = struct('tag','fev2024_cg',  'LOG_FILE',FEV, 't_trains',{tr_fev}, 't_val',val_fev, ...
    'AERO_FREE',OFF, 'P0_aero',zeros(7,1), 'campaign','fev2024', 'kT_bounds',[0.4 2.0], ...
    'nota','fev/2024, 15 parâmetros, CG da campanha (estimado pelo equilíbrio)');
% Janelas escolhidas pela SATURAÇÃO: a bateria vai caindo e o motor 2 encosta no
% teto de 1950 µs em fração crescente do voo (0% em 141–151, 78% em 191–196).
% Só os trechos com pouca saturação servem para identificar.
tr_lim  = {[141,151];[161,166];[171,176]};   % 0%, 12% e 4% de saturação
val_lim = [206, 211];                         % 16%, fim do voo, não visto no treino
C{end+1} = struct('tag','fev2024_lim', 'LOG_FILE',FEV, 't_trains',{tr_lim}, 't_val',val_lim, ...
    'AERO_FREE',OFF, 'P0_aero',zeros(7,1), 'campaign','fev2024', 'kT_bounds',[0.4 2.0], ...
    'nota','fev/2024, janelas sem saturação, CG da campanha');
C{end+1} = struct('tag','fev2024_cg_aero', 'LOG_FILE',FEV, 't_trains',{tr_fev}, 't_val',val_fev, ...
    'AERO_FREE',ON,  'P0_aero',AVL, 'campaign','fev2024', 'kT_bounds',[0.4 2.0], ...
    'nota','fev/2024, modelo estendido, CG da campanha');
C{end+1} = struct('tag','fev2024_aero',   'LOG_FILE',FEV, 't_trains',{tr_fev}, 't_val',val_fev, ...
    'AERO_FREE',ON,  'P0_aero',AVL, 'nota','fev/2024, modelo estendido');
% Os motores de fev/2024 não são os de 2026: com a curva de bancada de 2026 o
% empuxo no PWM voado fica ~25% abaixo do peso, e k_T satura em 1,40. Duas saídas:
%   kT   — mesma função de custo, k_T livre até 2,5 (absorve a escala de empuxo)
%   rot  — custo só rotacional: momentos dependem de DIFERENÇAS de empuxo, não da
%          escala absoluta, então a incerteza do motor não contamina p, q, r
C{end+1} = struct('tag','fev2024_kT',  'LOG_FILE',FEV, 't_trains',{tr_fev}, 't_val',val_fev, ...
    'AERO_FREE',OFF, 'P0_aero',zeros(7,1), 'kT_bounds',[0.4 2.5], ...
    'nota','fev/2024, 15 parâmetros, k_T livre até 2,5');
C{end+1} = struct('tag','fev2024_rot', 'LOG_FILE',FEV, 't_trains',{tr_fev}, 't_val',val_fev, ...
    'AERO_FREE',OFF, 'P0_aero',zeros(7,1), 'kT_bounds',[0.4 2.5], 'COST_MODE','rotational', ...
    'nota','fev/2024, custo só rotacional, k_T livre até 2,5');

for i = 1:numel(C)
    CFG = C{i};
    d = fullfile(paths.outputs, 'runs', CFG.tag);
    if ~exist(d,'dir'), mkdir(d); end
    save(fullfile(d,'cfg.mat'), 'CFG');
    fprintf('  %-16s %s\n', CFG.tag, CFG.nota);
end
fprintf('\n  %d configurações em %s\n', numel(C), fullfile(paths.outputs,'runs'));
