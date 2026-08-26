% identify_plant.m - VTOL multi-maneuver identification (15 params)
%
% Localização: 3_identification/
% Substitui o antigo new_identification.m.
%
% Multi-maneuver: usa múltiplas janelas de treino simultaneamente.
% O custo é a soma dos resíduos de todos os segmentos, conforme
% recomendado por Klein & Morelli e Jategaonkar.
%
% Fase A: EEM rotacional (algébrico, dados concatenados) → P_eem
% Fase B: OEM progressivo (1s → 2s → 5s → full, multi-segmento) → P_final
%
% P(1:4)   = inércias [Jx, Jy, Jz, Jxz] → G1-G8, InvJy computados internamente
% P(5:8)   = k_T1..k_T4
% P(9:12)  = k_Q1..k_Q4
% P(13:15) = Dp, Dq, Dr
% (Bp/Bq/Br, Xu/Yv/Zw, Bz REMOVIDOS — modelo é rotacional puro)
% (Xu_m, Yv_m, Zw_m REMOVIDOS — não-identificáveis sem ground truth de u,v,w)
% (dx_cg, dy_cg REMOVIDOS — CG oficial = ponto de referência do CAD)

%% 0. Setup de paths (adiciona subpastas ao MATLAB path)
addpath(fileparts(fileparts(mfilename('fullpath'))));   % raiz do projeto
paths = setup_paths();

%% 1. Configuração de janelas e log
% EDITE AQUI: escolha o log e as janelas de treino/validação.
%
% LOG_FILE:  nome do arquivo dentro de 1_data/
%   - 'log_data.mat'                              → log antigo (struct ATT/IMU/RCOU)
%   - '4 25-05-2026 09-31-48.log-132954.mat'      → log novo (mp_export, motores novos)
%   - '5 25-05-2026 10-27-52.bin-155264.mat'      → log novo
%   - '6 31-12-1979 21-00-00.bin-116760.mat'      → log novo
%   - '7 31-12-1979 21-00-00.bin-113223.mat'      → log novo

LOG_FILE = 'logs_concat.mat';   % ◄── escolher aqui

% ---- CONFIGURAÇÃO EXTERNA (opcional) -------------------------------------
% Se a variável de ambiente IDP_CFG_FILE apontar para um .mat contendo um
% struct CFG, os campos dele SOBRESCREVEM os defaults deste bloco:
%   CFG.LOG_FILE   nome do log em 1_data/
%   CFG.t_trains   cell de janelas de treino     CFG.t_val  janela de validação
%   CFG.AERO_FREE  struct de flags (ver seção 5) CFG.P0_aero  chute de P(16:22)
%   CFG.tag        subpasta de saída em outputs/runs/<tag>/ e images/runs/<tag>/
% Assim dá para rodar outra campanha de voo, ou uma variante do modelo, sem
% editar o script (era o papel do make_diag.py, que gerava um clone do arquivo).
CFG = struct();
cfg_file = getenv('IDP_CFG_FILE');
if ~isempty(cfg_file) && exist(cfg_file, 'file')
    CFG = load(cfg_file);  if isfield(CFG,'CFG'), CFG = CFG.CFG; end
    fprintf('  >>> CFG externo: %s\n', cfg_file);
end
if isfield(CFG,'LOG_FILE'), LOG_FILE = CFG.LOG_FILE; end
% Campanha de voo: fixa a configuração da aeronave (CG, motores) para TODO o
% processo, inclusive nas funções internas que chamam parameters() sem argumento.
% Forma do amortecimento para esta rodada ('moment' | 'rate'), via override do
% parameters(). Sem o campo, vale o default do parameters.m.
if isfield(CFG,'damp_form')
    setappdata(0,'damp_form_override', CFG.damp_form);  clear aero_gains
    fprintf('  >>> Forma do amortecimento: %s\n', CFG.damp_form);
end
if isfield(CFG,'kv_thrust'),  setappdata(0,'kv_thrust', CFG.kv_thrust);  fprintf('  >>> k_v no empuxo: %.3f N/(m/s)\n', CFG.kv_thrust); end
if isfield(CFG,'imu_rx'),     setappdata(0,'imu_rx_override', CFG.imu_rx); fprintf('  >>> r_x da IMU: %.3f m\n', CFG.imu_rx); end
if isfield(CFG,'campaign')
    parameters('campaign', CFG.campaign);
    clear aero_gains          % o cache de aero_gains precisa reler a geometria
    fprintf('  >>> Campanha: %s\n', CFG.campaign);
end

% Janelas (ATENÇÃO: t em segundos, alinhado com TimeUS do log escolhido)
% Pro log "4 25-05-2026 09-31-48" (motores novos, voo entre 429-655 s):
t_trains = {[4, 24]; ...
            [25, 41]; ...
            [42, 62]; ...
            [63, 99]; ...
            [100, 125]};
% Janela de validação: usar 600-645s (voo principal, alta excitação de
% yaw — diagnose_mz confirmou correlação 0.96 nessa janela com Mz_target).
% A antiga 494-514 dava R² enganoso porque tem r oscilando só em ±0.1,
% e erro DC pequeno do Mz vira drift visível integrado open-loop.
t_val    = [610, 630];   % (antes 605–625) — janela de validação oficial
if isfield(CFG,'t_trains'), t_trains = CFG.t_trains(:); end
if isfield(CFG,'t_val'),    t_val    = CFG.t_val(:)';   end
RUN_TAG = '';  if isfield(CFG,'tag'), RUN_TAG = CFG.tag; end

% Janelas pro log antigo (mantidas como referência — comente acima e descomente aqui):
% t_trains = {[157, 167]; [167, 177]; [177, 187]};
% t_val    = [147, 157];

%% 2. Data loading and preprocessing (formato unificado via load_log_data.m)

log_path = LOG_FILE;                                   % aceita caminho absoluto
if ~exist(log_path, 'file'), log_path = fullfile(paths.data, LOG_FILE); end
L = load_log_data(log_path);

% Extrai variáveis do struct uniforme L (compatibilidade com código abaixo)
time_IMU = L.time_IMU;   time_ATT = L.time_ATT;   time_RCOU = L.time_RCOU;
gyrX_raw = L.gyrX_raw;   gyrY_raw = L.gyrY_raw;   gyrZ_raw = L.gyrZ_raw;
accX_raw = L.accX_raw;   accY_raw = L.accY_raw;   accZ_raw = L.accZ_raw;
pwm1_raw = L.pwm1_raw;   pwm2_raw = L.pwm2_raw;
pwm3_raw = L.pwm3_raw;   pwm4_raw = L.pwm4_raw;

% Grade comum de tempo (dt=0.1 s, igual antes)
% NÃO incluir GPS aqui: ele não é usado na identificação rotacional (pqr/acc)
% e no logs_concat o GPS termina em 268s enquanto IMU/ATT/RCOU vão até ~749s.
% Incluí-lo truncava a grade em 268s e barrava janelas válidas (ex.: 455-465).
all_starts = [min(time_IMU), min(time_ATT), min(time_RCOU)];
all_ends   = [max(time_IMU), max(time_ATT), max(time_RCOU)];
t_start  = max(all_starts);
t_end    = min(all_ends);
t_common = t_start:0.1:t_end;
dt = 0.1;

gyrX_interp = interp1(time_IMU, gyrX_raw, t_common, 'linear');
gyrY_interp = interp1(time_IMU, gyrY_raw, t_common, 'linear');
gyrZ_interp = interp1(time_IMU, gyrZ_raw, t_common, 'linear');
accX_interp = interp1(time_IMU, accX_raw, t_common, 'linear');
accY_interp = interp1(time_IMU, accY_raw, t_common, 'linear');
accZ_interp = interp1(time_IMU, accZ_raw, t_common, 'linear');
roll_interp  = interp1(time_ATT, L.roll_deg,  t_common, 'linear');
pitch_interp = interp1(time_ATT, L.pitch_deg, t_common, 'linear');
yaw_interp   = interp1(time_ATT, L.yaw_deg,   t_common, 'linear');
pwm1_interp = interp1(time_RCOU, pwm1_raw, t_common, 'linear');
pwm2_interp = interp1(time_RCOU, pwm2_raw, t_common, 'linear');
pwm3_interp = interp1(time_RCOU, pwm3_raw, t_common, 'linear');
pwm4_interp = interp1(time_RCOU, pwm4_raw, t_common, 'linear');

% ---- CADEIA DE ATUAÇÃO ÚNICA (2_model/motor_chain.m) ----------------------
% PWM(log) → atraso puro (parameters().motor.delay_s) → RPM_ss = CR·PWM+Ω_b →
% lag de 1ª ordem τ_m (solução exata p/ entrada linear por partes) → RPM.
% Daqui em diante os arrays "pwm*_interp" passam a conter RPM já atrasada e
% lagada, e func_T_ref/func_Q_ref agem sobre RPM (T = kT_rpm·RPM²). É a mesma
% cadeia usada em sim_window / results_id_figures / compare_*, portanto a
% validação vê exatamente o mesmo motor que a identificação.
% (Histórico: DELAY_PWM=1 + filtro discreto a 10 Hz + empuxo constante por passo
%  → substituído; ver identify_plant_prelag_backup.m.txt e outputs/diag/*.)
[rpm_lag, func_T_ref, func_Q_ref, mc_info] = motor_chain(t_common(:), ...
    [pwm1_interp(:), pwm2_interp(:), pwm3_interp(:), pwm4_interp(:)]);
pwm1_interp = rpm_lag(:,1)';  pwm2_interp = rpm_lag(:,2)';
pwm3_interp = rpm_lag(:,3)';  pwm4_interp = rpm_lag(:,4)';
fprintf('  >>> motor_chain: atraso %.2f s (%d amostra) | lag τ_m = %.3f s | T = kT·RPM²\n', ...
    mc_info.delay_s, mc_info.delay_n, mc_info.tau_m);
DELAY_PWM = mc_info.delay_n;   % mantido só para logs/plots que o citam
USE_MOTOR_MODEL = true;        % idem (compatibilidade)

% ---- VELOCIDADE ESTIMADA (escala a aerodinâmica da estrutura, P(17:22)) ----
% Sem GPS utilizável no ginásio: EKF VD (barômetro) + integração ancorada da
% aceleração quase-estática de atitude e empuxo (5_validation/estimate_velocity.m).
% É EXÓGENA (não é estado do modelo), e a MESMA série vai para o custo EEM, para
% o custo OEM e para a simulação de validação (appdata 'aero_vel' lido pelo
% vtol_dynamics), de modo que identificação e validação veem o mesmo V(t).
VEst = estimate_velocity(L, t_common(:));
V_grid = VEst.V;  v_grid = VEst.v;  w_grid = VEst.w;
V_grid(~isfinite(V_grid)) = 0;  v_grid(~isfinite(v_grid)) = 0;  w_grid(~isfinite(w_grid)) = 0;
setappdata(0, 'aero_vel', struct('t', t_common(:), 'V', V_grid, 'v', v_grid, 'w', w_grid));
fprintf('  >>> V estimada: média %.2f | p95 %.2f | máx %.2f m/s (escala P(17:22))\n', ...
    mean(V_grid), prctile(V_grid,95), max(V_grid));

% Sanity check: as janelas escolhidas estão dentro do log?
log_t_min = t_common(1);  log_t_max = t_common(end);
for s = 1:numel(t_trains)
    if t_trains{s}(1) < log_t_min || t_trains{s}(2) > log_t_max
        error('identify_plant:badWindow', ...
            'Janela treino %d [%g, %g] fora do log [%.1f, %.1f]', ...
            s, t_trains{s}(1), t_trains{s}(2), log_t_min, log_t_max);
    end
end
if t_val(1) < log_t_min || t_val(2) > log_t_max
    error('identify_plant:badWindow', 'Janela validação fora do log.');
end
fprintf('  Log "%s": t = %.1f a %.1f s (%s)\n', ...
    LOG_FILE, log_t_min, log_t_max, L.format);

%% 3. Motor reference models → definidos por motor_chain (T = kT_rpm·RPM², Q = kQ_rpm·RPM²)

%% 4. Extract training segments
n_seg = length(t_trains);
segs = cell(n_seg, 1);

fprintf('\n  Segmentos de treino: %d\n', n_seg);
for s = 1:n_seg
    seg = struct();   % zera (evita herdar tipo de variável 'seg' do workspace, ex.: cell do identify_v2)
    idx_s = (t_common >= t_trains{s}(1)) & (t_common <= t_trains{s}(2));
    seg.time = t_common(idx_s)';
    seg.N    = sum(idx_s);
    seg.pwm  = [pwm1_interp(idx_s)', pwm2_interp(idx_s)', ...
                pwm3_interp(idx_s)', pwm4_interp(idx_s)'];
    seg.pqr  = [gyrX_interp(idx_s)', gyrY_interp(idx_s)', gyrZ_interp(idx_s)'];
    seg.acc  = [accX_interp(idx_s)', accY_interp(idx_s)', accZ_interp(idx_s)'];
    seg.att  = [roll_interp(idx_s)', pitch_interp(idx_s)', yaw_interp(idx_s)'];
    seg.att_rad = deg2rad(seg.att);
    seg.V = V_grid(idx_s);  seg.v = v_grid(idx_s);  seg.w = w_grid(idx_s);
    seg.T_ref = zeros(seg.N, 4);
    seg.Q_ref = zeros(seg.N, 4);
    for j = 1:4
        seg.T_ref(:,j) = func_T_ref(seg.pwm(:,j));
        seg.Q_ref(:,j) = func_Q_ref(seg.pwm(:,j));
    end
    segs{s} = seg;
    fprintf('    [%d] %d-%ds  (%d pontos)\n', s, t_trains{s}(1), t_trains{s}(2), seg.N);
end

% Validation data
idx_val = (t_common >= t_val(1)) & (t_common <= t_val(2));
time_vl = t_common(idx_val)';
pwm_vl  = [pwm1_interp(idx_val)', pwm2_interp(idx_val)', ...
            pwm3_interp(idx_val)', pwm4_interp(idx_val)'];
pqr_vl  = [gyrX_interp(idx_val)', gyrY_interp(idx_val)', gyrZ_interp(idx_val)'];
acc_vl  = [accX_interp(idx_val)', accY_interp(idx_val)', accZ_interp(idx_val)'];
att_vl  = [roll_interp(idx_val)', pitch_interp(idx_val)', yaw_interp(idx_val)'];

%% 5. Parameters: initial guess and physical bounds (22 params)
%  (dx_cg, dy_cg REMOVIDOS — CG oficial é o ponto de referência do CAD)
% P(1:4)   = inércias [Jx, Jy, Jz, Jxz]
% P(5:8)   = k_T1..k_T4
% P(9:12)  = k_Q1..k_Q4
% P(13:15) = Dp, Dq, Dr
% (Bp/Bq/Br, Xu/Yv/Zw, Bz REMOVIDOS — modelo é rotacional puro)
% (translacionais drag removidos — só Bz)

% Vem de parameters.m (fonte única)
proj_params = parameters();
P0 = proj_params.P0_J;
lb = proj_params.bounds.lb;
ub = proj_params.bounds.ub;
param_names = proj_params.param_names;

n_params = numel(P0);          % 22 no modelo estendido (15 rotacionais + 7 aerodinâmicos)
n_rot    = proj_params.n_rot_only;   % 15 — bloco puramente rotacional (J, k_T, k_Q, D)

% ---- QUAIS COEFICIENTES AERODINÂMICOS FICAM LIVRES --------------------------
% Decidido pelo 3_identification/identifiability_check.m (Jacobiano em Θ₀, EEM):
% no envelope voado (V p95 ≈ 1,3 m/s) as cotas de Cramér-Rao em Θ₀ dão
%   Cl_p 33% | Cl_β 61% | Cm_q 58% | Cm_α 64% | Cn_r 306% | Cn_β 21%
% e as correlações rotor↔estrutura ρ(D_p,Cl_p)=0,53, ρ(D_q,Cm_q)=0,79,
% ρ(D_r,Cn_r)=0,85. Ou seja: guinada não separa amortecimento de rotor do
% aerodinâmico (canal fraco, r pequeno), então Cn_r fica TRAVADO no valor da AVL.
% Os demais ficam livres. C_d entra pelo canal translacional (acelerômetro).
% Ponha todos em true para o teste "aerodinâmica totalmente livre".
% DEFAULT: todos TRAVADOS nos valores da AVL. Liberá-los foi testado e é
% sobreajuste (runs/aero_free_2026: rolagem sobe para 0,80 e arfagem cai para
% 0,52, com C_lp e C_mq indo a zero e C_mα a 4× a AVL). Travados, eles entram
% como física sem grau de liberdade novo e melhoram a guinada de 0,73 para 0,81.
AERO_FREE = struct('C_d', false, 'Cl_p', false, 'Cl_b', false, ...
                   'Cm_q', false, 'Cm_a', false, 'Cn_r', false, 'Cn_b', false);
if isfield(CFG,'AERO_FREE'), AERO_FREE = CFG.AERO_FREE; end
if isfield(CFG,'P0_aero'),   P0(16:22) = CFG.P0_aero(:); end
% Rotor, momento por velocidade P(23:25) = L_v, M_w, N_v: default TRAVADOS em 0.
% CFG.rotvel_free = true libera (teste NASA C_lv, C_mw, C_nv).
if ~(isfield(CFG,'rotvel_free') && CFG.rotvel_free)
    lb(23:25) = P0(23:25);  ub(23:25) = P0(23:25);
end
% Forma do Hu et al. (2024): P(23:25) = c_l, c_m, c_n, termo −½ρV²S·(b,c̄,b)·c.
% Substitui os seis coeficientes de asa P(17:22), que ficam zerados.
if isfield(CFG,'hu_form') && CFG.hu_form
    setappdata(0,'hu_form',true);
    P0(17:22) = 0;  lb(17:22) = 0;  ub(17:22) = 0;
    lb(23:25) = -2;  ub(23:25) = 2;
    fprintf('  >>> Forma do Hu: c_l, c_m, c_n livres; asa P(17:22) zerada\n');
else
    if isappdata(0,'hu_form'), rmappdata(0,'hu_form'); end
end
aero_map = {'C_d',16; 'Cl_p',17; 'Cl_b',18; 'Cm_q',19; 'Cm_a',20; 'Cn_r',21; 'Cn_b',22};
for a = 1:size(aero_map,1)
    ia = aero_map{a,2};
    if ~AERO_FREE.(aero_map{a,1}), lb(ia) = P0(ia);  ub(ia) = P0(ia); end
end
fprintf('  >>> Aerodinâmicos livres: %s | travados: %s\n', ...
    strjoin(aero_map(cellfun(@(n) AERO_FREE.(n), aero_map(:,1)),1), ', '), ...
    strjoin(aero_map(cellfun(@(n) ~AERO_FREE.(n), aero_map(:,1)),1), ', '));

% FLAGS de "trava": lb==ub==valor → lsqnonlin não otimiza o parâmetro.
LOCK_MOTORS = false;   % true: trava k_T/k_Q=1.0 | false: livre [0.40,1.40]

% Regularização de pares CW/CCW (identifiability): yaw só vê Q1+Q2-Q3-Q4, logo
% k_Q individual não é único → penaliza a diferença dentro do par (M1+M2, M3+M4).
% Guia de λ: 1≈0%, 100≈50%, 1000≈trava do peso do data fit.
REG_LAMBDA.kt_pair = 100;   % penaliza |k_T1-k_T2| e |k_T3-k_T4|
REG_LAMBDA.kq_pair = 100;   % penaliza |k_Q1-k_Q2| e |k_Q3-k_Q4|
% Restrição física: desigualdade triangular dos momentos de inércia
%   (Jz≤Jx+Jy, Jy≤Jx+Jz, Jx≤Jy+Jz). One-sided: só pune violação.
%   Necessária porque a folga triangular do CAD é ~1.2% e yaw é inobservável.
REG_LAMBDA.tri = 1e4;       % grande → optimizer nunca sai do físico
if isfield(CFG,'tri_lambda'), REG_LAMBDA.tri = CFG.tri_lambda; fprintf('  >>> tri_lambda (CFG): %.1e\n', REG_LAMBDA.tri); end

% Modo do custo (quais sinais no resíduo):
%   'rotational' → só p,q,r (k_T solto em magnitude → AccZ diverge ~20% no hover)
%   'translational' → só acc (calibra IMU offset isolado)
%   'full' → tudo (default); k_T ancorado por AccZ≈-g no hover.
COST_MODE = 'full';
if isfield(CFG,'COST_MODE'), COST_MODE = CFG.COST_MODE; end
% Limites de k_T por campanha: os motores de fev/2024 não são os de 2026 e a
% curva de bancada disponível é a dos novos, logo k_T precisa de faixa maior
% para absorver a diferença de escala de empuxo.
if isfield(CFG,'kT_bounds'), lb(5:8) = CFG.kT_bounds(1); ub(5:8) = CFG.kT_bounds(2); end
if isfield(CFG,'kQ_bounds'), lb(9:12) = CFG.kQ_bounds(1); ub(9:12) = CFG.kQ_bounds(2); end
% Amortecimento: chute e trava. Serve para testar se os valores FÍSICOS
% (2_model/damping_budget.m) sustentam a validação tão bem quanto os efetivos.
if isfield(CFG,'P0_damp'), P0(13:15) = CFG.P0_damp(:); end
if isfield(CFG,'lock_damp') && CFG.lock_damp
    lb(13:15) = P0(13:15);  ub(13:15) = P0(13:15);
    fprintf('  >>> c_p, c_q, c_r TRAVADOS em [%.3f %.3f %.3f]\n', P0(13:15));
end
% Limites dos aerodinâmicos: quando c_p, c_q, c_r são zerados, os termos ∝ V
% precisam de faixa muito maior para tentar substituí-los.
% J_xz livre: ele é quem dita Γ4, e portanto a intensidade do acoplamento
% cruzado da forma de momento. Travado no CAD ele tem CR de 74%.
% Limites do amortecimento: dependem da FORMA (N·m·s na de momento, 1/s na de taxa)
if isfield(CFG,'damp_bounds')
    lb(13:15) = CFG.damp_bounds(:,1);  ub(13:15) = CFG.damp_bounds(:,2);
end
if isfield(CFG,'Jxz_bounds'), lb(4) = CFG.Jxz_bounds(1); ub(4) = CFG.Jxz_bounds(2); end
if isfield(CFG,'Jz_bounds'),  lb(3) = CFG.Jz_bounds(1);  ub(3) = CFG.Jz_bounds(2);  end
if isfield(CFG,'aero_bounds')
    lb(16:22) = CFG.aero_bounds(:,1);  ub(16:22) = CFG.aero_bounds(:,2);
end
if isfield(CFG,'J_bounds')                       % fator multiplicativo sobre o CAD
    lb(1:2) = proj_params.P0_J(1:2)*CFG.J_bounds(1);
    ub(1:2) = proj_params.P0_J(1:2)*CFG.J_bounds(2);
end

if LOCK_MOTORS
    P0(5:12) = 1.0;
    lb(5:12) = 1.0;
    ub(5:12) = 1.0;
    fprintf('  >>> k_T(1:4), k_Q(1:4) TRAVADOS em 1.0 (motores idênticos)\n');
end
fprintf('  >>> Regularização: λ_kT_pair=%.1e  λ_kQ_pair=%.1e\n', ...
    REG_LAMBDA.kt_pair, REG_LAMBDA.kq_pair);

% NOTA: pra TESTAR parâmetros específicos SEM rodar otimização, use o
% script separado: 5_validation/validate_params.m

R2_func = @(y_e, y_s) 1 - sum((y_e - y_s).^2) / max(sum((y_e - mean(y_e)).^2), 1e-12);
constants_sim.m = proj_params.m;
constants_sim.g = proj_params.g;
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);  % Mesmo que Simulink

% Criar pasta de imagens se não existir (centralizado em outputs/images/)
img_dir = paths.images;
out_dir = paths.outputs;
if ~isempty(RUN_TAG)
    img_dir = fullfile(paths.images, 'runs', RUN_TAG);
    out_dir = fullfile(paths.outputs, 'runs', RUN_TAG);
end
if ~exist(img_dir, 'dir'), mkdir(img_dir); end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ========================================================================
%  6. VALIDAÇÃO COM CHUTE INICIAL (P0) - primeiro segmento + validação
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  VALIDAÇÃO DO MODELO COM P0 (chute inicial)\n');
fprintf('==========================================================\n');

seg1 = segs{1};

% === DIAGNÓSTICO: comparar valores com Simulink ===
diag_P = P0;
diag_Jx = diag_P(1); diag_Jy = diag_P(2); diag_Jz = diag_P(3); diag_Jxz = diag_P(4);
diag_gamma = diag_Jx*diag_Jz - diag_Jxz^2;
fprintf('\n  --- Constantes do modelo (compare com Simulink) ---\n');
fprintf('  Jx=%.6f  Jy=%.6f  Jz=%.6f  Jxz=%.6f\n', diag_Jx, diag_Jy, diag_Jz, diag_Jxz);
fprintf('  gamma0 = %.6f\n', diag_gamma);
fprintf('  G1=%.6f  G2=%.6f  G3=%.6f  G4=%.6f\n', ...
    diag_Jxz*(diag_Jx-diag_Jy+diag_Jz)/diag_gamma, ...
    (diag_Jz*(diag_Jz-diag_Jy)+diag_Jxz^2)/diag_gamma, ...
    diag_Jz/diag_gamma, diag_Jxz/diag_gamma);
fprintf('  G5=%.6f  G6=%.6f  G7=%.6f  G8=%.6f\n', ...
    (diag_Jz-diag_Jx)/diag_Jy, diag_Jxz/diag_Jy, ...
    (diag_Jx*(diag_Jx-diag_Jy)+diag_Jxz^2)/diag_gamma, diag_Jx/diag_gamma);
fprintf('  invJy = %.6f\n', 1/diag_Jy);
fprintf('  k_T = [%.4f, %.4f, %.4f, %.4f]\n', diag_P(5:8));
fprintf('  k_Q = [%.4f, %.4f, %.4f, %.4f]\n', diag_P(9:12));
fprintf('  Dp=%.4f  Dq=%.4f  Dr=%.4f\n', diag_P(13:15));
% (Bp, Bq, Br removidos do modelo)
fprintf('  Braços (de parameters.m): Lx_r=%.6f  Lx_l=%.6f  Ly_f=%.6f  Ly_r=%.6f\n', ...
    proj_params.arms.Lx_r, proj_params.arms.Lx_l, ...
    proj_params.arms.Ly_f,    proj_params.arms.Ly_r);

% Derivadas no instante t=0 (para comparação pontual com Simulink)
t0 = seg1.time(1);
y0_rot = seg1.pqr(1,:)';
dy0 = vtol_dynamics(t0, y0_rot, P0, seg1.time, seg1.pwm, func_T_ref, func_Q_ref, constants_sim);
pwm0 = seg1.pwm(1,:);
T0 = diag_P(5:8)' .* [func_T_ref(pwm0(1)), func_T_ref(pwm0(2)), func_T_ref(pwm0(3)), func_T_ref(pwm0(4))];
Q0 = diag_P(9:12)' .* [func_Q_ref(pwm0(1)), func_Q_ref(pwm0(2)), func_Q_ref(pwm0(3)), func_Q_ref(pwm0(4))];
Lx_r = proj_params.arms.Lx_r; Lx_l = proj_params.arms.Lx_l;
Ly_f = proj_params.arms.Ly_f;    Ly_r = proj_params.arms.Ly_r;
Mx0 = -(Lx_r*T0(1) - Lx_l*T0(2) - Lx_l*T0(3) + Lx_r*T0(4));
% My — ArduPilot QuadX (M1,M3=FRONT; M2,M4=REAR)
My0 = Ly_f*T0(1) - Ly_r*T0(2) + Ly_f*T0(3) - Ly_r*T0(4);
Mz0 = Q0(1) + Q0(2) - Q0(3) - Q0(4);
fprintf('\n  --- Estado inicial (t=%.1fs) ---\n', t0);
fprintf('  PWM = [%.0f, %.0f, %.0f, %.0f] us\n', pwm0);
fprintf('  T_ref(PWM) = [%.4f, %.4f, %.4f, %.4f] N\n', ...
    func_T_ref(pwm0(1)), func_T_ref(pwm0(2)), func_T_ref(pwm0(3)), func_T_ref(pwm0(4)));
fprintf('  T_mr (k_T*T_ref) = [%.4f, %.4f, %.4f, %.4f] N\n', T0);
fprintf('  Q_mr (k_Q*Q_ref) = [%.6f, %.6f, %.6f, %.6f] N·m\n', Q0);
fprintf('  Mx=%.6f  My=%.6f  Mz=%.6f  N·m\n', Mx0, My0, Mz0);
fprintf('  p0=%.4f  q0=%.4f  r0=%.4f  rad/s\n', y0_rot);
fprintf('  p_dot=%.6f  q_dot=%.6f  r_dot=%.6f  rad/s²\n', dy0);
fprintf('  -----------------------------------------------\n\n');

[res_P0] = simulate_full(P0, seg1.time, seg1.pwm, seg1.pqr, seg1.acc, seg1.att, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, ode_opts);

print_R2('P0', res_P0, seg1.pqr, pqr_vl, seg1.acc, acc_vl, t_trains{1}, t_val, R2_func);
plot_all_results('CHUTE INICIAL (P0)', res_P0, seg1.time, seg1.pqr, seg1.acc, ...
    time_vl, pqr_vl, acc_vl, R2_func, t_trains{1}, t_val, att_vl);

%% ========================================================================
%  7. FASE A: EEM rotacional (dados concatenados de todos os segmentos)
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  FASE A: EEM rotacional (%d segmentos, %d params)\n', n_seg, n_rot);
fprintf('==========================================================\n');

smooth_win = 5;
p_all = []; q_all = []; r_all = [];
pd_all = []; qd_all = []; rd_all = [];
Tr_all = []; Qr_all = [];
V_all = []; vb_all = []; wb_all = [];

for s = 1:n_seg
    sg = segs{s};
    V_all = [V_all; sg.V(:)]; vb_all = [vb_all; sg.v(:)]; wb_all = [wb_all; sg.w(:)]; %#ok<AGROW>
    p_s = sg.pqr(:,1); q_s = sg.pqr(:,2); r_s = sg.pqr(:,3);
    p_all = [p_all; p_s]; %#ok<AGROW>
    q_all = [q_all; q_s]; %#ok<AGROW>
    r_all = [r_all; r_s]; %#ok<AGROW>
    pd_all = [pd_all; gradient(movmean(p_s, smooth_win), dt)]; %#ok<AGROW>
    qd_all = [qd_all; gradient(movmean(q_s, smooth_win), dt)]; %#ok<AGROW>
    rd_all = [rd_all; gradient(movmean(r_s, smooth_win), dt)]; %#ok<AGROW>
    Tr_all = [Tr_all; sg.T_ref]; %#ok<AGROW>
    Qr_all = [Qr_all; sg.Q_ref]; %#ok<AGROW>
end

var_pd = var(pd_all); if var_pd < 1e-12, var_pd = 1; end
var_qd = var(qd_all); if var_qd < 1e-12, var_qd = 1; end
var_rd = var(rd_all); if var_rd < 1e-12, var_rd = 1; end
weights_eem = [1/var_pd; 1/var_qd; 1/var_rd];

aero_eem = struct('V', V_all, 'v', vb_all, 'w', wb_all);
cost_eem = @(Pv) eem_cost_function(ones(n_params,1), Pv, weights_eem, ...
    p_all, q_all, r_all, pd_all, qd_all, rd_all, Tr_all, Qr_all, REG_LAMBDA, aero_eem);

% O EEM vê momentos, não forças: C_d (P16) fica travado nesta fase e só é
% identificado no OEM, pelo resíduo do acelerômetro.
lb_eem = lb;  ub_eem = ub;  lb_eem(16) = P0(16);  ub_eem(16) = P0(16);

opts_eem = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter', ...
    'MaxIterations', 1500, ...
    'MaxFunctionEvaluations', 60000, ...
    'StepTolerance', 1e-14, ...
    'FunctionTolerance', 1e-14);

[P_eem, rn_eem, ~, ef_eem] = lsqnonlin(cost_eem, P0, lb_eem, ub_eem, opts_eem);
P_eem_rot = P_eem;   % (nome legado usado nos prints abaixo)

fprintf('\n  EEM Resnorm: %.4f  |  Exit flag: %d\n', rn_eem, ef_eem);

% Breakdown reg vs data — pra calibrar λ
e_eem_final = cost_eem(P_eem_rot);
n_reg = 7;  % 4 pares (kt12,kt34,kq12,kq34) + 3 triangulares (inércia)
rn_data = sum(e_eem_final(1:end-n_reg).^2);
rn_reg  = sum(e_eem_final(end-n_reg+1:end).^2);
fprintf('  Breakdown: data=%.4f  reg=%.4f  (reg/data = %.2f%%)\n', ...
    rn_data, rn_reg, 100*rn_reg/max(rn_data, 1e-12));
fprintf('  k_T par: [%.3f,%.3f] vs [%.3f,%.3f]   Δ=%.3f,%.3f\n', ...
    P_eem_rot(5), P_eem_rot(6), P_eem_rot(7), P_eem_rot(8), ...
    P_eem_rot(5)-P_eem_rot(6), P_eem_rot(7)-P_eem_rot(8));
fprintf('  k_Q par: [%.3f,%.3f] vs [%.3f,%.3f]   Δ=%.3f,%.3f\n', ...
    P_eem_rot(9), P_eem_rot(10), P_eem_rot(11), P_eem_rot(12), ...
    P_eem_rot(9)-P_eem_rot(10), P_eem_rot(11)-P_eem_rot(12));
for i = 1:n_params
    fprintf('    %-6s: %12.6f  (chute: %12.6f)\n', param_names{i}, P_eem(i), P0(i));
end

%% ========================================================================
%  8. FASE B: OEM Progressivo Multi-Segmento (1s → 2s → 3s → 5s → 10s → 20s)
%  ========================================================================

% Pesos OEM: computados sobre dados concatenados de todos os segmentos
pqr_cat = cell2mat(cellfun(@(s) s.pqr, segs, 'UniformOutput', false));
acc_cat = cell2mat(cellfun(@(s) s.acc, segs, 'UniformOutput', false));

var_p = var(pqr_cat(:,1)); if var_p < 1e-12, var_p = 1; end
var_q = var(pqr_cat(:,2)); if var_q < 1e-12, var_q = 1; end
var_r = var(pqr_cat(:,3)); if var_r < 1e-12, var_r = 1; end
weights_pqr = [1/var_p; 1/var_q; 1/var_r];

var_ax = var(acc_cat(:,1)); if var_ax < 1e-12, var_ax = 1; end
var_ay = var(acc_cat(:,2)); if var_ay < 1e-12, var_ay = 1; end
var_az = var(acc_cat(:,3)); if var_az < 1e-12, var_az = 1; end
weights_acc = [1/var_ax; 1/var_ay; 1/var_az];

% Multijanelamento (multiple shooting): cada estágio reinicia a integração a
% cada win_sec segundos, usando o estado medido. Janela curta perdoa erro de
% fase e pesa o transitório; janela longa cobra o comportamento acumulado e é
% mais sensível a erro de ganho e a viés. Percorrer de 1 a 20 s dá vistas
% diferentes do mesmo dado, e o estágio vencedor é escolhido pelo R² médio.
win_durations = [1.0, 2.0, 3.0, 5.0, 10.0, 20.0];
if isfield(CFG,'win_durations'), win_durations = CFG.win_durations(:)'; end
P_current = P_eem;
best_R2_mean = -Inf;
P_best = P_eem;
best_stage_name = 'EEM';

for stage = 1:length(win_durations)
    win_sec = win_durations(stage);
    if isinf(win_sec)
        stage_name = 'FULL';
    else
        stage_name = sprintf('%.0fs', win_sec);
    end

    fprintf('\n==========================================================\n');
    fprintf('  FASE B.%d: OEM [%s] (%d segmentos)\n', stage, stage_name, n_seg);
    fprintf('==========================================================\n');

    cost_oem = @(P) oem_multi_seg_cost(P, segs, win_sec, ...
        constants_sim.m, constants_sim.g, dt, weights_pqr, weights_acc, ...
        COST_MODE, REG_LAMBDA);

    max_iter = 500;
    if isinf(win_sec), max_iter = 1000; end

    opts_oem = optimoptions('lsqnonlin', ...
        'Algorithm', 'trust-region-reflective', ...
        'Display', 'iter', ...
        'MaxIterations', max_iter, ...
        'MaxFunctionEvaluations', 80000, ...
        'StepTolerance', 1e-14, ...
        'FunctionTolerance', 1e-14);

    [P_current, rn, ~, ef] = lsqnonlin(cost_oem, P_current, lb, ub, opts_oem);
    fprintf('  [%s] Resnorm: %.4f | Exit: %d\n', stage_name, rn, ef);

    % Avaliar validação para escolher o melhor estágio
    sg_eval = segs{1};
    res_eval = simulate_full(P_current, sg_eval.time, sg_eval.pwm, sg_eval.pqr, sg_eval.acc, sg_eval.att, ...
        time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
        func_T_ref, func_Q_ref, constants_sim, ode_opts);

    if res_eval.pqr_vl_ok
        R2_p = R2_func(pqr_vl(:,1), res_eval.p_s_vl);
        R2_q = R2_func(pqr_vl(:,2), res_eval.q_s_vl);
        R2_r = R2_func(pqr_vl(:,3), res_eval.r_s_vl);
        R2_mean_stage = mean([R2_p, R2_q, R2_r]);
        fprintf('  [%s] Val R² p=%.4f | q=%.4f | r=%.4f | média=%.4f\n', ...
            stage_name, R2_p, R2_q, R2_r, R2_mean_stage);

        if R2_mean_stage > best_R2_mean
            best_R2_mean = R2_mean_stage;
            P_best = P_current;
            best_stage_name = stage_name;
        end
    else
        fprintf('  [%s] Val -> DIVERGIU\n', stage_name);
    end
end

fprintf('\n  >>> Melhor estágio: [%s] com R² médio = %.4f\n', best_stage_name, best_R2_mean);

P_final = P_best;

fprintf('\n  --- Parâmetros Finais ---\n');
for i = 1:n_params
    fprintf('    %-6s: %12.6f  (EEM: %12.6f)  (P0: %12.6f)\n', ...
        param_names{i}, P_final(i), P_eem(i), P0(i));
end

% Mostrar G's resultantes das inércias identificadas
Jx_f = P_final(1); Jy_f = P_final(2); Jz_f = P_final(3); Jxz_f = P_final(4);
gam0_f = Jx_f*Jz_f - Jxz_f^2;
fprintf('\n  --- Constantes G (derivadas das inércias) ---\n');
fprintf('    G1=%8.4f  G2=%8.4f  G3=%8.4f  G4=%8.4f\n', ...
    Jxz_f*(Jx_f-Jy_f+Jz_f)/gam0_f, (Jz_f*(Jz_f-Jy_f)+Jxz_f^2)/gam0_f, ...
    Jz_f/gam0_f, Jxz_f/gam0_f);
fprintf('    G5=%8.4f  G6=%8.4f  G7=%8.4f  G8=%8.4f  InvJy=%8.4f\n', ...
    (Jz_f-Jx_f)/Jy_f, Jxz_f/Jy_f, (Jx_f*(Jx_f-Jy_f)+Jxz_f^2)/gam0_f, ...
    Jx_f/gam0_f, 1/Jy_f);

%% ========================================================================
%  8a. QUALIDADE DOS PARÂMETROS — erro-padrão via Jacobiano (≈ Cramér-Rao)
%  Cov(P) ≈ σ²·(JᵀJ)⁻¹,  σ² = resnorm/(N−p).  std% alto = parâmetro pouco
%  observável (ex.: Jz com yaw fraco). Fecha a história da identifiabilidade.
%  ========================================================================
try
    cost_se = @(P) oem_multi_seg_cost(P, segs, Inf, ...
        constants_sim.m, constants_sim.g, dt, weights_pqr, weights_acc, ...
        COST_MODE, REG_LAMBDA);
    opts_se = optimoptions('lsqnonlin', 'Display','off', ...
        'MaxIterations',0, 'MaxFunctionEvaluations',1);
    [~, rn_se, res_se, ~, ~, ~, Jse] = lsqnonlin(cost_se, P_final, lb, ub, opts_se);
    Nres  = numel(res_se);
    dof   = max(Nres - n_params, 1);
    sig2  = rn_se / dof;
    Cov   = sig2 * pinv(full(Jse' * Jse));
    se    = sqrt(max(diag(Cov), 0));
    at_lb = abs(P_final(:) - lb(:)) <= 1e-6*max(abs(lb(:)),1);
    at_ub = abs(P_final(:) - ub(:)) <= 1e-6*max(abs(ub(:)),1);
    fprintf('\n  --- Qualidade dos parâmetros (erro-padrão ≈ Cramér-Rao) ---\n');
    fprintf('  %-6s %12s %12s %9s   obs\n', 'param', 'valor', 'std', 'std%%');
    for i = 1:n_params
        relse = 100*se(i) / max(abs(P_final(i)), 1e-9);
        if     at_ub(i), tag = '[no UB]';
        elseif at_lb(i), tag = '[no LB]';
        elseif relse > 50, tag = '<< MAL identificado';
        elseif relse > 20, tag = '[CRB>20%% — Sato FALHA]';
        else,  tag = '[CRB ok]';     % Sato: erro-padrao < 20%%
        end
        fprintf('  %-6s %12.6f %12.6f %8.1f%%   %s\n', ...
            param_names{i}, P_final(i), se(i), relse, tag);
    end
    fprintf('  (Sato: CRB<20%% = bem identificado. std%% alto/"no LB/UB" = pouco\n');
    fprintf('   observável — tipicamente yaw→Jz/k_Q/Jxz/Dr)\n');

    % --- Viabilidade física: desigualdade triangular dos momentos de inércia ---
    %  Jz<=Jx+Jy (e permutações) DEVE valer para QUALQUER corpo rígido. Se violar,
    %  o "modelo" não corresponde a nenhum tensor de inércia físico.
    Jxv=P_final(1); Jyv=P_final(2); Jzv=P_final(3);
    tri = [Jzv, Jxv+Jyv; Jyv, Jxv+Jzv; Jxv, Jyv+Jzv];
    nm  = {'Jz<=Jx+Jy','Jy<=Jx+Jz','Jx<=Jy+Jz'};
    fprintf('\n  --- Viabilidade física (desigualdade triangular) ---\n');
    for i = 1:3
        folga = 100*(tri(i,2)-tri(i,1))/max(tri(i,1),1e-9);
        if tri(i,1) <= tri(i,2)
            fprintf('   %-11s  %.4f <= %.4f   folga %+5.1f%%   OK\n', ...
                nm{i}, tri(i,1), tri(i,2), folga);
        else
            fprintf('   %-11s  %.4f  > %.4f   folga %+5.1f%%   VIOLA <<<\n', ...
                nm{i}, tri(i,1), tri(i,2), folga);
        end
    end
catch ME_se
    fprintf('\n  (erro-padrão não calculado: %s)\n', ME_se.message);
end

%% ========================================================================
%  8b. EXPORTAR PARÂMETROS PARA simulate.m
%  ========================================================================
params_file = fullfile(out_dir, 'P_identified.mat');
save(params_file, 'P_final', 'P0', 'param_names', 'P_eem', 'LOG_FILE', 't_trains', 't_val', 'AERO_FREE');
fprintf('\n  Parâmetros exportados para: %s\n', params_file);

%% ========================================================================
%  9. VALIDAÇÃO COM P_final
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  VALIDAÇÃO COM P_final (Rot + Trans separados)\n');
fprintf('==========================================================\n');

sg = segs{1};
[res_final] = simulate_full(P_final, sg.time, sg.pwm, sg.pqr, sg.acc, sg.att, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, ode_opts);

print_R2('P_final', res_final, sg.pqr, pqr_vl, sg.acc, acc_vl, ...
    t_trains{1}, t_val, R2_func);

% Resumo comparável entre rodadas (outputs/runs/<tag>/summary.mat)
summary = struct('tag', RUN_TAG, 'log', LOG_FILE, 't_val', t_val, ...
    'P_final', P_final, 'P0', P0, 'P_eem', P_eem, 'param_names', {param_names}, ...
    'AERO_FREE', AERO_FREE, 'V_p95', prctile(V_grid,95), 'V_max', max(V_grid));
summary.R2_val = [R2_func(pqr_vl(:,1), res_final.p_s_vl), ...
                  R2_func(pqr_vl(:,2), res_final.q_s_vl), ...
                  R2_func(pqr_vl(:,3), res_final.r_s_vl)];
summary.R2_acc = [R2_func(acc_vl(:,1), res_final.accX_s_vl), ...
                  R2_func(acc_vl(:,2), res_final.accY_s_vl), ...
                  R2_func(acc_vl(:,3), res_final.accZ_s_vl)];
try, summary.CR_pct = 100*se(:)'./max(abs(P_final(:)'),1e-9); catch, end
save(fullfile(out_dir,'summary.mat'), 'summary');

plot_all_results('P_final', ...
    res_final, sg.time, sg.pqr, sg.acc, time_vl, pqr_vl, acc_vl, R2_func, t_trains{1}, t_val, att_vl);

% Diagnóstico: Torques Mx, My, Mz na janela de validação
plot_torques(P_final, time_vl, pwm_vl, func_T_ref, func_Q_ref, pqr_vl, t_val);

% Diagnóstico: Forças Fx, Fy, Fz na janela de validação
plot_forces(P_final, time_vl, pwm_vl, att_vl, acc_vl, func_T_ref, constants_sim, t_val, res_final);

%% ========================================================================
%  9b. VALIDAÇÃO SEMI-ACOPLADA (consistente com a identificação)
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  VALIDAÇÃO SEMI-ACOPLADA (att medida + RK4 translacional)\n');
fprintf('==========================================================\n');

res_sc_final = simulate_full_semicoupled(P_final, sg.time, sg.pwm, sg.pqr, sg.acc, sg.att, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, []);

print_R2('P_final (semi-acoplado)', res_sc_final, sg.pqr, pqr_vl, sg.acc, acc_vl, ...
    t_trains{1}, t_val, R2_func);

plot_all_results('P_final (semi-acoplado)', ...
    res_sc_final, sg.time, sg.pqr, sg.acc, time_vl, pqr_vl, acc_vl, R2_func, t_trains{1}, t_val, att_vl);

res_sc_P0 = simulate_full_semicoupled(P0, sg.time, sg.pwm, sg.pqr, sg.acc, sg.att, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, []);

%% ========================================================================
%  9c. VALIDAÇÃO HÍBRIDA (att medida, pqr integrado HONESTO)
%      Isola o modelo ROTACIONAL: testa pqr sem o drift da att integrada.
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  VALIDAÇÃO HÍBRIDA (pqr integrado + att medida)\n');
fprintf('==========================================================\n');

res_hyb_final = simulate_full_hybrid(P_final, sg.time, sg.pwm, sg.pqr, sg.acc, sg.att, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, []);

print_R2('P_final (hibrido)', res_hyb_final, sg.pqr, pqr_vl, sg.acc, acc_vl, ...
    t_trains{1}, t_val, R2_func);

plot_all_results('P_final (hibrido)', ...
    res_hyb_final, sg.time, sg.pqr, sg.acc, time_vl, pqr_vl, acc_vl, R2_func, t_trains{1}, t_val, att_vl);

%% ========================================================================
%  9d. TESTE DE BRANCURA DO RESÍDUO (taxas, validação híbrida)
%  Resíduo BRANCO (ACF dentro das bandas 95%, exceto lag 0) ⇒ erro irredutível
%  → modelo no limite dos dados. COLORIDO (ACF estoura em lags>0) ⇒ dinâmica
%  não modelada (aero de rotor, lag de motor, asa...) → ainda melhorável.
%  Mesmo teste/figura do identify_v2, sobre a validação híbrida (rotacional puro).
%  ========================================================================
fprintf('\n--- Teste de brancura do resíduo (taxas, híbrido) ---\n');
maxlag_w = 30;
ewn = {'p', pqr_vl(:,1)-res_hyb_final.p_s_vl; ...
       'q', pqr_vl(:,2)-res_hyb_final.q_s_vl; ...
       'r', pqr_vl(:,3)-res_hyb_final.r_s_vl};
figW = figure('Name','identify_plant brancura','Position',[80 80 1100 700]);
for cW = 1:3
    eW = ewn{cW,2};
    [acfW, lagsW] = acf_white(eW, maxlag_w);
    NW = numel(eW);  bW = 1.96/sqrt(NW);
    insideW = mean(abs(acfW(2:end)) <= bW) * 100;
    if insideW >= 90, verdictW = '~BRANCO (irredutível)'; else, verdictW = 'COLORIDO (sobra dinâmica)'; end
    fprintf('  %-3s: %3.0f%% dos lags dentro de ±%.3f (95%%)  → %s\n', ewn{cW,1}, insideW, bW, verdictW);
    subplot(3,1,cW); hold on; grid on;
    stem(lagsW(2:end), acfW(2:end), 'filled', 'MarkerSize', 3);
    yline(bW,'r--'); yline(-bW,'r--');
    ylabel(['ACF ' ewn{cW,1}]); ylim([-0.5 0.5]);
    if cW==1, title(sprintf('Brancura do resíduo [%g-%gs] — dentro das bandas = branco', t_val(1), t_val(2))); end
    if cW==3, xlabel('lag (amostras)'); end
end
saveas(figW, fullfile(img_dir,'identify_plant_whiteness.png'));
fprintf('  (≥90%% dentro da banda ≈ branco. Figura: %s)\n', fullfile(img_dir,'identify_plant_whiteness.png'));

%% 10. RESUMO FINAL: P0 vs P_final (apenas validação)
fprintf('\n==========================================================\n');
fprintf('  RESUMO FINAL — Validação (%d-%ds)\n', t_val(1), t_val(2));
fprintf('==========================================================\n');
fprintf('  [P0] Chute Inicial:\n');
if res_P0.pqr_vl_ok
    fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
        R2_func(pqr_vl(:,1), res_P0.p_s_vl), R2_func(pqr_vl(:,2), res_P0.q_s_vl), R2_func(pqr_vl(:,3), res_P0.r_s_vl));
else, fprintf('    p,q,r -> DIVERGIU\n');
end
if res_P0.full_vl_ok
    fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
        R2_func(acc_vl(:,1), res_P0.accX_s_vl), R2_func(acc_vl(:,2), res_P0.accY_s_vl), R2_func(acc_vl(:,3), res_P0.accZ_s_vl));
end
fprintf('  [P_final] Após Otimização:\n');
if res_final.pqr_vl_ok
    fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
        R2_func(pqr_vl(:,1), res_final.p_s_vl), R2_func(pqr_vl(:,2), res_final.q_s_vl), R2_func(pqr_vl(:,3), res_final.r_s_vl));
else, fprintf('    p,q,r -> DIVERGIU\n');
end
if res_final.full_vl_ok
    fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
        R2_func(acc_vl(:,1), res_final.accX_s_vl), R2_func(acc_vl(:,2), res_final.accY_s_vl), R2_func(acc_vl(:,3), res_final.accZ_s_vl));
end
fprintf('\n  --- Semi-acoplado (att medida + RK4 translacional) ---\n');
fprintf('  [P0] Chute Inicial:\n');
fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
    R2_func(acc_vl(:,1), res_sc_P0.accX_s_vl), R2_func(acc_vl(:,2), res_sc_P0.accY_s_vl), R2_func(acc_vl(:,3), res_sc_P0.accZ_s_vl));
fprintf('  [P_final] Após Otimização:\n');
fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
    R2_func(acc_vl(:,1), res_sc_final.accX_s_vl), R2_func(acc_vl(:,2), res_sc_final.accY_s_vl), R2_func(acc_vl(:,3), res_sc_final.accZ_s_vl));
fprintf('==========================================================\n');

disp('Script finalizado.');

%% ========================================================================
%  FUNÇÕES LOCAIS
%  ========================================================================

function [aX, aY, aZ] = accel_eval(p, q, r, u, v, w, pwm, k_T, func_T_ref, m, r_imu, dt, kdm_P)
%ACCEL_EVAL  Força específica da IMU (modelo de sensor), unificando o bloco antes
%  repetido nas simulate_*: T=Σ k_T·fT(pwm); α = d/dt(p,q,r); accelerometer_model
%  com IMU offset. Saída: accX/accY/accZ.
    T = zeros(numel(p), 1);
    for j = 1:4, T = T + k_T(j) .* func_T_ref(pwm(:,j)); end
    if nargin < 13, kdm_P = []; end
    [aX, aY, aZ] = accelerometer_model(p, q, r, u, v, w, T./m, ...
        gradient(p, dt), gradient(q, dt), gradient(r, dt), r_imu, kdm_P);
end

function [acf, lags] = acf_white(e, maxlag)
% Autocorrelação normalizada (sem toolbox). acf(1)=lag 0=1.
    e = e(:) - mean(e);
    N = numel(e);  c0 = sum(e.^2) + 1e-12;
    acf = zeros(maxlag+1,1);
    for k = 0:maxlag
        acf(k+1) = sum(e(1:N-k) .* e(1+k:N)) / c0;
    end
    lags = (0:maxlag)';
end

function e = oem_multi_seg_cost(P, segs, win_sec, m, g, dt, weights_pqr, weights_acc, cost_mode, reg_lambda)
% Custo OEM multi-segmento: concatena resíduos de cada segmento.
% cost_mode: 'full' (padrão), 'rotational', 'translational'
% reg_lambda: struct .kt_pair .kq_pair (penaliza diferenças dentro de par CW/CCW)
    if nargin < 9, cost_mode = 'full'; end
    if nargin < 10, reg_lambda = struct('kt_pair', 0, 'kq_pair', 0); end
    if ~isfield(reg_lambda, 'kt_pair'), reg_lambda.kt_pair = 0; end
    if ~isfield(reg_lambda, 'kq_pair'), reg_lambda.kq_pair = 0; end

    e_all = [];
    for s = 1:length(segs)
        sg = segs{s};
        N = sg.N;

        if isinf(win_sec)
            win_len = N;
        else
            win_len = round(win_sec / dt);
        end

        win_starts = 1:win_len:N;
        if length(win_starts) > 1 && win_starts(end) == N
            win_starts(end) = [];
        end
        win_ends = [win_starts(2:end)-1, N];

        aero_s = [];
        if isfield(sg,'V'), aero_s = struct('V', sg.V(:), 'v', sg.v(:), 'w', sg.w(:)); end
        e_seg = oem_ms_cost_func(P, sg.pqr, sg.acc, sg.att_rad, ...
            sg.T_ref, sg.Q_ref, m, g, dt, N, ...
            win_starts, win_ends, weights_pqr, weights_acc, cost_mode, aero_s);

        e_all = [e_all; e_seg]; %#ok<AGROW>
    end

    % Regularização: penaliza diferenças dentro do mesmo par CW/CCW
    k_T = P(5:8);  k_Q = P(9:12);
    e_reg = [sqrt(reg_lambda.kt_pair) * (k_T(1) - k_T(2)); ...
             sqrt(reg_lambda.kt_pair) * (k_T(3) - k_T(4)); ...
             sqrt(reg_lambda.kq_pair) * (k_Q(1) - k_Q(2)); ...
             sqrt(reg_lambda.kq_pair) * (k_Q(3) - k_Q(4))];

    % Restrição física: desigualdade triangular dos momentos de inércia.
    %  ANTES só estava no EEM → o OEM ficava livre e violava Jz≤Jx+Jy. Agora
    %  o OEM também pune violação (one-sided: zero quando físico).
    tri_lam = 0; if isfield(reg_lambda,'tri'), tri_lam = reg_lambda.tri; end
    Jx=P(1); Jy=P(2); Jz=P(3);
    e_tri = sqrt(tri_lam) * [max(0, Jz-(Jx+Jy)); ...
                            max(0, Jy-(Jx+Jz)); ...
                            max(0, Jx-(Jy+Jz))];

    e = [e_all; e_reg; e_tri];
end

function e = oem_ms_cost_func(P, pqr, acc, att_rad, T_ref, Q_ref, m, g, dt, N, ...
    win_starts, win_ends, weights_pqr, weights_acc, cost_mode, aero)

    if nargin < 15, cost_mode = 'full'; end
    if nargin < 16, aero = []; end

    % Obter handles centralizados do vtol_dynamics (edite apenas lá!)
    dyn_h = vtol_dynamics('get_handles');
    trans_dot_fn = dyn_h.trans_dot;
    % Acelerômetro: arquivo separado (modelo de sensor, não de planta).
    % IMU offset r_imu lido de parameters().

    % Inércias → constantes G (corpo rígido, consistência garantida)
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
    gamma0 = Jx*Jz - Jxz^2;
    G1 = Jxz*(Jx - Jy + Jz) / gamma0;
    G2 = (Jz*(Jz - Jy) + Jxz^2) / gamma0;
    G3 = Jz / gamma0;
    G4 = Jxz / gamma0;
    G5 = (Jz - Jx) / Jy;
    G6 = Jxz / Jy;
    G7 = (Jx*(Jx - Jy) + Jxz^2) / gamma0;
    G8 = Jx / gamma0;
    invJy = 1 / Jy;

    k_T = P(5:8); k_Q = P(9:12);
    Dp = P(13); Dq = P(14); Dr = P(15);

    % Braços + IMU offset (de parameters() — fonte única)
    proj_p_oem = parameters();
    % FORMA DO AMORTECIMENTO (parameters().damp_form). Em 'moment' os P(13:15)
    % sao L_p, M_q, N_r [N.m.s] e entram DENTRO do vetor de momentos, junto com
    % a aerodinamica; em 'rate' sao c_p, c_q, c_r [1/s] subtraidos da derivada.
    Lpm = 0; Mqm = 0; Nrm = 0;  Dpr = Dp; Dqr = Dq; Drr = Dr;
    if isfield(proj_p_oem,'damp_form') && strcmp(proj_p_oem.damp_form,'moment')
        Lpm = Dp; Mqm = Dq; Nrm = Dr;  Dpr = 0; Dqr = 0; Drr = 0;
    end
    Lx_r = proj_p_oem.arms.Lx_r;
    Lx_l = proj_p_oem.arms.Lx_l;
    Ly_f = proj_p_oem.arms.Ly_f;
    Ly_r = proj_p_oem.arms.Ly_r;
    r_imu = proj_p_oem.imu_offset;

    % --- AERODINÂMICA DA ESTRUTURA (P(17:22), ∝ V) — coeficientes por amostra ---
    % aLp(k) = ¼ρSb²·Cl_p·V(k) multiplica o p SIMULADO (é amortecimento, entra
    % dentro do RK4); aLb(k) = ½ρSb·Cl_β·V(k)·v(k) é exógeno (parcela de
    % derrapagem). Idem para arfagem e guinada. V, v, w vêm de estimate_velocity.
    use_aero = ~isempty(aero) && numel(P) >= 22;
    if use_aero
        kA = aero_gains(proj_p_oem);
        Va = aero.V(:);  vb = aero.v(:);  wb = aero.w(:);
        Va(~isfinite(Va)) = 0;  vb(~isfinite(vb)) = 0;  wb(~isfinite(wb)) = 0;
        aLp = kA.Lp*P(17)*Va;   aLb = kA.Lb*P(18)*Va.*vb;
        aMq = kA.Mq*P(19)*Va;   aMa = kA.Ma*P(20)*Va.*wb;
        aNr = kA.Nr*P(21)*Va;   aNb = kA.Nb*P(22)*Va.*vb;
        if numel(P) >= 25
            if isappdata(0,'hu_form')   % forma do Hu: −½ρV²S·(b,c̄,b)·(c_l,c_m,c_n)
                qS = 0.5*kA.rho*Va.^2*kA.S;
                aLb = aLb - qS*kA.b*P(23);  aMa = aMa - qS*kA.c*P(24);  aNb = aNb - qS*kA.b*P(25);
            else                        % rotor: −L_v·v, −M_w·w, −N_v·v
                aLb = aLb - P(23)*vb;  aMa = aMa - P(24)*wb;  aNb = aNb - P(25)*vb;
            end
        end
    else
        aLp = zeros(N,1); aLb = aLp; aMq = aLp; aMa = aLp; aNr = aLp; aNb = aLp;
    end

    n_sub = 5;
    dt_sub = dt / n_sub;
    h2 = dt_sub / 2;
    MAX_VAL = 50;

    p_sim = zeros(N, 1);
    q_sim = zeros(N, 1);
    r_sim = zeros(N, 1);

    for w = 1:length(win_starts)
        i_s = win_starts(w);
        i_e = win_ends(w);
        ps = pqr(i_s,1); qs = pqr(i_s,2); rs = pqr(i_s,3);
        p_sim(i_s) = ps; q_sim(i_s) = qs; r_sim(i_s) = rs;

        for k = i_s:i_e-1
            % Momentos no início (k) e no fim (k+1) do passo — ArduPilot QuadX
            % (M1,M3=FRONT; M2,M4=REAR). O empuxo vem de motor_chain (RPM contínua
            % interpolada linearmente), então M(t) é interpolado LINEARMENTE dentro
            % do passo, igual ao que ode45/vtol_dynamics fazem na validação.
            % (Antes: M constante por 0,1 s → ≈50 ms de atraso efetivo a mais.)
            Tmr0 = k_T(:)' .* T_ref(k,:);    Qmr0 = k_Q(:)' .* Q_ref(k,:);
            Tmr1 = k_T(:)' .* T_ref(k+1,:);  Qmr1 = k_Q(:)' .* Q_ref(k+1,:);
            Mx0 = -(Lx_r*Tmr0(1) - Lx_l*Tmr0(2) - Lx_l*Tmr0(3) + Lx_r*Tmr0(4));
            My0 =  Ly_f*Tmr0(1) - Ly_r*Tmr0(2) + Ly_f*Tmr0(3) - Ly_r*Tmr0(4);
            Mz0 = Qmr0(1) + Qmr0(2) - Qmr0(3) - Qmr0(4);
            Mx1 = -(Lx_r*Tmr1(1) - Lx_l*Tmr1(2) - Lx_l*Tmr1(3) + Lx_r*Tmr1(4));
            My1 =  Ly_f*Tmr1(1) - Ly_r*Tmr1(2) + Ly_f*Tmr1(3) - Ly_r*Tmr1(4);
            Mz1 = Qmr1(1) + Qmr1(2) - Qmr1(3) - Qmr1(4);
            % parcelas aerodinâmicas EXÓGENAS (derrapagem/incidência) somam-se aos
            % momentos; as parcelas de AMORTECIMENTO (∝ V·p, V·q, V·r) dependem do
            % estado e entram estágio a estágio, logo abaixo.
            Mx0 = Mx0 + aLb(k);  Mx1 = Mx1 + aLb(k+1);
            My0 = My0 + aMa(k);  My1 = My1 + aMa(k+1);
            Mz0 = Mz0 + aNb(k);  Mz1 = Mz1 + aNb(k+1);
            dMx = Mx1-Mx0; dMy = My1-My0; dMz = Mz1-Mz0;
            aLp0 = aLp(k); daLp = aLp(k+1)-aLp(k);
            aMq0 = aMq(k); daMq = aMq(k+1)-aMq(k);
            aNr0 = aNr(k); daNr = aNr(k+1)-aNr(k);

            for si = 1:n_sub
                % RK4 sub-stepping com M(t) linear no passo: frações dos estágios
                fa = (si-1)/n_sub;  fb = fa + 0.5/n_sub;  fc = fa + 1/n_sub;
                Mxa = Mx0+fa*dMx; Mya = My0+fa*dMy; Mza = Mz0+fa*dMz;
                Mxb = Mx0+fb*dMx; Myb = My0+fb*dMy; Mzb = Mz0+fb*dMz;
                Mxc = Mx0+fc*dMx; Myc = My0+fc*dMy; Mzc = Mz0+fc*dMz;
                % coeficientes que multiplicam o estado DENTRO de M: aerodinamica
                % (aLp etc., ja negativa) menos o amortecimento em forma de momento
                La = aLp0+fa*daLp - Lpm; Lb_ = aLp0+fb*daLp - Lpm; Lc = aLp0+fc*daLp - Lpm;
                Qa = aMq0+fa*daMq - Mqm; Qb = aMq0+fb*daMq - Mqm; Qc = aMq0+fc*daMq - Mqm;
                Ra = aNr0+fa*daNr - Nrm; Rb = aNr0+fb*daNr - Nrm; Rc = aNr0+fc*daNr - Nrm;
                % k1
                MxA = Mxa + La*ps;  MyA = Mya + Qa*qs;  MzA = Mza + Ra*rs;
                pd1 = G1*ps*qs - G2*qs*rs + G3*MxA + G4*MzA - Dpr*ps;
                qd1 = G5*ps*rs - G6*(ps^2 - rs^2) + invJy*MyA - Dqr*qs;
                rd1 = G7*ps*qs - G1*qs*rs + G4*MxA + G8*MzA - Drr*rs;
                % k2
                p2 = ps+h2*pd1; q2 = qs+h2*qd1; r2 = rs+h2*rd1;
                MxB = Mxb + Lb_*p2;  MyB = Myb + Qb*q2;  MzB = Mzb + Rb*r2;
                pd2 = G1*p2*q2 - G2*q2*r2 + G3*MxB + G4*MzB - Dpr*p2;
                qd2 = G5*p2*r2 - G6*(p2^2 - r2^2) + invJy*MyB - Dqr*q2;
                rd2 = G7*p2*q2 - G1*q2*r2 + G4*MxB + G8*MzB - Drr*r2;
                % k3
                p3 = ps+h2*pd2; q3 = qs+h2*qd2; r3 = rs+h2*rd2;
                MxC = Mxb + Lb_*p3;  MyC = Myb + Qb*q3;  MzC = Mzb + Rb*r3;
                pd3 = G1*p3*q3 - G2*q3*r3 + G3*MxC + G4*MzC - Dpr*p3;
                qd3 = G5*p3*r3 - G6*(p3^2 - r3^2) + invJy*MyC - Dqr*q3;
                rd3 = G7*p3*q3 - G1*q3*r3 + G4*MxC + G8*MzC - Drr*r3;
                % k4
                p4 = ps+dt_sub*pd3; q4 = qs+dt_sub*qd3; r4 = rs+dt_sub*rd3;
                MxD = Mxc + Lc*p4;  MyD = Myc + Qc*q4;  MzD = Mzc + Rc*r4;
                pd4 = G1*p4*q4 - G2*q4*r4 + G3*MxD + G4*MzD - Dpr*p4;
                qd4 = G5*p4*r4 - G6*(p4^2 - r4^2) + invJy*MyD - Dqr*q4;
                rd4 = G7*p4*q4 - G1*q4*r4 + G4*MxD + G8*MzD - Drr*r4;
                % update
                ps = ps + dt_sub/6*(pd1 + 2*pd2 + 2*pd3 + pd4);
                qs = qs + dt_sub/6*(qd1 + 2*qd2 + 2*qd3 + qd4);
                rs = rs + dt_sub/6*(rd1 + 2*rd2 + 2*rd3 + rd4);
            end

            if ~isfinite(ps), ps = MAX_VAL; end
            if ~isfinite(qs), qs = MAX_VAL; end
            if ~isfinite(rs), rs = MAX_VAL; end
            ps = max(min(ps, MAX_VAL), -MAX_VAL);
            qs = max(min(qs, MAX_VAL), -MAX_VAL);
            rs = max(min(rs, MAX_VAL), -MAX_VAL);

            p_sim(k+1) = ps; q_sim(k+1) = qs; r_sim(k+1) = rs;
        end
    end

    sw_r = sqrt(weights_pqr(:));
    e_rot = [sw_r(1)*(pqr(:,1) - p_sim); ...
             sw_r(2)*(pqr(:,2) - q_sim); ...
             sw_r(3)*(pqr(:,3) - r_sim)];

    T_total = k_T(1)*T_ref(:,1) + k_T(2)*T_ref(:,2) + ...
              k_T(3)*T_ref(:,3) + k_T(4)*T_ref(:,4);
    % influxo no empuxo total: T ← T − 4·k_v·w (k_v medido, w da velocidade estimada)
    kvT = proj_p_oem.k_v_thrust;
    if kvT > 0 && use_aero, T_total = T_total - 4*kvT*wb; end

    phi_r = att_rad(:,1); theta_r = att_rad(:,2);
    gx = -g * sin(theta_r);
    gy = g * cos(theta_r) .* sin(phi_r);
    gz = g * cos(theta_r) .* cos(phi_r);

    % Integração translacional RK4 com sub-stepping (consistente c/ rotacional)
    n_sub_t = n_sub;        % mesmo número de sub-passos
    dt_sub_t = dt / n_sub_t;
    h2_t = dt_sub_t / 2;

    % Arrasto translacional por amostra, Nx3 [1/s]:
    %   x, y : arrasto induzido de rotor (g·C_d, constante) + estrutura ∝ V
    %   z    : só estrutura (o C_d de Beard age em x,y)
    kdm_P = [];
    if numel(P) >= 16, kdm_P = repmat(g*P(16)*[1 1 0], N, 1); end
    if numel(P) >= 25 && any(P(23:25) ~= 0) && use_aero
        kF = aero_gains(proj_p_oem);
        kdm_P = kdm_P + (kF.F/m) * (Va(:) * P(23:25)');
    end

    u_int = zeros(N,1); v_int = zeros(N,1); w_int = zeros(N,1);
    for k = 1:N-1
        pk = pqr(k,1); qk = pqr(k,2); rk = pqr(k,3);
        gxk = gx(k); gyk = gy(k); gzk = gz(k);
        Tk_m = T_total(k)/m;

        us = u_int(k); vs = v_int(k); ws = w_int(k);
        kdm_k = [];  if ~isempty(kdm_P), kdm_k = kdm_P(k,:); end
        for si = 1:n_sub_t
            % k1
            [ud1,vd1,wd1] = trans_dot_fn(pk,qk,rk, us,vs,ws, gxk,gyk,gzk, Tk_m, kdm_k);
            % k2
            u2=us+h2_t*ud1; v2=vs+h2_t*vd1; w2=ws+h2_t*wd1;
            [ud2,vd2,wd2] = trans_dot_fn(pk,qk,rk, u2,v2,w2, gxk,gyk,gzk, Tk_m, kdm_k);
            % k3
            u3=us+h2_t*ud2; v3=vs+h2_t*vd2; w3=ws+h2_t*wd2;
            [ud3,vd3,wd3] = trans_dot_fn(pk,qk,rk, u3,v3,w3, gxk,gyk,gzk, Tk_m, kdm_k);
            % k4
            u4=us+dt_sub_t*ud3; v4=vs+dt_sub_t*vd3; w4=ws+dt_sub_t*wd3;
            [ud4,vd4,wd4] = trans_dot_fn(pk,qk,rk, u4,v4,w4, gxk,gyk,gzk, Tk_m, kdm_k);
            % update
            us = us + dt_sub_t/6*(ud1 + 2*ud2 + 2*ud3 + ud4);
            vs = vs + dt_sub_t/6*(vd1 + 2*vd2 + 2*vd3 + vd4);
            ws = ws + dt_sub_t/6*(wd1 + 2*wd2 + 2*wd3 + wd4);
        end

        if ~isfinite(us), us = MAX_VAL; end
        if ~isfinite(vs), vs = MAX_VAL; end
        if ~isfinite(ws), ws = MAX_VAL; end
        us = max(min(us, MAX_VAL), -MAX_VAL);
        vs = max(min(vs, MAX_VAL), -MAX_VAL);
        ws = max(min(ws, MAX_VAL), -MAX_VAL);

        u_int(k+1) = us; v_int(k+1) = vs; w_int(k+1) = ws;
    end

    % Saída do modelo (acelerômetro como sensor — arquivo separado).
    % α via derivada numérica de p, q, r (no caso da cost: pqr medido).
    p_dot_sig = gradient(pqr(:,1), dt);
    q_dot_sig = gradient(pqr(:,2), dt);
    r_dot_sig = gradient(pqr(:,3), dt);
    [accX_m, accY_m, accZ_m] = accelerometer_model( ...
        pqr(:,1), pqr(:,2), pqr(:,3), ...
        u_int, v_int, w_int, ...
        T_total/m, ...
        p_dot_sig, q_dot_sig, r_dot_sig, ...
        r_imu, kdm_P);

    sw_a = sqrt(weights_acc(:));
    e_acc = [sw_a(1)*(acc(:,1) - accX_m); ...
             sw_a(2)*(acc(:,2) - accY_m); ...
             sw_a(3)*(acc(:,3) - accZ_m)];

    if strcmp(cost_mode, 'rotational')
        e = e_rot;
    elseif strcmp(cost_mode, 'translational')
        e = e_acc;
    else
        e = [e_rot; e_acc];
    end
end

function res = simulate_full(P, time_tr, pwm_tr, pqr_tr, acc_tr, att_tr, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, ode_opts)
% Simulação completa model-driven usando integração unificada de 9 estados.
% Replica exatamente o Simulink: todos os subsistemas (rotacional, cinemática,
% translacional) são integrados simultaneamente no mesmo solver.
% Nenhum dado medido alimenta a simulação — totalmente model-driven.

    N_tr = length(time_tr);
    N_vl = length(time_vl);
    g = constants_sim.g;

    res = simulate_one_window(P, time_tr, pwm_tr, pqr_tr, att_tr, N_tr, ...
        func_T_ref, func_Q_ref, constants_sim, ode_opts, g, 'tr');
    res2 = simulate_one_window(P, time_vl, pwm_vl, pqr_vl, att_vl, N_vl, ...
        func_T_ref, func_Q_ref, constants_sim, ode_opts, g, 'vl');

    % Mesclar resultados validação no struct
    fnames = fieldnames(res2);
    for i = 1:length(fnames)
        res.(fnames{i}) = res2.(fnames{i});
    end
end

function res = simulate_one_window(P, time, pwm, pqr, att, N, ...
    func_T_ref, func_Q_ref, constants_sim, ode_opts, g, suffix)
% Integração unificada de 9 estados via vtol_dynamics.
%
%   y = [p; q; r; phi; theta; psi; u; v; w]
%
% Todos os subsistemas (rotacional, cinemática de Euler, translacional) são
% integrados simultaneamente pelo mesmo solver — idêntico ao Simulink.
% Nenhuma interpolação entre etapas, nenhuma reamostragem.

    nan3 = NaN(N,1);
    kdm_P = kdm_of(P, time, constants_sim.m);   % Nx3: rotor + estrutura ∝ V

    % Condição inicial: [p,q,r] e [phi,theta,psi] dos dados medidos, [u,v,w]=0
    att_rad0 = deg2rad(att(1,:));
    y0 = [pqr(1,:)'; att_rad0(:); 0; 0; 0];

    % Integração unificada de 9 estados
    ode_func = @(t,y) vtol_dynamics(t, y, P, time, pwm, func_T_ref, func_Q_ref, constants_sim);
    try
        [t_s, y_s] = ode45(ode_func, time, y0, ode_opts);
        y_out = interp1(t_s, y_s, time, 'linear', 'extrap');

        % Extrair estados
        res.(['p_s_' suffix]) = y_out(:,1);
        res.(['q_s_' suffix]) = y_out(:,2);
        res.(['r_s_' suffix]) = y_out(:,3);
        res.(['pqr_' suffix '_ok']) = true;

        res.(['phi_s_' suffix])   = rad2deg(y_out(:,4));
        res.(['theta_s_' suffix]) = rad2deg(y_out(:,5));
        res.(['psi_s_' suffix])   = rad2deg(y_out(:,6));

        % Acelerômetro (sensor): pqr e u,v,w SIMULADOS (9 estados).
        [accX, accY, accZ] = accel_eval(y_out(:,1), y_out(:,2), y_out(:,3), ...
            y_out(:,7), y_out(:,8), y_out(:,9), pwm, P(5:8), func_T_ref, ...
            constants_sim.m, parameters().imu_offset, time(2)-time(1), kdm_P);
        res.(['accX_s_' suffix]) = accX;
        res.(['accY_s_' suffix]) = accY;
        res.(['accZ_s_' suffix]) = accZ;
        res.(['full_' suffix '_ok']) = true;

    catch ME
        fprintf('  [%s] ode45 9-estados falhou: %s\n', suffix, ME.message);
        res.(['p_s_' suffix]) = nan3; res.(['q_s_' suffix]) = nan3; res.(['r_s_' suffix]) = nan3;
        res.(['phi_s_' suffix]) = nan3; res.(['theta_s_' suffix]) = nan3; res.(['psi_s_' suffix]) = nan3;
        res.(['accX_s_' suffix]) = nan3; res.(['accY_s_' suffix]) = nan3; res.(['accZ_s_' suffix]) = nan3;
        res.(['pqr_' suffix '_ok']) = false;
        res.(['full_' suffix '_ok']) = false;
    end
end

function res = simulate_full_semicoupled(P, time_tr, pwm_tr, pqr_tr, acc_tr, att_tr, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, ~)
% Validação semi-acoplada: usa p,q,r e atitude MEDIDOS, integra apenas u,v,w.
% Consistente com oem_ms_cost_func (mesma lógica que a identificação usa).
% Retorna struct com mesmos campos que simulate_full (compatível com plots).

    dyn_h = vtol_dynamics('get_handles');
    trans_dot_fn = dyn_h.trans_dot;
    % Acelerômetro: arquivo separado (sensor)

    g = constants_sim.g;
    m = constants_sim.m;
    kdm_P = [];  if numel(P) >= 16, kdm_P = g * P(16) * [1 1 0]; end   % C_d → g·C_d (x,y)
    kdm_k = kdm_P;   % simuladores auxiliares: sem as forças ∝ V da estrutura
    k_T = P(5:8);

    proj_p_sc = parameters();
    r_imu_sc = proj_p_sc.imu_offset;

    res = struct();
    datasets = {{'tr', time_tr, pwm_tr, pqr_tr, acc_tr, att_tr}, ...
                {'vl', time_vl, pwm_vl, pqr_vl, acc_vl, att_vl}};

    for d = 1:2
        suffix = datasets{d}{1};
        time   = datasets{d}{2};
        pwm    = datasets{d}{3};
        pqr    = datasets{d}{4};
        att    = datasets{d}{6};
        N      = length(time);
        dt     = time(2) - time(1);

        % p,q,r e atitude = dados medidos (sem modelo rotacional)
        res.(['p_s_' suffix]) = pqr(:,1);
        res.(['q_s_' suffix]) = pqr(:,2);
        res.(['r_s_' suffix]) = pqr(:,3);
        res.(['pqr_' suffix '_ok']) = true;

        res.(['phi_s_' suffix])   = att(:,1);
        res.(['theta_s_' suffix]) = att(:,2);
        res.(['psi_s_' suffix])   = att(:,3);

        % Empuxo total
        T_total = zeros(N,1);
        for k = 1:N
            for j = 1:4
                T_total(k) = T_total(k) + k_T(j) * func_T_ref(pwm(k,j));
            end
        end

        % Gravidade projetada a partir de atitude MEDIDA
        att_rad = deg2rad(att);
        phi_m = att_rad(:,1); theta_m = att_rad(:,2);
        gx = -g * sin(theta_m);
        gy =  g * cos(theta_m) .* sin(phi_m);
        gz =  g * cos(theta_m) .* cos(phi_m);

        % RK4 translacional com sub-stepping (u,v,w apenas)
        n_sub = 5;
        dt_sub = dt / n_sub;
        h2 = dt_sub / 2;

        u_int = zeros(N,1); v_int = zeros(N,1); w_int = zeros(N,1);
        for k = 1:N-1
            pk = pqr(k,1); qk = pqr(k,2); rk = pqr(k,3);
            gxk = gx(k); gyk = gy(k); gzk = gz(k);
            Tk_m = T_total(k)/m;

            us = u_int(k); vs = v_int(k); ws = w_int(k);
            for si = 1:n_sub
                [ud1,vd1,wd1] = trans_dot_fn(pk,qk,rk, us,vs,ws, gxk,gyk,gzk, Tk_m, kdm_k);
                u2=us+h2*ud1; v2=vs+h2*vd1; w2=ws+h2*wd1;
                [ud2,vd2,wd2] = trans_dot_fn(pk,qk,rk, u2,v2,w2, gxk,gyk,gzk, Tk_m, kdm_k);
                u3=us+h2*ud2; v3=vs+h2*vd2; w3=ws+h2*wd2;
                [ud3,vd3,wd3] = trans_dot_fn(pk,qk,rk, u3,v3,w3, gxk,gyk,gzk, Tk_m, kdm_k);
                u4=us+dt_sub*ud3; v4=vs+dt_sub*vd3; w4=ws+dt_sub*wd3;
                [ud4,vd4,wd4] = trans_dot_fn(pk,qk,rk, u4,v4,w4, gxk,gyk,gzk, Tk_m, kdm_k);
                us = us + dt_sub/6*(ud1+2*ud2+2*ud3+ud4);
                vs = vs + dt_sub/6*(vd1+2*vd2+2*vd3+vd4);
                ws = ws + dt_sub/6*(wd1+2*wd2+2*wd3+wd4);
            end
            u_int(k+1) = us; v_int(k+1) = vs; w_int(k+1) = ws;
        end

        % Acelerômetro (sensor): pqr MEDIDO + u,v,w integrados.
        [accX, accY, accZ] = accel_eval(pqr(:,1), pqr(:,2), pqr(:,3), ...
            u_int, v_int, w_int, pwm, k_T, func_T_ref, m, r_imu_sc, dt, kdm_P);

        res.(['accX_s_' suffix]) = accX;
        res.(['accY_s_' suffix]) = accY;
        res.(['accZ_s_' suffix]) = accZ;
        res.(['full_' suffix '_ok']) = true;
    end
end

function res = simulate_full_hybrid(P, time_tr, pwm_tr, pqr_tr, acc_tr, att_tr, ...
    time_vl, pwm_vl, pqr_vl, acc_vl, att_vl, ...
    func_T_ref, func_Q_ref, constants_sim, ~)
% Validação HÍBRIDA: integra p,q,r normalmente (testa rotacional) mas
% usa atitude MEDIDA do EKF para gravidade no translacional (sem drift).
% Diferente do semi-acoplado: pqr aqui é INTEGRADO, não passthrough.

    dyn_h = vtol_dynamics('get_handles');
    trans_dot_fn = dyn_h.trans_dot;
    % Acelerômetro: arquivo separado (sensor)

    g = constants_sim.g;
    m = constants_sim.m;
    kdm_P = [];  if numel(P) >= 16, kdm_P = g * P(16) * [1 1 0]; end   % C_d → g·C_d (x,y)
    kdm_k = kdm_P;   % simuladores auxiliares: sem as forças ∝ V da estrutura

    % Parâmetros (igual ao oem_ms_cost_func)
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
    gamma0 = Jx*Jz - Jxz^2;
    G1 = Jxz*(Jx - Jy + Jz) / gamma0;
    G2 = (Jz*(Jz - Jy) + Jxz^2) / gamma0;
    G3 = Jz / gamma0;
    G4 = Jxz / gamma0;
    G5 = (Jz - Jx) / Jy;
    G6 = Jxz / Jy;
    G7 = (Jx*(Jx - Jy) + Jxz^2) / gamma0;
    G8 = Jx / gamma0;
    invJy = 1 / Jy;

    k_T = P(5:8); k_Q = P(9:12);
    Dp = P(13); Dq = P(14); Dr = P(15);

    % Braços + IMU offset (de parameters() — fonte única)
    proj_p_sim = parameters();
    % FORMA DO AMORTECIMENTO (parameters().damp_form). Em 'moment' os P(13:15)
    % sao L_p, M_q, N_r [N.m.s] e entram DENTRO do vetor de momentos, junto com
    % a aerodinamica; em 'rate' sao c_p, c_q, c_r [1/s] subtraidos da derivada.
    Lpm = 0; Mqm = 0; Nrm = 0;  Dpr = Dp; Dqr = Dq; Drr = Dr;
    if isfield(proj_p_sim,'damp_form') && strcmp(proj_p_sim.damp_form,'moment')
        Lpm = Dp; Mqm = Dq; Nrm = Dr;  Dpr = 0; Dqr = 0; Drr = 0;
    end
    Lx_r = proj_p_sim.arms.Lx_r; Lx_l = proj_p_sim.arms.Lx_l;
    Ly_f = proj_p_sim.arms.Ly_f; Ly_r = proj_p_sim.arms.Ly_r;
    r_imu_hyb = proj_p_sim.imu_offset;

    res = struct();
    datasets = {{'tr', time_tr, pwm_tr, pqr_tr, acc_tr, att_tr}, ...
                {'vl', time_vl, pwm_vl, pqr_vl, acc_vl, att_vl}};

    for d = 1:2
        suffix = datasets{d}{1};
        time   = datasets{d}{2};
        pwm    = datasets{d}{3};
        pqr    = datasets{d}{4};
        att    = datasets{d}{6};
        N = length(time);
        dt = time(2) - time(1);

        % T_ref e Q_ref
        T_ref = zeros(N, 4);
        Q_ref = zeros(N, 4);
        for j = 1:4
            T_ref(:,j) = func_T_ref(pwm(:,j));
            Q_ref(:,j) = func_Q_ref(pwm(:,j));
        end

        % ============================================================
        % INTEGRAR p, q, r (RK4 sub-steps) - igual ao oem_ms_cost_func
        % ============================================================
        n_sub = 5;
        dt_sub = dt / n_sub;
        h2 = dt_sub / 2;
        MAX_VAL = 50;

        p_sim = zeros(N, 1); q_sim = zeros(N, 1); r_sim = zeros(N, 1);
        ps = pqr(1,1); qs = pqr(1,2); rs = pqr(1,3);   % IC do primeiro medido
        p_sim(1) = ps; q_sim(1) = qs; r_sim(1) = rs;
        [aLp, aLb, aMq, aMa, aNr, aNb] = aero_terms_local(P, time);

        for k = 1:N-1
            Tmr = k_T(:)' .* T_ref(k,:);
            Qmr = k_Q(:)' .* Q_ref(k,:);
            % Momentos — ArduPilot QuadX (M1,M3=FRONT; M2,M4=REAR)
            Mx = -(Lx_r*Tmr(1) - Lx_l*Tmr(2) - Lx_l*Tmr(3) + Lx_r*Tmr(4));
            My =  Ly_f*Tmr(1) - Ly_r*Tmr(2) + Ly_f*Tmr(3) - Ly_r*Tmr(4);
            Mz = Qmr(1) + Qmr(2) - Qmr(3) - Qmr(4);
            % aerodinâmica da estrutura (constante no passo, como os momentos)
            Mx = Mx + aLb(k);  My = My + aMa(k);  Mz = Mz + aNb(k);
            cLp = aLp(k) - Lpm; cMq = aMq(k) - Mqm; cNr = aNr(k) - Nrm;

            for si = 1:n_sub
                MxA = Mx + cLp*ps;  MyA = My + cMq*qs;  MzA = Mz + cNr*rs;
                pd1 = G1*ps*qs - G2*qs*rs + G3*MxA + G4*MzA - Dpr*ps;
                qd1 = G5*ps*rs - G6*(ps^2 - rs^2) + invJy*MyA - Dqr*qs;
                rd1 = G7*ps*qs - G1*qs*rs + G4*MxA + G8*MzA - Drr*rs;
                p2 = ps+h2*pd1; q2 = qs+h2*qd1; r2 = rs+h2*rd1;
                MxB = Mx + cLp*p2;  MyB = My + cMq*q2;  MzB = Mz + cNr*r2;
                pd2 = G1*p2*q2 - G2*q2*r2 + G3*MxB + G4*MzB - Dpr*p2;
                qd2 = G5*p2*r2 - G6*(p2^2 - r2^2) + invJy*MyB - Dqr*q2;
                rd2 = G7*p2*q2 - G1*q2*r2 + G4*MxB + G8*MzB - Drr*r2;
                p3 = ps+h2*pd2; q3 = qs+h2*qd2; r3 = rs+h2*rd2;
                MxC = Mx + cLp*p3;  MyC = My + cMq*q3;  MzC = Mz + cNr*r3;
                pd3 = G1*p3*q3 - G2*q3*r3 + G3*MxC + G4*MzC - Dpr*p3;
                qd3 = G5*p3*r3 - G6*(p3^2 - r3^2) + invJy*MyC - Dqr*q3;
                rd3 = G7*p3*q3 - G1*q3*r3 + G4*MxC + G8*MzC - Drr*r3;
                p4 = ps+dt_sub*pd3; q4 = qs+dt_sub*qd3; r4 = rs+dt_sub*rd3;
                MxD = Mx + cLp*p4;  MyD = My + cMq*q4;  MzD = Mz + cNr*r4;
                pd4 = G1*p4*q4 - G2*q4*r4 + G3*MxD + G4*MzD - Dpr*p4;
                qd4 = G5*p4*r4 - G6*(p4^2 - r4^2) + invJy*MyD - Dqr*q4;
                rd4 = G7*p4*q4 - G1*q4*r4 + G4*MxD + G8*MzD - Drr*r4;
                ps = ps + dt_sub/6 * (pd1 + 2*pd2 + 2*pd3 + pd4);
                qs = qs + dt_sub/6 * (qd1 + 2*qd2 + 2*qd3 + qd4);
                rs = rs + dt_sub/6 * (rd1 + 2*rd2 + 2*rd3 + rd4);
            end
            if ~isfinite(ps), ps = MAX_VAL; end
            if ~isfinite(qs), qs = MAX_VAL; end
            if ~isfinite(rs), rs = MAX_VAL; end
            ps = max(min(ps, MAX_VAL), -MAX_VAL);
            qs = max(min(qs, MAX_VAL), -MAX_VAL);
            rs = max(min(rs, MAX_VAL), -MAX_VAL);
            p_sim(k+1) = ps; q_sim(k+1) = qs; r_sim(k+1) = rs;
        end

        res.(['p_s_' suffix]) = p_sim;
        res.(['q_s_' suffix]) = q_sim;
        res.(['r_s_' suffix]) = r_sim;
        res.(['pqr_' suffix '_ok']) = true;

        % Atitude: pra plot, mostra a MEDIDA (já que não é integrada)
        res.(['phi_s_' suffix])   = att(:,1);
        res.(['theta_s_' suffix]) = att(:,2);
        res.(['psi_s_' suffix])   = att(:,3);

        % ============================================================
        % INTEGRAR u, v, w usando pqr SIMULADO + att MEDIDA pra gravidade
        % ============================================================
        T_total = zeros(N,1);
        for k = 1:N
            T_total(k) = sum(k_T(:)' .* T_ref(k,:));
        end

        att_rad = deg2rad(att);
        phi_m = att_rad(:,1); theta_m = att_rad(:,2);
        gx = -g * sin(theta_m);
        gy =  g * cos(theta_m) .* sin(phi_m);
        gz =  g * cos(theta_m) .* cos(phi_m);

        n_sub_t = n_sub;
        dt_sub_t = dt / n_sub_t;
        h2_t = dt_sub_t / 2;

        u_int = zeros(N,1); v_int = zeros(N,1); w_int = zeros(N,1);
        for k = 1:N-1
            pk = p_sim(k); qk = q_sim(k); rk = r_sim(k);  % ← simulado, não medido
            gxk = gx(k); gyk = gy(k); gzk = gz(k);
            Tk_m = T_total(k)/m;

            us = u_int(k); vs = v_int(k); ws = w_int(k);
            for si = 1:n_sub_t
                [ud1,vd1,wd1] = trans_dot_fn(pk,qk,rk, us,vs,ws, gxk,gyk,gzk, Tk_m, kdm_k);
                u2=us+h2_t*ud1; v2=vs+h2_t*vd1; w2=ws+h2_t*wd1;
                [ud2,vd2,wd2] = trans_dot_fn(pk,qk,rk, u2,v2,w2, gxk,gyk,gzk, Tk_m, kdm_k);
                u3=us+h2_t*ud2; v3=vs+h2_t*vd2; w3=ws+h2_t*wd2;
                [ud3,vd3,wd3] = trans_dot_fn(pk,qk,rk, u3,v3,w3, gxk,gyk,gzk, Tk_m, kdm_k);
                u4=us+dt_sub_t*ud3; v4=vs+dt_sub_t*vd3; w4=ws+dt_sub_t*wd3;
                [ud4,vd4,wd4] = trans_dot_fn(pk,qk,rk, u4,v4,w4, gxk,gyk,gzk, Tk_m, kdm_k);
                us = us + dt_sub_t/6*(ud1+2*ud2+2*ud3+ud4);
                vs = vs + dt_sub_t/6*(vd1+2*vd2+2*vd3+vd4);
                ws = ws + dt_sub_t/6*(wd1+2*wd2+2*wd3+wd4);
            end
            u_int(k+1) = us; v_int(k+1) = vs; w_int(k+1) = ws;
        end

        % Acelerômetro (sensor): pqr SIMULADO + u,v,w integrados.
        [accX, accY, accZ] = accel_eval(p_sim, q_sim, r_sim, ...
            u_int, v_int, w_int, pwm, k_T, func_T_ref, m, r_imu_hyb, dt, kdm_P);

        res.(['accX_s_' suffix]) = accX;
        res.(['accY_s_' suffix]) = accY;
        res.(['accZ_s_' suffix]) = accZ;
        res.(['full_' suffix '_ok']) = true;
    end
end

function print_R2(label, res, ~, pqr_vl, ~, acc_vl, ~, t_val, R2_func)
    fprintf('  [%s] Validação (%d-%ds):\n', label, t_val(1), t_val(2));
    if res.pqr_vl_ok
        fprintf('    R² p=%.4f | q=%.4f | r=%.4f\n', ...
            R2_func(pqr_vl(:,1), res.p_s_vl), R2_func(pqr_vl(:,2), res.q_s_vl), R2_func(pqr_vl(:,3), res.r_s_vl));
    else, fprintf('    p,q,r -> DIVERGIU\n');
    end
    if res.full_vl_ok
        fprintf('    R² AccX=%.4f | AccY=%.4f | AccZ=%.4f\n', ...
            R2_func(acc_vl(:,1), res.accX_s_vl), R2_func(acc_vl(:,2), res.accY_s_vl), R2_func(acc_vl(:,3), res.accZ_s_vl));
    end
    % Métricas extras: RMSE físico + TIC (Theil) + decomposição bias/var/cov
    rows = {};
    if res.pqr_vl_ok
        rows = [rows; {'p',pqr_vl(:,1),res.p_s_vl}; {'q',pqr_vl(:,2),res.q_s_vl}; {'r',pqr_vl(:,3),res.r_s_vl}];
    end
    if res.full_vl_ok
        rows = [rows; {'accX',acc_vl(:,1),res.accX_s_vl}; {'accY',acc_vl(:,2),res.accY_s_vl}; {'accZ',acc_vl(:,3),res.accZ_s_vl}];
    end
    if ~isempty(rows)
        fprintf('    %-5s %9s %7s   Theil bias/var/cov\n', 'canal', 'RMSE', 'TIC');
        for i = 1:size(rows,1)
            fm = fit_metrics(rows{i,2}, rows{i,3});
            fprintf('    %-5s %9.4f %7.3f   %.2f/%.2f/%.2f\n', ...
                rows{i,1}, fm.RMSE, fm.TIC, fm.U_bias, fm.U_var, fm.U_cov);
        end
    end
end

function plot_ms(t, meas, sim, R2f, ylab, ttl, leg1, last)
%PLOT_MS  Subplot medido (azul) vs simulado (vermelho tracejado) + R² no título.
    plot(t, meas, 'b-', t, sim, 'r--', 'LineWidth', 1.3);
    legend(leg1, 'Sim'); ylabel(ylab); grid on;
    title(sprintf('%s (R^2=%.3f)', ttl, R2f(meas, sim)));
    if last, xlabel('Tempo (s)'); end
end

function plot_all_results(titulo, res, ~, ~, ~, ...
    time_vl, pqr_vl, acc_vl, R2_func, ~, t_val_range, att_vl)

    lbl_vl = sprintf('VALIDAÇÃO [%d-%ds]', t_val_range(1), t_val_range(2));
    prefix = stage_prefix(titulo);   % '01_P0_', '02_Pfinal_' ou '03_Pfinal_semi_'

    img = setup_paths().images;
    fig1 = figure('Name', ['pqr - ' titulo], 'Position', [80 50 900 700], 'Visible', 'off');
    subplot(3,1,1); plot_ms(time_vl, pqr_vl(:,1), res.p_s_vl, R2_func, 'p (rad/s)', [lbl_vl ' — p'], 'Exp', false);
    subplot(3,1,2); plot_ms(time_vl, pqr_vl(:,2), res.q_s_vl, R2_func, 'q (rad/s)', [lbl_vl ' — q'], 'Exp', false);
    subplot(3,1,3); plot_ms(time_vl, pqr_vl(:,3), res.r_s_vl, R2_func, 'r (rad/s)', [lbl_vl ' — r'], 'Exp', true);
    sgtitle(['Vel. Angulares (Val) — ' titulo]); saveas(fig1, fullfile(img, [prefix 'pqr.png']));

    fig2 = figure('Name', ['Acc - ' titulo], 'Position', [120 90 900 700], 'Visible', 'off');
    subplot(3,1,1); plot_ms(time_vl, acc_vl(:,1), res.accX_s_vl, R2_func, 'AccX (m/s²)', [lbl_vl ' — AccX'], 'IMU', false);
    subplot(3,1,2); plot_ms(time_vl, acc_vl(:,2), res.accY_s_vl, R2_func, 'AccY (m/s²)', [lbl_vl ' — AccY'], 'IMU', false);
    subplot(3,1,3); plot_ms(time_vl, acc_vl(:,3), res.accZ_s_vl, R2_func, 'AccZ (m/s²)', [lbl_vl ' — AccZ'], 'IMU', true);
    sgtitle(['Acelerações (Val) — ' titulo]); saveas(fig2, fullfile(img, [prefix 'acc.png']));

    % --- Gráfico de Atitude (phi, theta, psi) ---
    if nargin >= 12 && ~isempty(att_vl) && isfield(res, 'phi_s_vl') && res.pqr_vl_ok
        fig3 = figure('Name', ['Atitude - ' titulo], 'Position', [160 130 900 700], 'Visible', 'off');
        subplot(3,1,1); plot_ms(time_vl, att_vl(:,1), res.phi_s_vl,   R2_func, '\phi (°)',   [lbl_vl ' — \phi'],   'EKF', false);
        subplot(3,1,2); plot_ms(time_vl, att_vl(:,2), res.theta_s_vl, R2_func, '\theta (°)', [lbl_vl ' — \theta'], 'EKF', false);
        subplot(3,1,3); plot_ms(time_vl, att_vl(:,3), res.psi_s_vl,   R2_func, '\psi (°)',   [lbl_vl ' — \psi'],   'EKF', true);
        sgtitle(['Atitude (Val) — ' titulo]); saveas(fig3, fullfile(img, [prefix 'att.png']));
    end
end

function plot_torques(P, time_vl, pwm_vl, func_T_ref, func_Q_ref, pqr_vl, t_val)
    k_T = P(5:8);
    k_Q = P(9:12);
    N_vl = length(time_vl);

    T_ref_vl = zeros(N_vl, 4);
    Q_ref_vl = zeros(N_vl, 4);
    for j = 1:4
        T_ref_vl(:,j) = func_T_ref(pwm_vl(:,j));
        Q_ref_vl(:,j) = func_Q_ref(pwm_vl(:,j));
    end

    Tmr = T_ref_vl .* k_T';
    Qmr = Q_ref_vl .* k_Q';

    % Braços (lidos de parameters() — fonte única de verdade)
    proj_p_pt = parameters();
    Lx_r = proj_p_pt.arms.Lx_r; Lx_l = proj_p_pt.arms.Lx_l;
    Ly_f = proj_p_pt.arms.Ly_f;    Ly_r = proj_p_pt.arms.Ly_r;

    Mx = -(Lx_r*Tmr(:,1) - Lx_l*Tmr(:,2) - Lx_l*Tmr(:,3) + Lx_r*Tmr(:,4));
    % My — ArduPilot QuadX (M1,M3=FRONT; M2,M4=REAR)
    My =  Ly_f*Tmr(:,1) - Ly_r*Tmr(:,2) + Ly_f*Tmr(:,3) - Ly_r*Tmr(:,4);
    % Mz: pares DIAGONAIS padrão ArduPilot (M1+M2 CCW, M3+M4 CW)
    Mz = Qmr(:,1) + Qmr(:,2) - Qmr(:,3) - Qmr(:,4);

    fig3 = figure('Name', 'Torques Validação', 'Position', [80 50 1000 900], 'Visible', 'off');

    subplot(4,1,1);
    plot(time_vl, Tmr(:,1), time_vl, Tmr(:,2), time_vl, Tmr(:,3), time_vl, Tmr(:,4), 'LineWidth', 1.2);
    legend('T_1','T_2','T_3','T_4'); ylabel('Empuxo (N)'); grid on;
    title(sprintf('Empuxo individual (k_T=[%.3f, %.3f, %.3f, %.3f])', k_T(1), k_T(2), k_T(3), k_T(4)));

    subplot(4,1,2);
    yyaxis left;
    plot(time_vl, Mx, 'r-', 'LineWidth', 1.3);
    ylabel('Mx (N·m)'); set(gca,'YColor','r');
    yyaxis right;
    plot(time_vl, pqr_vl(:,1), 'b-', 'LineWidth', 1);
    ylabel('p (rad/s)'); set(gca,'YColor','b');
    grid on; legend({'Mx modelo','p medido'}, 'Location','best');
    title(sprintf('Roll: Mx (modelo)  vs  p (medido)   (Lx=%.3f)', Lx_r));
    align_yyaxis_zero(gca);

    subplot(4,1,3);
    yyaxis left;
    plot(time_vl, My, 'r-', 'LineWidth', 1.3);
    ylabel('My (N·m)'); set(gca,'YColor','r');
    yyaxis right;
    plot(time_vl, pqr_vl(:,2), 'b-', 'LineWidth', 1);
    ylabel('q (rad/s)'); set(gca,'YColor','b');
    grid on; legend({'My modelo','q medido'}, 'Location','best');
    title(sprintf('Pitch: My (modelo)  vs  q (medido)   (Ly_f=%.3f, Ly_r=%.3f)', Ly_f, Ly_r));
    align_yyaxis_zero(gca);

    subplot(4,1,4);
    yyaxis left;
    plot(time_vl, Mz, 'r-', 'LineWidth', 1.3);
    ylabel('Mz (N·m)'); set(gca,'YColor','r');
    yyaxis right;
    plot(time_vl, pqr_vl(:,3), 'b-', 'LineWidth', 1);
    ylabel('r (rad/s)'); set(gca,'YColor','b');
    grid on; legend({'Mz modelo','r medido'}, 'Location','best');
    xlabel('Tempo (s)');
    title(sprintf('Yaw: Mz (modelo)  vs  r (medido)   (k_Q=[%.3f, %.3f, %.3f, %.3f])', k_Q(1), k_Q(2), k_Q(3), k_Q(4)));
    align_yyaxis_zero(gca);

    sgtitle(sprintf('Torques na Validação [%d-%ds]  (CG fixo do CAD)', t_val(1), t_val(2)));
    saveas(fig3, fullfile(setup_paths().images, '05_torques.png'));

    fprintf('\n  --- Torques (validação) | CG fixo do CAD ---\n');
    fprintf('    Braços: Lx_r=%.4f Lx_l=%.4f Ly_f=%.4f Ly_r=%.4f\n', Lx_r, Lx_l, Ly_f, Ly_r);
    fprintf('    Mx: min=%.4f  max=%.4f  std=%.4f\n', min(Mx), max(Mx), std(Mx));
    fprintf('    My: min=%.4f  max=%.4f  std=%.4f\n', min(My), max(My), std(My));
    fprintf('    Mz: min=%.4f  max=%.4f  std=%.4f\n', min(Mz), max(Mz), std(Mz));
end

function plot_forces(P, time_vl, pwm_vl, att_vl, acc_vl, func_T_ref, constants_sim, t_val, res)
    k_T = P(5:8);
    % drag e Bz removidos do modelo
    Xu_m = 0; Yv_m = 0; Zw_m = 0; Bz = 0;
    m = constants_sim.m;
    g = constants_sim.g;
    N_vl = length(time_vl);

    % Atitude em radianos
    phi   = deg2rad(att_vl(:,1));
    theta = deg2rad(att_vl(:,2));

    % Empuxo individual e total
    Ti = zeros(N_vl, 4);
    for j = 1:4
        Ti(:,j) = k_T(j) * func_T_ref(pwm_vl(:,j));
    end
    T_total = sum(Ti, 2);

    % Componentes da gravidade no body frame (por unidade de massa)
    gx = -g * sin(theta);
    gy =  g * cos(theta) .* sin(phi);
    gz =  g * cos(theta) .* cos(phi);

    % Forças por unidade de massa (acelerações)
    Fz_m = -T_total/m + gz;       % thrust + gravidade z

    % =====================================================================
    %  Figura 1: Forças translacionais com modelo simulado
    % =====================================================================
    fig4 = figure('Name', 'Forças Translacionais', 'Position', [80 50 1000 900], 'Visible', 'off');

    % --- Subplot 1: Empuxo total vs peso ---
    subplot(4,1,1);
    plot(time_vl, T_total, 'r-', 'LineWidth', 1.3);
    hold on;
    yline(m*g, 'k--', 'LineWidth', 1.2);
    hold off;
    legend('T_{total}', sprintf('mg = %.2f N', m*g));
    ylabel('Força (N)'); grid on;
    title(sprintf('Empuxo total vs Peso  (T_{med}=%.2f N, mg=%.2f N, ratio=%.3f)', ...
        mean(T_total), m*g, mean(T_total)/(m*g)));

    % --- Subplot 2: AccX medido vs simulado ---
    subplot(4,1,2);
    plot(time_vl, acc_vl(:,1), 'b-', 'LineWidth', 1.0);
    hold on;
    if res.full_vl_ok
        plot(time_vl, res.accX_s_vl, 'r--', 'LineWidth', 1.3);
    end
    plot(time_vl, gx, 'Color', [1.0 0.6 0.0], 'LineStyle', ':', 'LineWidth', 1.2);
    hold off;
    if res.full_vl_ok
        legend('AccX_{IMU} (medido)', 'AccX_{modelo}', 'gx = -g sin\theta');
    else
        legend('AccX_{IMU} (medido)', 'gx = -g sin\theta');
    end
    ylabel('m/s^2'); grid on;
    title('Eixo X: força específica medida vs modelo');

    % --- Subplot 3: AccY medido vs simulado ---
    subplot(4,1,3);
    plot(time_vl, acc_vl(:,2), 'b-', 'LineWidth', 1.0);
    hold on;
    if res.full_vl_ok
        plot(time_vl, res.accY_s_vl, 'r--', 'LineWidth', 1.3);
    end
    plot(time_vl, gy, 'Color', [1.0 0.6 0.0], 'LineStyle', ':', 'LineWidth', 1.2);
    hold off;
    if res.full_vl_ok
        legend('AccY_{IMU} (medido)', 'AccY_{modelo}', 'gy = g cos\theta sin\phi');
    else
        legend('AccY_{IMU} (medido)', 'gy = g cos\theta sin\phi');
    end
    ylabel('m/s^2'); grid on;
    title('Eixo Y: força específica medida vs modelo');

    % --- Subplot 4: AccZ medido vs simulado ---
    subplot(4,1,4);
    plot(time_vl, acc_vl(:,3), 'b-', 'LineWidth', 1.0);
    hold on;
    if res.full_vl_ok
        plot(time_vl, res.accZ_s_vl, 'r--', 'LineWidth', 1.3);
    end
    plot(time_vl, Fz_m, 'Color', [1.0 0.6 0.0], 'LineStyle', ':', 'LineWidth', 1.2);
    yline(-g, 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 0.8);
    hold off;
    if res.full_vl_ok
        legend('AccZ_{IMU} (medido)', 'AccZ_{modelo}', 'Fz/m = -T/m+gz', '-g');
    else
        legend('AccZ_{IMU} (medido)', 'Fz/m = -T/m+gz', '-g');
    end
    ylabel('m/s^2'); xlabel('Tempo (s)'); grid on;
    title(sprintf('Eixo Z: força específica medida vs modelo  (Bz=%.3f, Zw=%.3f)', Bz, Zw_m));

    sgtitle(sprintf('Forças Translacionais — Validação [%d-%ds]  (m=%.1f kg)', t_val(1), t_val(2), m));
    saveas(fig4, fullfile(setup_paths().images, '06_forcas.png'));

    % =====================================================================
    %  Figura 2: Análise dos PWMs por motor
    % =====================================================================
    fig5 = figure('Name', 'Análise PWM', 'Position', [100 50 1000 700], 'Visible', 'off');

    % Subplot 1: PWM de cada motor ao longo do tempo
    subplot(3,1,1);
    plot(time_vl, pwm_vl(:,1), 'r-', 'LineWidth', 1.0); hold on;
    plot(time_vl, pwm_vl(:,2), 'b-', 'LineWidth', 1.0);
    plot(time_vl, pwm_vl(:,3), 'g-', 'LineWidth', 1.0);
    plot(time_vl, pwm_vl(:,4), 'm-', 'LineWidth', 1.0);
    hold off;
    legend('M1 (FR,CW)', 'M2 (RL,CW)', 'M3 (FL,CCW)', 'M4 (RR,CCW)');
    ylabel('PWM (\mus)'); grid on;
    title('PWM por motor');

    % Subplot 2: Empuxo individual (com k_T)
    subplot(3,1,2);
    plot(time_vl, Ti(:,1), 'r-', 'LineWidth', 1.0); hold on;
    plot(time_vl, Ti(:,2), 'b-', 'LineWidth', 1.0);
    plot(time_vl, Ti(:,3), 'g-', 'LineWidth', 1.0);
    plot(time_vl, Ti(:,4), 'm-', 'LineWidth', 1.0);
    yline(m*g/4, 'k--', 'LineWidth', 1.0);
    hold off;
    legend(sprintf('T1 (k_T=%.3f)', k_T(1)), sprintf('T2 (k_T=%.3f)', k_T(2)), ...
           sprintf('T3 (k_T=%.3f)', k_T(3)), sprintf('T4 (k_T=%.3f)', k_T(4)), 'mg/4');
    ylabel('Empuxo (N)'); grid on;
    title('Empuxo por motor (T_i = k_{Ti} \cdot T_{ref}(PWM_i))');

    % Subplot 3: T_ref cru (sem k_T) — para ver o que o FC realmente comanda
    T_ref_raw = zeros(N_vl, 4);
    for j = 1:4
        T_ref_raw(:,j) = func_T_ref(pwm_vl(:,j));
    end
    subplot(3,1,3);
    plot(time_vl, T_ref_raw(:,1), 'r-', 'LineWidth', 1.0); hold on;
    plot(time_vl, T_ref_raw(:,2), 'b-', 'LineWidth', 1.0);
    plot(time_vl, T_ref_raw(:,3), 'g-', 'LineWidth', 1.0);
    plot(time_vl, T_ref_raw(:,4), 'm-', 'LineWidth', 1.0);
    yline(m*g/4, 'k--', 'LineWidth', 1.0);
    hold off;
    T_ref_total = sum(T_ref_raw, 2);
    legend(sprintf('T_{ref1} (med=%.2f)', mean(T_ref_raw(:,1))), ...
           sprintf('T_{ref2} (med=%.2f)', mean(T_ref_raw(:,2))), ...
           sprintf('T_{ref3} (med=%.2f)', mean(T_ref_raw(:,3))), ...
           sprintf('T_{ref4} (med=%.2f)', mean(T_ref_raw(:,4))), 'mg/4');
    ylabel('Empuxo ref (N)'); xlabel('Tempo (s)'); grid on;
    title(sprintf('T_{ref}(PWM) cru (sem k_T) — Total medio=%.2f N vs mg=%.2f N, ratio=%.3f', ...
        mean(T_ref_total), m*g, mean(T_ref_total)/(m*g)));

    sgtitle(sprintf('Análise PWM e Empuxo — Validação [%d-%ds]', t_val(1), t_val(2)));
    saveas(fig5, fullfile(setup_paths().images, '04_pwm.png'));

    % Estatísticas
    fprintf('\n  --- Forças (validação) ---\n');
    fprintf('    T_total: média=%.3f N  (mg=%.3f N, ratio=%.4f)\n', mean(T_total), m*g, mean(T_total)/(m*g));
    fprintf('    Ti médios: [%.2f, %.2f, %.2f, %.2f] N\n', mean(Ti(:,1)), mean(Ti(:,2)), mean(Ti(:,3)), mean(Ti(:,4)));
    fprintf('    T_ref cru (sem k_T): [%.2f, %.2f, %.2f, %.2f] N  (total=%.2f)\n', ...
        mean(T_ref_raw(:,1)), mean(T_ref_raw(:,2)), mean(T_ref_raw(:,3)), mean(T_ref_raw(:,4)), mean(T_ref_total));
    fprintf('    PWM médios: [%.0f, %.0f, %.0f, %.0f] us\n', mean(pwm_vl(:,1)), mean(pwm_vl(:,2)), mean(pwm_vl(:,3)), mean(pwm_vl(:,4)));
    fprintf('    PWM std:    [%.1f, %.1f, %.1f, %.1f] us\n', std(pwm_vl(:,1)), std(pwm_vl(:,2)), std(pwm_vl(:,3)), std(pwm_vl(:,4)));
    fprintf('    Fz/m média=%.4f  (esperado ~0 em hover: -T/m+gz ~0)\n', mean(Fz_m));
    fprintf('    Xu_m=%.4f  Yv_m=%.4f  Zw_m=%.4f  Bz=%.4f\n', Xu_m, Yv_m, Zw_m, Bz);

    % Análise: assimetria constante (CG/motor) vs variável (vento)
    fprintf('\n  --- Análise de assimetria PWM ---\n');
    pwm_mean = mean(pwm_vl);
    pwm_std  = std(pwm_vl);
    cv = pwm_std ./ pwm_mean * 100;  % coeficiente de variação (%)
    fprintf('    Coeficiente de variação PWM: [%.1f%%, %.1f%%, %.1f%%, %.1f%%]\n', cv(1), cv(2), cv(3), cv(4));
    fprintf('    (CV baixo = patamar constante → diferença CG/motor)\n');
    fprintf('    (CV alto  = variação temporal → vento/perturbação)\n');

    % Assimetria roll: (M1+M4) vs (M2+M3)
    pwm_right = pwm_vl(:,1) + pwm_vl(:,4);  % motores direita
    pwm_left  = pwm_vl(:,2) + pwm_vl(:,3);  % motores esquerda
    roll_asym = pwm_right - pwm_left;
    fprintf('    Assimetria roll (dir-esq): média=%.1f us  std=%.1f us\n', mean(roll_asym), std(roll_asym));

    % Assimetria pitch: (M1+M3) vs (M2+M4)
    pwm_front = pwm_vl(:,1) + pwm_vl(:,3);  % motores frente
    pwm_rear  = pwm_vl(:,2) + pwm_vl(:,4);  % motores traseira
    pitch_asym = pwm_front - pwm_rear;
    fprintf('    Assimetria pitch (frente-trás): média=%.1f us  std=%.1f us\n', mean(pitch_asym), std(pitch_asym));

    % Assimetria yaw: (M1+M2) CW vs (M3+M4) CCW
    pwm_cw  = pwm_vl(:,1) + pwm_vl(:,2);
    pwm_ccw = pwm_vl(:,3) + pwm_vl(:,4);
    yaw_asym = pwm_cw - pwm_ccw;
    fprintf('    Assimetria yaw (CW-CCW): média=%.1f us  std=%.1f us\n', mean(yaw_asym), std(yaw_asym));

    if mean(abs(std(roll_asym))) < 5 && mean(abs(std(pitch_asym))) < 5
        fprintf('    → Assimetrias ESTÁVEIS: predomina CG offset / diferença de motores\n');
    else
        fprintf('    → Assimetrias VARIÁVEIS: predomina vento / perturbação externa\n');
    end
end

function p = stage_prefix(titulo)
% Mapeia título da etapa para prefixo numerado dos arquivos de imagem.
%   titulo = 'CHUTE_INICIAL_P0'             -> '01_P0_'
%   titulo = 'P_final'                      -> '02_Pfinal_'
%   titulo = 'P_final (hibrido)'            -> '03_Pfinal_hib_'
%   titulo = 'P_final (semi-acoplado)'      -> '04_Pfinal_semi_'
% Padrão: <ordem>_<stage>_<variável>.png
    t = lower(titulo);
    if contains(t, 'p0') || contains(t, 'chute')
        p = '01_P0_';
    elseif contains(t, 'hibrido') || contains(t, 'híbrido')
        p = '03_Pfinal_hib_';
    elseif contains(t, 'semi')
        p = '04_Pfinal_semi_';
    else
        p = '02_Pfinal_';
    end
end

function [aLp, aLb, aMq, aMa, aNr, aNb] = aero_terms_local(P, time)
%AERO_TERMS_LOCAL  Coeficientes aerodinâmicos por amostra para os simuladores
% auxiliares (híbrido, semi-acoplado). V, v, w vêm do appdata 'aero_vel' posto
% no início do script (velocidade estimada de bordo). Sem appdata, ou com P de
% 15 elementos, devolve zeros — o modelo volta a ser o rotacional puro.
    N = numel(time);
    z = zeros(N,1);  aLp=z; aLb=z; aMq=z; aMa=z; aNr=z; aNb=z;
    av = getappdata(0, 'aero_vel');
    if numel(P) < 22 || isempty(av) || ~isstruct(av), return; end
    V = interp1(av.t, av.V, time(:), 'linear', 0);
    v = interp1(av.t, av.v, time(:), 'linear', 0);
    w = interp1(av.t, av.w, time(:), 'linear', 0);
    V(~isfinite(V)) = 0;  v(~isfinite(v)) = 0;  w(~isfinite(w)) = 0;
    kA = aero_gains();
    aLp = kA.Lp*P(17)*V;   aLb = kA.Lb*P(18)*V.*v;
    aMq = kA.Mq*P(19)*V;   aMa = kA.Ma*P(20)*V.*w;
    aNr = kA.Nr*P(21)*V;   aNb = kA.Nb*P(22)*V.*v;
    if numel(P) >= 25
        if isappdata(0,'hu_form')
            qS = 0.5*kA.rho*V.^2*kA.S;
            aLb = aLb - qS*kA.b*P(23);  aMa = aMa - qS*kA.c*P(24);  aNb = aNb - qS*kA.b*P(25);
        else
            aLb = aLb - P(23)*v;  aMa = aMa - P(24)*w;  aNb = aNb - P(25)*v;
        end
    end
end
