"""make_diag.py — gera identify_plant_diag.m a partir de identify_plant.m.

O script de diagnóstico é uma CÓPIA do oficial com overrides (struct DIAG) e
saídas redirecionadas para outputs/diag/<tag>/. Nunca edite identify_plant_diag.m
à mão: edite o oficial e rode  python3 make_diag.py  para regenerar.

Cada âncora abaixo deve ocorrer exatamente `count` vezes no oficial (assert), o
que protege contra edições silenciosas no arquivo-fonte.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "identify_plant.m")
DST = os.path.join(HERE, "identify_plant_diag.m")
s = open(SRC, encoding="utf-8").read()

def rep(anchor, new, count=1):
    global s
    n = s.count(anchor)
    assert n == count, f"âncora encontrada {n}x (esperado {count}): {anchor[:70]!r}"
    s = s.replace(anchor, new)

HEADER = """% identify_plant_diag.m — CÓPIA DE DIAGNÓSTICO de identify_plant.m  (GERADO por make_diag.py)
% =========================================================================
% Mesmo pipeline (EEM → OEM progressivo → CRB → validação → brancura), mas
% (i) aceita overrides pela struct DIAG (definida no workspace antes de rodar) e
% (ii) grava TUDO em outputs/diag/<DIAG.tag>/ — NUNCA sobrescreve os outputs
% oficiais (P_identified.mat, outputs/images/*).
%
% Campos de DIAG (todos opcionais):
%   tag        nome da rodada (pasta de saída)                 default 'base'
%   tau_m      lag 1ª ordem do motor [s]      (override de parameters().motor.tau_m)
%   delay_pwm  atraso puro em AMOSTRAS (0,1 s) (override de parameters().motor.delay_s)
%   P0_damp    chute inicial [cp cq cr]
%   lock_damp  TRAVA [cp cq cr] nesses valores (lb=ub)
%   lock_kT    TRAVA k_T1..4 nesses valores (escalar ou 4×1)  — testa "cT>1 é empuxo ou momento?"
%   lb_J/ub_J  limites das inércias [Jx Jy Jz Jxz] (libera Jz/Jxz, hoje travados)
%   t_trains   janelas de treino {[a b];...}   t_val  janela de validação [a b]
%   drag_cd    C_d do arrasto induzido do rotor (Beard 14.4.2): f_d = −T·C_d·[u;v;0].
%              Entra na EOM translacional e no acelerômetro. 0 = desliga (padrão).
%   damp_mode  'const' (padrão) | 'aero' | 'hybrid' — TESTE: amortecimento −(c0 + D·V(t))·p com V(t)
%              da estimate_velocity (atitude+barômetro). 'aero': c0 = 0, D identificado (Cl_p, Cm_q,
%              Cn_r ficam em summary.aero_coef); 'hybrid': c0 = a priori de rotor (prior_damping), D idem.
%
% Ao final salva outputs/diag/<tag>/summary.mat com P0, P_eem, P_final, CRB
% (nos 3 pontos, pelo custo OEM e pelo custo EEM), R² de validação, brancura, tempo.
% Driver: run_diag.m.  Relatório: diag_report.m.
% =========================================================================
"""
s = HEADER + s

rep("paths = setup_paths();\n",
    "paths = setup_paths();\n"
    "\n%% 0b. DIAG — overrides e redirecionamento de saídas (ver cabeçalho)\n"
    "if ~exist('DIAG','var') || ~isstruct(DIAG), DIAG = struct(); end\n"
    "if ~isfield(DIAG,'tag'), DIAG.tag = 'base'; end\n"
    "diag_dir = fullfile(paths.outputs, 'diag', DIAG.tag);\n"
    "if ~exist(diag_dir,'dir'), mkdir(diag_dir); end\n"
    "paths.images = diag_dir;  setappdata(0,'diag_img_dir',diag_dir);\n"
    "t_diag = tic;  white_pct = nan(1,3);  se = nan(15,1);  se_P0 = nan(15,1);  se_eem = nan(15,1);  se_P0_e = nan(15,1);  se_eem_e = nan(15,1);\n"
    "fprintf('\\n>>> DIAG rodada \"%s\" → %s\\n', DIAG.tag, diag_dir); disp(DIAG);\n", count=1)

rep("setup_paths().images", "getappdata(0,'diag_img_dir')", count=4)

rep("t_val    = [610, 630];",
    "t_val    = [610, 630];\n"
    "if isfield(DIAG,'t_trains'), t_trains = DIAG.t_trains; end   % DIAG\n"
    "if isfield(DIAG,'t_val'),    t_val    = DIAG.t_val;    end   % DIAG\n"
    "% (linha original:) t_val")

# cadeia de motor: overrides de tau_m / delay
rep("[rpm_lag, func_T_ref, func_Q_ref, mc_info] = motor_chain(t_common(:), ...\n"
    "    [pwm1_interp(:), pwm2_interp(:), pwm3_interp(:), pwm4_interp(:)]);",
    "mc_args = {};                                                     % DIAG\n"
    "if isfield(DIAG,'tau_m'),     mc_args = [mc_args, {'tau_m',   DIAG.tau_m}];         end\n"
    "if isfield(DIAG,'delay_pwm'), mc_args = [mc_args, {'delay_s', DIAG.delay_pwm*dt}]; end\n"
    "[rpm_lag, func_T_ref, func_Q_ref, mc_info] = motor_chain(t_common(:), ...\n"
    "    [pwm1_interp(:), pwm2_interp(:), pwm3_interp(:), pwm4_interp(:)], mc_args{:});")

rep("param_names = proj_params.param_names;\n",
    "param_names = proj_params.param_names;\n"
    "if isfield(DIAG,'P0_damp')                                    % DIAG: chute inicial\n"
    "    P0(13:15) = DIAG.P0_damp(:);\n"
    "    fprintf('  >>> DIAG: P0(cp,cq,cr) = [%.3f %.3f %.3f]\\n', P0(13:15));\n"
    "end\n"
    "if isfield(DIAG,'lock_damp')                                  % DIAG: trava cp,cq,cr\n"
    "    P0(13:15) = DIAG.lock_damp(:); lb(13:15) = DIAG.lock_damp(:); ub(13:15) = DIAG.lock_damp(:);\n"
    "    fprintf('  >>> DIAG: cp,cq,cr TRAVADOS em [%.3f %.3f %.3f]\\n', P0(13:15));\n"
    "end\n"
    "if isfield(DIAG,'lock_kT')                                    % DIAG: trava k_T1..4\n"
    "    v = DIAG.lock_kT(:); if numel(v)==1, v = v*ones(4,1); end\n"
    "    P0(5:8) = v; lb(5:8) = v; ub(5:8) = v;\n"
    "    fprintf('  >>> DIAG: k_T TRAVADOS em [%.3f %.3f %.3f %.3f]\\n', v);\n"
    "end\n"
    "if isfield(DIAG,'lb_J'), lb(1:4) = DIAG.lb_J(:); end          % DIAG: limites das inércias\n"
    "if isfield(DIAG,'ub_J'), ub(1:4) = DIAG.ub_J(:); end\n"
    "if isfield(DIAG,'lb_J') || isfield(DIAG,'ub_J')\n"
    "    fprintf('  >>> DIAG: limites J  lb = [%.5f %.5f %.5f %.5f]\\n', lb(1:4));\n"
    "    fprintf('                       ub = [%.5f %.5f %.5f %.5f]\\n', ub(1:4));\n"
    "end\n")

rep("%% ========================================================================\n%  6. VALIDAÇÃO COM CHUTE INICIAL (P0)",
    "%% 5b. DIAG — amortecimento ∝ V (damp_mode)\n"
    "if isfield(DIAG,'drag_cd') && DIAG.drag_cd > 0\n"
    "    setappdata(0,'diag_drag', DIAG.drag_cd);\n"
    "    fprintf('  >>> DIAG: arrasto induzido do rotor C_d = %.4f (g·C_d = %.3f 1/s)\\n', DIAG.drag_cd, 9.81*DIAG.drag_cd);\n"
    "elseif isappdata(0,'diag_drag'), rmappdata(0,'diag_drag');\n"
    "end\n"
    "damp_mode = 'const'; if isfield(DIAG,'damp_mode'), damp_mode = DIAG.damp_mode; end\n"
    "if any(strcmp(damp_mode,{'aero','hybrid'}))\n"
    "    VEst = estimate_velocity(L, t_common(:));                 % velocidade estimada (atitude+baro)\n"
    "    Vgrid = VEst.V;  Vgrid(~isfinite(Vgrid)) = 0;\n"
    "    if strcmp(damp_mode,'hybrid'), c0d = [1.19; 1.20; 0.08]; else, c0d = [0;0;0]; end   % a priori de rotor\n"
    "    for s = 1:n_seg, segs{s}.V = interp1(t_common(:), Vgrid, segs{s}.time(:), 'linear', 'extrap'); end\n"
    "    V_eem = []; for s = 1:n_seg, V_eem = [V_eem; segs{s}.V(:)]; end %#ok<AGROW>\n"
    "    setappdata(0,'diag_damp', struct('mode',damp_mode,'t',t_common(:),'V',Vgrid,'V_eem',V_eem,'c0',c0d));\n"
    "    fprintf('  >>> DIAG: damp_mode = %s | V estimada p50 %.2f p95 %.2f m/s | c0 = [%.2f %.2f %.2f]\\n', damp_mode, median(Vgrid(Vgrid>0)), prctile(Vgrid(Vgrid>0),95), c0d);\n"
    "    if ~isfield(DIAG,'P0_damp')   % chute do coeficiente ∝ V pela AVL: D = ρ S b² |Cl_p|/(4Jx) etc.\n"
    "        rho_ = 1.225; S_ = proj_params.wing.S; b_ = proj_params.wing.b; c_ = proj_params.wing.c;\n"
    "        P0(13:15) = [rho_*S_*b_^2*0.406/(4*P0(1)); rho_*S_*c_^2*8.96/(2*P0(2)); rho_*S_*b_^2*0.070/(4*P0(3))];\n"
    "        fprintf('  >>> DIAG: P0 do amortecimento ∝ V pela AVL = [%.3f %.3f %.3f]\\n', P0(13:15));\n"
    "    end\n"
    "else\n"
    "    if isappdata(0,'diag_damp'), rmappdata(0,'diag_damp'); end\n"
    "end\n\n"
    "%% ========================================================================\n%  6. VALIDAÇÃO COM CHUTE INICIAL (P0)")

rep("params_file = fullfile(paths.outputs, 'P_identified.mat');",
    "params_file = fullfile(diag_dir, 'P_identified.mat');   % DIAG: nunca o oficial")

rep("        e_seg = oem_ms_cost_func(P, sg.pqr, sg.acc, sg.att_rad, ...\n"
    "            sg.T_ref, sg.Q_ref, m, g, dt, N, ...\n"
    "            win_starts, win_ends, weights_pqr, weights_acc, cost_mode);",
    "        Vseg = ones(N,1); if isfield(sg,'V'), Vseg = sg.V(:); end     % DIAG damp_mode\n"
    "        e_seg = oem_ms_cost_func(P, sg.pqr, sg.acc, sg.att_rad, ...\n"
    "            sg.T_ref, sg.Q_ref, m, g, dt, N, ...\n"
    "            win_starts, win_ends, weights_pqr, weights_acc, cost_mode, Vseg);")

rep("function e = oem_ms_cost_func(P, pqr, acc, att_rad, T_ref, Q_ref, m, g, dt, N, ...\n"
    "    win_starts, win_ends, weights_pqr, weights_acc, cost_mode)\n\n"
    "    if nargin < 15, cost_mode = 'full'; end",
    "function e = oem_ms_cost_func(P, pqr, acc, att_rad, T_ref, Q_ref, m, g, dt, N, ...\n"
    "    win_starts, win_ends, weights_pqr, weights_acc, cost_mode, Vseg)\n\n"
    "    if nargin < 15, cost_mode = 'full'; end\n"
    "    if nargin < 16 || isempty(Vseg), Vseg = ones(N,1); end          % DIAG damp_mode\n"
    "    dd_ = getappdata(0,'diag_damp');  c0_ = [0;0;0];\n"
    "    if ~isempty(dd_) && isstruct(dd_) && any(strcmp(dd_.mode,{'aero','hybrid'})), c0_ = dd_.c0(:); else, Vseg = ones(N,1); end")

# há DUAS ocorrências: (1) no oem_ms_cost_func e (2) no simulate_full_hybrid
rep("    Dp = P(13); Dq = P(14); Dr = P(15);\n    % Bp/Bq/Br, Xu/Yv/Zw, Bz removidos do modelo",
    "    Dp0 = P(13); Dq0 = P(14); Dr0 = P(15);   % (DIAG: escalados por V dentro do loop)\n"
    "    Dp = Dp0; Dq = Dq0; Dr = Dr0;\n"
    "    % Bp/Bq/Br, Xu/Yv/Zw, Bz removidos do modelo", count=2)
# no simulate_full_hybrid o loop é 'for k = 1:N-1' com Tmr = ...; escala por V(time(k)) via appdata
rep("        for k = 1:N-1\n            Tmr = k_T(:)' .* T_ref(k,:);\n            Qmr = k_Q(:)' .* Q_ref(k,:);\n            % Momentos — ArduPilot QuadX (M1,M3=FRONT; M2,M4=REAR)",
    "        ddh_ = getappdata(0,'diag_damp');   % DIAG damp_mode (validação híbrida)\n"
    "        useV_ = ~isempty(ddh_) && isstruct(ddh_) && any(strcmp(ddh_.mode,{'aero','hybrid'}));\n"
    "        for k = 1:N-1\n"
    "            if useV_, Vk_ = interp1(ddh_.t, ddh_.V, time(k), 'linear', 'extrap'); Dp = ddh_.c0(1)+Dp0*Vk_; Dq = ddh_.c0(2)+Dq0*Vk_; Dr = ddh_.c0(3)+Dr0*Vk_; end\n"
    "            Tmr = k_T(:)' .* T_ref(k,:);\n            Qmr = k_Q(:)' .* Q_ref(k,:);\n            % Momentos — ArduPilot QuadX (M1,M3=FRONT; M2,M4=REAR)")

rep("        for k = i_s:i_e-1\n            % Momentos no início (k) e no fim (k+1) do passo",
    "        for k = i_s:i_e-1\n"
    "            Dp = c0_(1) + Dp0*Vseg(k);  Dq = c0_(2) + Dq0*Vseg(k);  Dr = c0_(3) + Dr0*Vseg(k);   % DIAG damp_mode\n"
    "            % Momentos no início (k) e no fim (k+1) do passo")

rep("    if insideW >= 90, verdictW",
    "    white_pct(cW) = insideW;   % DIAG\n"
    "    if insideW >= 90, verdictW")

# CR nos 3 pontos (custo OEM e custo EEM)
rep("    se    = sqrt(max(diag(Cov), 0));\n    at_lb = abs(P_final(:) - lb(:)) <= 1e-6*max(abs(lb(:)),1);",
    "    se    = sqrt(max(diag(Cov), 0));\n"
    "    % --- DIAG: mesmo cálculo em Θ0 e Θ_EEM (custo OEM) e pelo custo do EEM ---\n"
    "    se_P0  = diag_se_at(cost_se, P0,    lb, ub, opts_se, n_params);\n"
    "    se_eem = diag_se_at(cost_se, P_eem, lb, ub, opts_se, n_params);\n"
    "    cr = @(sv,Pv) 100*sv(:)./max(abs(Pv(:)),1e-9);\n"
    "    CR_P0 = cr(se_P0,P0);  CR_eem = cr(se_eem,P_eem);  CR_oem = cr(se,P_final);\n"
    "    fprintf('\\n  --- Tabela 5.3 reproduzida (valor | CR%% pelo custo OEM nos 3 pontos) ---\\n');\n"
    "    fprintf('  %-6s | %10s %6s | %10s %6s | %10s %6s\\n', 'param', 'Θ0','CR%','Θ_EEM','CR%','Θ_OEM','CR%');\n"
    "    for i = 1:n_params\n"
    "        fprintf('  %-6s | %10.5f %6.1f | %10.5f %6.1f | %10.5f %6.1f\\n', param_names{i}, ...\n"
    "            P0(i), CR_P0(i), P_eem(i), CR_eem(i), P_final(i), CR_oem(i));\n"
    "    end\n"
    "    se_P0_e  = diag_se_at(cost_eem, P0(1:n_rot),    lb(1:n_rot), ub(1:n_rot), opts_se, n_rot);\n"
    "    se_eem_e = diag_se_at(cost_eem, P_eem(1:n_rot), lb(1:n_rot), ub(1:n_rot), opts_se, n_rot);\n"
    "    CR_P0_e = cr(se_P0_e,P0(1:n_rot));  CR_eem_e = cr(se_eem_e,P_eem(1:n_rot));\n"
    "    fprintf('\\n  --- CR%% pelo custo do EEM (Θ0 | Θ_EEM) e pelo custo do OEM (Θ_OEM) ---\\n');\n"
    "    for i = 1:n_params\n"
    "        fprintf('  %-6s | %10.5f %6.1f | %10.5f %6.1f | %10.5f %6.1f\\n', param_names{i}, ...\n"
    "            P0(i), CR_P0_e(i), P_eem(i), CR_eem_e(i), P_final(i), CR_oem(i));\n"
    "    end\n"
    "    at_lb = abs(P_final(:) - lb(:)) <= 1e-6*max(abs(lb(:)),1);")

rep("disp('Script finalizado.');\n",
    "%% 11. DIAG — resumo da rodada\n"
    "r2v = @(res, ok) iff_r2(res, ok, pqr_vl, acc_vl, R2_func);\n"
    "summary = struct();\n"
    "summary.tag = DIAG.tag;  summary.DIAG = DIAG;  summary.motor = mc_info;\n"
    "summary.param_names = {param_names{:}};\n"
    "summary.P0 = P0(:);  summary.P_eem = P_eem(:);  summary.P_final = P_final(:);\n"
    "summary.lb = lb(:);  summary.ub = ub(:);  summary.se = se(:);\n"
    "summary.CR_pct = 100*se(:)./max(abs(P_final(:)),1e-9);\n"
    "summary.CR_pct_P0  = 100*se_P0(:)./max(abs(P0(:)),1e-9);\n"
    "summary.CR_pct_eem = 100*se_eem(:)./max(abs(P_eem(:)),1e-9);\n"
    "summary.CR_pct_P0_eemcost  = 100*se_P0_e(:)./max(abs(P0(:)),1e-9);\n"
    "summary.CR_pct_eem_eemcost = 100*se_eem_e(:)./max(abs(P_eem(:)),1e-9);\n"
    "summary.best_stage = best_stage_name;\n"
    "summary.R2_P0    = r2v(res_P0, 'pqr');      summary.R2acc_P0    = r2v(res_P0, 'acc');\n"
    "summary.R2_final = r2v(res_final, 'pqr');   summary.R2acc_final = r2v(res_final, 'acc');\n"
    "summary.R2_hyb   = r2v(res_hyb_final, 'pqr');\n"
    "summary.white_pct = white_pct;  summary.t_trains = t_trains;  summary.t_val = t_val;\n"
    "summary.damp_mode = damp_mode;\n"
    "summary.drag_cd = 0; if isfield(DIAG,'drag_cd'), summary.drag_cd = DIAG.drag_cd; end\n"
    "if isappdata(0,'diag_drag'), rmappdata(0,'diag_drag'); end\n"
    "if any(strcmp(damp_mode,{'aero','hybrid'}))   % D → coeficientes adimensionais (AVL: Cl_p −0,406, Cm_q −8,96, Cn_r −0,070)\n"
    "    rho_ = 1.225; S_ = proj_params.wing.S; b_ = proj_params.wing.b; c_ = proj_params.wing.c;\n"
    "    summary.aero_coef = struct('Clp', -4*P_final(1)*P_final(13)/(rho_*S_*b_^2), ...\n"
    "        'Cmq', -2*P_final(2)*P_final(14)/(rho_*S_*c_^2), 'Cnr', -4*P_final(3)*P_final(15)/(rho_*S_*b_^2), 'c0', c0d(:)');\n"
    "    fprintf('  >>> DIAG damp_mode=%s: Cl_p = %.3f | Cm_q = %.2f | Cn_r = %.3f  (AVL: -0.406 | -8.96 | -0.070)\\n', damp_mode, summary.aero_coef.Clp, summary.aero_coef.Cmq, summary.aero_coef.Cnr);\n"
    "    rmappdata(0,'diag_damp');\n"
    "end\n"
    "summary.elapsed_s = toc(t_diag);\n"
    "save(fullfile(diag_dir,'summary.mat'), 'summary');\n"
    "fprintf('\\n>>> DIAG \"%s\": cp=%.3f cq=%.3f cr=%.3f | R2 val p/q/r = %.3f/%.3f/%.3f | acc %.3f/%.3f/%.3f | %.0f s\\n', ...\n"
    "    DIAG.tag, P_final(13), P_final(14), P_final(15), summary.R2_final, summary.R2acc_final, summary.elapsed_s);\n"
    "disp('Script finalizado.');\n")

rep("function [aX, aY, aZ] = accel_eval(",
    "function se = diag_se_at(cost_fun, Pv, lb, ub, opts, n_params)\n"
    "%DIAG_SE_AT  Erro-padrão (≈ Cramér-Rao) de um custo avaliado em Pv, sem otimizar.\n"
    "    se = nan(n_params,1);\n"
    "    try\n"
    "        [~, rn, res, ~, ~, ~, J] = lsqnonlin(cost_fun, Pv(:), lb, ub, opts);\n"
    "        dof = max(numel(res) - n_params, 1);\n"
    "        se  = sqrt(max(diag((rn/dof) * pinv(full(J'*J))), 0));\n"
    "    catch ME\n"
    "        fprintf('  diag_se_at: falhou (%s)\\n', ME.message);\n"
    "    end\n"
    "end\n\n"
    "function r2 = iff_r2(res, kind, pqr_vl, acc_vl, R2f)\n"
    "%IFF_R2  R² de validação [1x3] ou NaN se a simulação divergiu (DIAG)\n"
    "    r2 = nan(1,3);\n"
    "    if strcmp(kind,'pqr') && res.pqr_vl_ok\n"
    "        r2 = [R2f(pqr_vl(:,1),res.p_s_vl), R2f(pqr_vl(:,2),res.q_s_vl), R2f(pqr_vl(:,3),res.r_s_vl)];\n"
    "    elseif strcmp(kind,'acc') && isfield(res,'full_vl_ok') && res.full_vl_ok\n"
    "        r2 = [R2f(acc_vl(:,1),res.accX_s_vl), R2f(acc_vl(:,2),res.accY_s_vl), R2f(acc_vl(:,3),res.accZ_s_vl)];\n"
    "    end\n"
    "end\n\n"
    "function [aX, aY, aZ] = accel_eval(")

open(DST, "w", encoding="utf-8").write(s)
print("ok:", DST, len(s.splitlines()), "linhas")
