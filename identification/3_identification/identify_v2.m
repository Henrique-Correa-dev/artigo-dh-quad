% identify_v2.m — PIPELINE DE IDENTIFICAÇÃO ORGANIZADO (VTOL quad)
% =========================================================================
% Processo limpo e bem-posto, com a filosofia "mede o que dá pra medir,
% identifica só o que o voo observa".
%
%   FIXO (medido / CAD):   massa (balança), inércias (CAD), braços, r_imu,
%                          curva fT/fQ (bancada)
%   IDENTIFICA (voo):      k_T(4), k_Q(4), Dp, Dq, Dr   — motor + amortecimento
%   NÃO identifica:        atitude (integrador) → valida em pqr; acc = diagnóstico
%
% POR QUE FIXAR A INÉRCIA: em roll, ṗ = (1/Jx)·braço·k_T·ΔfT − Dp·p, então
%   k_T e Jx são DEGENERADOS se ambos livres (vimos std=337% no teste). Fixar
%   J no CAD torna k_T identificável, elimina a penalidade triangular e o
%   problema de observabilidade do yaw (Jz). Processo fica bem-posto.
%
% MÉTODO (hierarquia clássica — Klein & Morelli / Jategaonkar):
%   1. Pré-proc: filtro Savitzky-Golay CONSISTENTE (pqr, acc, pwm) + derivada SG
%   2. EEM (equation-error): regressão robusta nas taxas → estimativa inicial
%   3. OEM (output-error) na malha de TAXA (estável, sem deriva) → refino final
%   4. Validação (janela held-out rica): R², RMSE, TIC + decomposição de Theil
%   5. Qualidade: erro-padrão de parâmetro (≈ Cramér-Rao) só nos parâmetros livres
%
% Uso:  >> identify_v2

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ===================== 0. CONFIG =====================
LOG_FILE = 'logs_concat.mat';
t_trains = {[4, 24]; [42, 62]; [100, 125]};   % janelas de treino (base concat)
t_val    = [605, 625];                          % validação (rica em p,q,r,acc)
dt       = 0.1;                                 % passo da grade comum (10 Hz)
% O PWM (RCOU) é logado a ~10 Hz NA FONTE (ArduPilot) → esta é a taxa do canal
% de entrada. Grids mais finos (ex.: 0.02) só INTERPOLAM o PWM: não trazem info
% e deixam o CRB otimista (resíduos correlacionados). Testado: TIC ~igual.
% Para subir a banda de verdade → logar o RCOU mais rápido no próximo voo.

% Savitzky-Golay (suaviza + deriva de forma consistente). frame ímpar > order.
SG_ORDER = 2;
SG_FRAME = 7;

% Atraso modelo×medido (lag_diagnostic deu ~1 amostra). DOIS modelos pra ele,
% que DISTINGUEM a causa (rode os dois e compare TIC, overshoot e brancura):
%   LAG_TAU > 0  → lag de 1ª ordem no PWM (passa-baixa): DESLOCA + ARREDONDA
%                  os picos. Se isso casar melhor (mata overshoot) → LAG DE MOTOR.
%   DELAY_PWM>0  → atraso PURO (shift de N amostras): só desloca, não arredonda.
%                  Se o shift puro já bastar → mais SYNC (artefato de timestamp).
% Use UM por vez (LAG_TAU tem prioridade). 0 = desliga.
% >>> lag_sweep CONFIRMOU erro de fase em p,q: TIC mínimo em τ≈150-200 ms
%     (p 0.55→0.43, q 0.46→0.32) e PSD do resíduo com pico em 0.65 Hz. r NÃO
%     melhora com lag (é ganho/k_Q, 0.10 Hz). Re-identifica COM o lag ligado pra
%     o k_T/Dp re-ajustarem (espera-se Dp cair do 7.2 inflado).
DELAY_PWM = 0;         % amostras (atraso puro) — desligado, agora usamos LAG
LAG_TAU   = 0.15;      % s — lag de 1ª ordem (motor/rotor). τ* do sweep ≈ 0.15-0.20

% Regularização de pares de motor (yaw fraco → k_Q individual mal-posto)
REG.kt_pair = 100;
REG.kq_pair = 100;
REG.tri     = 0;       % inércia FIXA no CAD já é física → sem penalidade triangular

DO_OEM = true;         % refino output-error na malha de taxa

%% ===================== 1. CARREGAR + FILTRAR =====================
proj = parameters();
constants = struct('m', proj.m, 'g', proj.g);
[fT, fQ] = motor_models();

L = load_log_data(fullfile(paths.data, LOG_FILE));
% grade comum (SEM GPS — ele trunca o concat; não é usado aqui)
t0 = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t1 = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t0:dt:t1)';

ip = @(tt,xx) interp1(tt, xx, tg, 'linear');
raw.p = ip(L.time_IMU, L.gyrX_raw);  raw.q = ip(L.time_IMU, L.gyrY_raw);  raw.r = ip(L.time_IMU, L.gyrZ_raw);
raw.ax= ip(L.time_IMU, L.accX_raw);  raw.ay= ip(L.time_IMU, L.accY_raw);  raw.az= ip(L.time_IMU, L.accZ_raw);
raw.w1= ip(L.time_RCOU, L.pwm1_raw); raw.w2= ip(L.time_RCOU, L.pwm2_raw);
raw.w3= ip(L.time_RCOU, L.pwm3_raw); raw.w4= ip(L.time_RCOU, L.pwm4_raw);
raw.roll = ip(L.time_ATT, L.roll_deg);  raw.pitch = ip(L.time_ATT, L.pitch_deg);  raw.yaw = ip(L.time_ATT, L.yaw_deg);

% FILTRAGEM CONSISTENTE: o MESMO SG em entrada (pwm) e saída (pqr, acc).
% Filtro linear comuta com sistema linear → preserva a dinâmica I/O.
sg = @(x) sgolayfilt(x, SG_ORDER, SG_FRAME);
f.p = sg(raw.p); f.q = sg(raw.q); f.r = sg(raw.r);
f.ax= sg(raw.ax); f.ay= sg(raw.ay); f.az= sg(raw.az);
f.w1= sg(raw.w1); f.w2= sg(raw.w2); f.w3= sg(raw.w3); f.w4= sg(raw.w4);
% Atraso modelo×medido — lag de 1ª ordem (arredonda+atrasa) OU shift puro.
if LAG_TAU > 0
    a = exp(-dt/LAG_TAU);                      % passa-baixa discreto
    lp = @(x) filter(1-a, [1 -a], x);
    f.w1=lp(f.w1); f.w2=lp(f.w2); f.w3=lp(f.w3); f.w4=lp(f.w4);
    fprintf('  >>> PWM com LAG de 1ª ordem τ=%.3f s (arredonda picos)\n', LAG_TAU);
elseif DELAY_PWM > 0
    dl = @(x) [repmat(x(1),DELAY_PWM,1); x(1:end-DELAY_PWM)];
    f.w1=dl(f.w1); f.w2=dl(f.w2); f.w3=dl(f.w3); f.w4=dl(f.w4);
    fprintf('  >>> PWM com atraso PURO de %d amostra(s) (só desloca)\n', DELAY_PWM);
end
% derivadas das taxas (para o EEM) — gradiente sobre o sinal já suavizado
f.pd = gradient(f.p, dt);  f.qd = gradient(f.q, dt);  f.rd = gradient(f.r, dt);
% atitude: CRUA (não filtra — yaw faz wrap em ±180° e o SG corromperia o salto).
% É só a referência medida usada na validação 'hybrid'.
f.roll = raw.roll;  f.pitch = raw.pitch;  f.yaw = raw.yaw;

fprintf('=== identify_v2 | log "%s" [%.1f, %.1f]s | SG(ord=%d,frame=%d) ===\n', ...
    LOG_FILE, t0, t1, SG_ORDER, SG_FRAME);

%% ===================== 2. SEGMENTOS DE TREINO =====================
seg = build_segs(t_trains, tg, f, fT, fQ);
fprintf('  %d segmentos de treino | validação [%g,%g]s\n', numel(seg), t_val(1), t_val(2));

%% ===================== 3. PARÂMETROS: P0 + bounds (inércia FIXA) =====================
P0 = proj.P0_J;             % [J(1:4); k_T(5:8); k_Q(9:12); D(13:15)]
lb = proj.bounds.lb;  ub = proj.bounds.ub;
% FIXA inércia no CAD (lb==ub ⇒ lsqnonlin não otimiza)
lb(1:4) = P0(1:4);  ub(1:4) = P0(1:4);
free = (ub > lb);           % máscara de parâmetros livres (k_T,k_Q,D)
fprintf('  Inércia FIXA no CAD: Jx=%.4f Jy=%.4f Jz=%.4f Jxz=%.5f\n', P0(1:4));
fprintf('  Parâmetros livres: %d (k_T,k_Q,Dp,Dq,Dr)\n', nnz(free));

%% ===================== 4. EEM (equation-error) =====================
% Concatena derivadas de todos os segmentos; pesos = 1/var de cada derivada.
P_all=[]; Q_all=[]; R_all=[]; Pd=[]; Qd=[]; Rd=[]; Tr=[]; Qr=[];
for s=1:numel(seg)
    sgs=seg{s};
    P_all=[P_all;sgs.p]; Q_all=[Q_all;sgs.q]; R_all=[R_all;sgs.r];
    Pd=[Pd;sgs.pd]; Qd=[Qd;sgs.qd]; Rd=[Rd;sgs.rd];
    Tr=[Tr;sgs.T_ref]; Qr=[Qr;sgs.Q_ref];
end
w_eem = [1/max(var(Pd),1e-9); 1/max(var(Qd),1e-9); 1/max(var(Rd),1e-9)];
cost_eem = @(P) eem_cost_function(ones(15,1), P, w_eem, P_all,Q_all,R_all, Pd,Qd,Rd, Tr,Qr, REG);
opt = optimoptions('lsqnonlin','Display','iter','MaxIterations',1000, ...
    'MaxFunctionEvaluations',60000,'FunctionTolerance',1e-12,'StepTolerance',1e-12);
fprintf('\n--- EEM ---\n');
P_eem = lsqnonlin(cost_eem, P0, lb, ub, opt);

%% ===================== 5. OEM (output-error, malha de TAXA) =====================
P_final = P_eem;
if DO_OEM
    fprintf('\n--- OEM (malha de taxa, 3 estados, RK4 passo-fixo) ---\n');
    % Pesos OEM = 1/var(pqr) (resíduo é em pqr, não na derivada)
    w_oem = [1/max(var(P_all),1e-9); 1/max(var(Q_all),1e-9); 1/max(var(R_all),1e-9)];
    cost_oem = @(P) oem_rate_cost(P, seg, fT, fQ, constants, w_oem, REG);
    % Opções próprias: RK4 fixo dá gradiente liso → poucas iterações bastam.
    opt_oem = optimoptions('lsqnonlin','Display','iter','MaxIterations',80, ...
        'MaxFunctionEvaluations',20000,'FunctionTolerance',1e-8,'StepTolerance',1e-10);
    P_final = lsqnonlin(cost_oem, P_eem, lb, ub, opt_oem);
end

%% ===================== 6. VALIDAÇÃO (janela held-out) =====================
v = window_data(t_val, tg, f);
res = sim_window('hybrid', P_final, v.time, v.pwm, v.pqr, v.att, constants);

fprintf('\n==========================================================\n');
fprintf('  VALIDAÇÃO [%g-%gs] — modo hybrid (pqr integrado + att medida)\n', t_val(1), t_val(2));
fprintf('==========================================================\n');
chans = {'p',v.pqr(:,1),res.p,'rad/s'; 'q',v.pqr(:,2),res.q,'rad/s'; 'r',v.pqr(:,3),res.r,'rad/s'; ...
         'accX',v.acc(:,1),res.accX,'m/s²'; 'accY',v.acc(:,2),res.accY,'m/s²'; 'accZ',v.acc(:,3),res.accZ,'m/s²'};
fprintf('  %-5s %8s %9s %7s %7s  Theil bias/var/cov\n','canal','R²','RMSE','TIC','');
for c=1:size(chans,1)
    fm = fit_metrics(chans{c,2}, chans{c,3});
    flag = ''; if fm.TIC>0.30, flag='  TIC>0.30'; end
    fprintf('  %-5s %8.3f %9.4f %7.3f %7s  %.2f/%.2f/%.2f%s\n', ...
        chans{c,1}, fm.R2, fm.RMSE, fm.TIC, '', fm.U_bias, fm.U_var, fm.U_cov, flag);
end
fprintf('  (taxas: TIC + cov→1 = resíduo branco = modelo certo; acc = diagnóstico)\n');

%% ===================== 6b. TESTE DE BRANCURA DO RESÍDUO =====================
% Resíduo BRANCO (ACF dentro das bandas 95%, exceto lag 0) ⇒ erro irredutível
% (ruído/banda) → modelo no limite do que os dados permitem.
% Resíduo COLORIDO (ACF estoura em lags>0) ⇒ dinâmica não modelada (melhorável).
fprintf('\n--- Teste de brancura do resíduo (taxas) ---\n');
maxlag = 30;
ewn = {'p', v.pqr(:,1)-res.p; 'q', v.pqr(:,2)-res.q; 'r', v.pqr(:,3)-res.r};
figW = figure('Name','identify_v2 brancura','Position',[80 80 1100 700]);
for c = 1:3
    e = ewn{c,2};
    [acf, lags] = acf_(e, maxlag);
    N = numel(e);  b = 1.96/sqrt(N);
    inside = mean(abs(acf(2:end)) <= b) * 100;     % % de lags (k≥1) dentro da banda
    if inside >= 90, verdict = '~BRANCO (irredutível)'; else, verdict = 'COLORIDO (sobra dinâmica)'; end
    fprintf('  %-3s: %3.0f%% dos lags dentro de ±%.3f (95%%)  → %s\n', ewn{c,1}, inside, b, verdict);
    subplot(3,1,c); hold on; grid on;
    stem(lags(2:end), acf(2:end), 'filled', 'MarkerSize', 3);
    yline(b,'r--'); yline(-b,'r--');
    ylabel(['ACF ' ewn{c,1}]); ylim([-0.5 0.5]);
    if c==1, title(sprintf('Brancura do resíduo [%g-%gs] — dentro das bandas = branco', t_val(1), t_val(2))); end
    if c==3, xlabel('lag (amostras)'); end
end
saveas(figW, fullfile(paths.images,'identify_v2_whiteness.png'));
fprintf('  (≥90%% dentro da banda ≈ branco. Figura: %s)\n', fullfile(paths.images,'identify_v2_whiteness.png'));

%% ===================== 7. QUALIDADE DOS PARÂMETROS (≈ CRB) =====================
fprintf('\n--- Qualidade dos parâmetros livres (erro-padrão ≈ Cramér-Rao) ---\n');
try
    [~, rn, rsd, ~, ~, ~, Jc] = lsqnonlin(cost_eem, P_final, lb, ub, ...
        optimoptions('lsqnonlin','Display','off','MaxIterations',0));
    idxf = find(free);
    Jf = full(Jc(:, idxf));                       % só colunas dos livres
    dof = max(numel(rsd) - numel(idxf), 1);
    Cov = (rn/dof) * pinv(Jf'*Jf);
    se  = sqrt(max(diag(Cov),0));
    fprintf('  %-6s %12s %12s %9s\n','param','valor','std','std%%');
    for i=1:numel(idxf)
        k=idxf(i); rel=100*se(i)/max(abs(P_final(k)),1e-9);
        tag=''; if rel>50, tag='  << fraco'; end
        fprintf('  %-6s %12.6f %12.6f %8.1f%%%s\n', proj.param_names{k}, P_final(k), se(i), rel, tag);
    end
catch ME
    fprintf('  (não calculado: %s)\n', ME.message);
end

%% ===================== 8. SALVAR + PLOT =====================
save(fullfile(paths.outputs,'P_identified_v2.mat'), 'P_final', 'P_eem', 'P0', 't_val', 't_trains');
fprintf('\n  Salvo: %s (P_identified_v2.mat)\n', paths.outputs);

fig=figure('Name','identify_v2 validação','Position',[60 60 1200 700]);
lab={'p (rad/s)','q (rad/s)','r (rad/s)'}; simv={res.p,res.q,res.r};
for i=1:3
    subplot(3,1,i); hold on; grid on;
    plot(v.time, v.pqr(:,i),'b-','LineWidth',1.3,'DisplayName','medido');
    plot(v.time, simv{i},'r--','LineWidth',1.3,'DisplayName','modelo');
    ylabel(lab{i}); if i==1, legend('Location','best'); title(sprintf('Validação [%g-%gs]',t_val(1),t_val(2))); end
    if i==3, xlabel('t (s)'); end
end
saveas(fig, fullfile(paths.images,'identify_v2_pqr.png'));
fprintf('  Figura: %s\n', fullfile(paths.images,'identify_v2_pqr.png'));


%% ===================== FUNÇÕES LOCAIS =====================
function segs = build_segs(t_trains, tg, f, fT, fQ)
    segs = cell(numel(t_trains),1);
    for s=1:numel(t_trains)
        m = (tg>=t_trains{s}(1)) & (tg<=t_trains{s}(2));
        d.time = tg(m);
        d.p=f.p(m); d.q=f.q(m); d.r=f.r(m);
        d.pd=f.pd(m); d.qd=f.qd(m); d.rd=f.rd(m);
        d.pqr=[f.p(m),f.q(m),f.r(m)];
        d.pwm=[f.w1(m),f.w2(m),f.w3(m),f.w4(m)];
        d.T_ref=[fT(d.pwm(:,1)),fT(d.pwm(:,2)),fT(d.pwm(:,3)),fT(d.pwm(:,4))];
        d.Q_ref=[fQ(d.pwm(:,1)),fQ(d.pwm(:,2)),fQ(d.pwm(:,3)),fQ(d.pwm(:,4))];
        segs{s}=d;
    end
end

function [acf, lags] = acf_(e, maxlag)
% Autocorrelação normalizada (sem toolbox). acf(1)=lag 0=1.
    e = e(:) - mean(e);
    N = numel(e);  c0 = sum(e.^2) + 1e-12;
    acf = zeros(maxlag+1,1);
    for k = 0:maxlag
        acf(k+1) = sum(e(1:N-k) .* e(1+k:N)) / c0;
    end
    lags = (0:maxlag)';
end

function v = window_data(t_win, tg, f)
    m=(tg>=t_win(1))&(tg<=t_win(2));
    v.time=tg(m);
    v.pqr=[f.p(m),f.q(m),f.r(m)];
    v.acc=[f.ax(m),f.ay(m),f.az(m)];
    v.pwm=[f.w1(m),f.w2(m),f.w3(m),f.w4(m)];
    v.att=[f.roll(m),f.pitch(m),f.yaw(m)];
end

function e = oem_rate_cost(P, segs, fT, fQ, constants, w, REG)
% Output-error SÓ na malha de taxa (3 estados, estável → sem deriva).
% Integração RK4 de PASSO FIXO (determinística): pequenas perturbações de
% parâmetro mudam a saída de forma suave → Jacobiano por diferenças finitas
% limpo (ode45 de passo variável dá gradiente ruidoso e trava o lsqnonlin).
    sw = sqrt(w(:));  e = [];
    for s=1:numel(segs)
        sg = segs{s};
        Y  = rk4_rate(P, sg.time, sg.pwm, fT, fQ, constants, sg.pqr(1,:)');
        r  = Y - sg.pqr;
        e  = [e; sw(1)*r(:,1); sw(2)*r(:,2); sw(3)*r(:,3)];
    end
    % MESMA regularização de pares do EEM — senão o OEM desfaz o vínculo e,
    % com yaw fraco, os k_Q individuais vagam pra valores absurdos (k_Q→0).
    if nargin >= 7 && ~isempty(REG)
        kT = P(5:8);  kQ = P(9:12);
        e = [e; sqrt(REG.kt_pair)*(kT(1)-kT(2)); ...
                sqrt(REG.kt_pair)*(kT(3)-kT(4)); ...
                sqrt(REG.kq_pair)*(kQ(1)-kQ(2)); ...
                sqrt(REG.kq_pair)*(kQ(3)-kQ(4))];
    end
end

function Y = rk4_rate(P, t, pwm, fT, fQ, constants, y0)
% RK4 passo-fixo do modelo rotacional (3 estados) sobre a grade t.
    N = numel(t);  Y = zeros(N,3);  Y(1,:) = y0(:)';
    fr = @(tt,yy) vtol_dynamics(tt, yy, P, t, pwm, fT, fQ, constants);
    for k = 1:N-1
        h  = t(k+1) - t(k);
        yk = Y(k,:)';
        k1 = fr(t(k),       yk);
        k2 = fr(t(k)+h/2,   yk + h/2*k1);
        k3 = fr(t(k)+h/2,   yk + h/2*k2);
        k4 = fr(t(k)+h,     yk + h*k3);
        Y(k+1,:) = (yk + h/6*(k1 + 2*k2 + 2*k3 + k4))';
    end
end
