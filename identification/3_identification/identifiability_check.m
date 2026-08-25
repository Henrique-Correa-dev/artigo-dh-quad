% identifiability_check.m — o dado enxerga os coeficientes aerodinâmicos?
% =========================================================================
% Teste barato (sem otimização) que responde, ANTES de identificar, quais dos
% parâmetros do modelo estendido de 22 são separáveis no envelope voado.
%
% Método: Jacobiano do resíduo de EQUAÇÃO (EEM, sem integração) no chute
% inicial Θ₀, com a velocidade estimada de bordo (estimate_velocity). Do
% Jacobiano saem três leituras:
%   [1] sensibilidade escalada  |∂e/∂θ|·|θ| — quanto do resíduo o parâmetro move
%   [2] cota de Cramér-Rao      CR% = 100·σ_θ/|θ|  em Θ₀
%   [3] matriz de correlação    pares |ρ| > 0,9 são o que NÃO dá para separar
% Mais dois diagnósticos diretos:
%   [4] colinearidade crua dos regressores: corr(p, V·p), corr(q, V·q), corr(r, V·r)
%   [5] varredura de velocidade: CR% se o mesmo voo tivesse V×2, ×5, ×10
%       (os termos aerodinâmicos escalam com V, os de rotor não) → diz a que
%       velocidade a aerodinâmica da estrutura passaria a ser identificável.
%
% Uso:  >> identifiability_check
% Saída: tabela no console + outputs/images/identifiability_check.png
%        + outputs/identifiability_check.mat
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();

LOG_FILE = 'logs_concat.mat';
T_TRAINS = {[4, 24]; [25, 41]; [42, 62]; [63, 99]; [100, 125]};   % mesmas do identify_plant
V_SCALES = [1 2 5 10];      % multiplicadores hipotéticos de velocidade

%% ---------------------------------------------------------------- dados
L = load_log_data(fullfile(paths.data, LOG_FILE));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';  dt = 0.1;
ip = @(tt,xx) interp1(tt, xx, tg, 'linear');

W4 = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ...
      ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
[rpm, func_T_ref, func_Q_ref, mc] = motor_chain(tg, W4);
fprintf('  motor_chain: atraso %.2f s | tau_m %.3f s\n', mc.delay_s, mc.tau_m);

gyr = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
VEst = estimate_velocity(L, tg);
Vg = VEst.V;  vg = VEst.v;  wg = VEst.w;
Vg(~isfinite(Vg)) = 0;  vg(~isfinite(vg)) = 0;  wg(~isfinite(wg)) = 0;

smooth_win = 5;
p_all=[]; q_all=[]; r_all=[]; pd_all=[]; qd_all=[]; rd_all=[];
Tr_all=[]; Qr_all=[]; V_all=[]; v_all=[]; w_all=[];
for s = 1:numel(T_TRAINS)
    ii = tg >= T_TRAINS{s}(1) & tg <= T_TRAINS{s}(2);
    ps = gyr(ii,1); qs = gyr(ii,2); rs = gyr(ii,3);
    p_all=[p_all;ps]; q_all=[q_all;qs]; r_all=[r_all;rs];                        %#ok<AGROW>
    pd_all=[pd_all; gradient(movmean(ps,smooth_win),dt)];                        %#ok<AGROW>
    qd_all=[qd_all; gradient(movmean(qs,smooth_win),dt)];                        %#ok<AGROW>
    rd_all=[rd_all; gradient(movmean(rs,smooth_win),dt)];                        %#ok<AGROW>
    Tr_all=[Tr_all; func_T_ref(rpm(ii,:))]; Qr_all=[Qr_all; func_Q_ref(rpm(ii,:))]; %#ok<AGROW>
    V_all=[V_all; Vg(ii)]; v_all=[v_all; vg(ii)]; w_all=[w_all; wg(ii)];         %#ok<AGROW>
end
N = numel(p_all);
fprintf('  Treino: %d amostras (%.0f s).  V: média %.2f | p95 %.2f | máx %.2f m/s\n', ...
    N, N*dt, mean(V_all), prctile(V_all,95), max(V_all));

weights = [1/var(pd_all); 1/var(qd_all); 1/var(rd_all)];
RL = struct('kt_pair',0, 'kq_pair',0, 'tri',0);   % CRB do ajuste puro, sem regularização

%% ------------------------------------------------- conjunto de análise
P0 = proj.P0_J(:);  lb = proj.bounds.lb(:);  ub = proj.bounds.ub(:);
names = proj.param_names;  np = numel(P0);
locked = abs(ub - lb) < 1e-12;                    % Jz, Jxz travados no CAD
ana = true(np,1);  ana(16) = false;               % C_d não entra nos momentos
ana(locked) = false;
idx = find(ana);  nA = numel(idx);
fprintf('  Analisados %d parâmetros (fora: %s).\n', nA, strjoin(names(~ana), ', '));

kA = aero_gains(proj);

%% ------------------------------------------------------------ Jacobiano
res_of = @(P, Vs) trim7(eem_cost_function(ones(np,1), P, weights, ...
    p_all, q_all, r_all, pd_all, qd_all, rd_all, Tr_all, Qr_all, RL, ...
    struct('V', Vs*V_all, 'v', v_all, 'w', w_all)));

sc = max(abs(P0), 1e-3);
CR = nan(nA, numel(V_SCALES));  SENS = nan(nA, numel(V_SCALES));
for js = 1:numel(V_SCALES)
    Vs = V_SCALES(js);
    e0 = res_of(P0, Vs);  Nres = numel(e0);
    J = zeros(Nres, nA);
    for k = 1:nA
        i = idx(k);  h = 1e-5*sc(i);
        Pp = P0; Pp(i) = Pp(i) + h;   Pm = P0; Pm(i) = Pm(i) - h;
        J(:,k) = (res_of(Pp,Vs) - res_of(Pm,Vs)) / (2*h);
    end
    % σ² é propriedade do SENSOR, não da velocidade hipotética: fixa no caso
    % real (V×1). Sem isso a varredura mede a piora do ajuste em Θ₀ (que cresce
    % com V porque os termos ∝ V ficam grandes) em vez da informação de Fisher.
    H = J.'*J;
    if Vs == V_SCALES(1), sigma2 = sum(e0.^2)/max(Nres - nA, 1); end
    C = sigma2 * pinv(H);
    se = sqrt(max(diag(C), 0));
    CR(:,js)   = 100*se ./ max(abs(P0(idx)), eps);
    SENS(:,js) = vecnorm(J,2,1)' .* abs(P0(idx)) / sqrt(Nres);   % resíduo RMS por 100% de θ
    if Vs == 1
        R1 = C ./ max(se*se.', eps);  H1 = H;  J1 = J;  cond1 = cond(H);
    end
end

%% ------------------------------------------------------------- relatório
fprintf('\n  ==========================================================================\n');
fprintf('   IDENTIFICABILIDADE EM Θ₀ (EEM, %d amostras de treino)   cond(JᵀJ) = %.1e\n', N, cond1);
fprintf('  ==========================================================================\n');
fprintf('  %-6s %10s %10s | %s\n', 'param', 'Θ₀', 'sens.', 'CR% com V×1   V×2   V×5   V×10');
for k = 1:nA
    fprintf('  %-6s %10.4f %10.2e |  %8.1f %8.1f %8.1f %8.1f\n', ...
        names{idx(k)}, P0(idx(k)), SENS(k,1), CR(k,1), CR(k,2), CR(k,3), CR(k,4));
end
fprintf('  (sens. = variação RMS do resíduo normalizado por 100%% de variação do parâmetro)\n');

fprintf('\n  Pares com |correlação| > 0,90 (não separáveis):\n');
found = false;
for a = 1:nA-1
    for b = a+1:nA
        if abs(R1(a,b)) > 0.90
            fprintf('    %-6s ↔ %-6s   ρ = %+.3f\n', names{idx(a)}, names{idx(b)}, R1(a,b));
            found = true;
        end
    end
end
if ~found, fprintf('    (nenhum)\n'); end

pos = @(nm) find(strcmp(names(idx), nm));
key = {'Dp','Cl_p'; 'Dq','Cm_q'; 'Dr','Cn_r'; 'Cl_b','Cn_b'; 'Dq','Cm_a'; 'Dp','Cl_b'};
fprintf('\n  Correlação dos pares rotor ↔ estrutura que interessam:\n');
for a = 1:size(key,1)
    ia = pos(key{a,1});  ib = pos(key{a,2});
    if ~isempty(ia) && ~isempty(ib)
        fprintf('    %-6s ↔ %-6s   ρ = %+.3f\n', key{a,1}, key{a,2}, R1(ia,ib));
    end
end

%% --------------------------------- colinearidade crua dos regressores
cc = @(a,b) corr_(a, b);
fprintf('\n  Colinearidade dos regressores de amortecimento (rotor × estrutura):\n');
fprintf('    corr(p, V·p) = %+.4f   corr(q, V·q) = %+.4f   corr(r, V·r) = %+.4f\n', ...
    cc(p_all, V_all.*p_all), cc(q_all, V_all.*q_all), cc(r_all, V_all.*r_all));

Jx = P0(1); Jy = P0(2); Jz = P0(3);
mag = @(x) std(x);
fprintf('\n  Peso relativo dos termos (desvio-padrão de cada contribuição a ṗ, q̇, ṙ):\n');
fprintf('    ṗ: rotor c_p·p = %.4f rad/s²   estrutura Cl_p = %.4f rad/s²  (%.1f%%)\n', ...
    mag(P0(13)*p_all), mag(kA.Lp*P0(17)*V_all.*p_all/Jx), 100*mag(kA.Lp*P0(17)*V_all.*p_all/Jx)/mag(P0(13)*p_all));
fprintf('    q̇: rotor c_q·q = %.4f rad/s²   estrutura Cm_q = %.4f rad/s²  (%.1f%%)\n', ...
    mag(P0(14)*q_all), mag(kA.Mq*P0(19)*V_all.*q_all/Jy), 100*mag(kA.Mq*P0(19)*V_all.*q_all/Jy)/mag(P0(14)*q_all));
fprintf('    ṙ: rotor c_r·r = %.4f rad/s²   estrutura Cn_r = %.4f rad/s²  (%.1f%%)\n', ...
    mag(P0(15)*r_all), mag(kA.Nr*P0(21)*V_all.*r_all/Jz), 100*mag(kA.Nr*P0(21)*V_all.*r_all/Jz)/mag(P0(15)*r_all));
fprintf('    (ruído do resíduo em Θ₀: RMS ṗ %.3f, q̇ %.3f, ṙ %.3f rad/s²)\n', ...
    rms_(pd_all), rms_(qd_all), rms_(rd_all));

save(fullfile(paths.outputs,'identifiability_check.mat'), 'CR','SENS','R1','idx','names','V_SCALES','P0','V_all','N');

%% ----------------------------------------------------------------- figura
f = figure('Position',[40 40 1300 820],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(2, 2, 'TileSpacing','compact','Padding','compact');

nexttile;                                    % matriz de correlação
imagesc(abs(R1), [0 1]); axis square; colormap(gca, flipud(gray)); colorbar;
set(gca,'XTick',1:nA,'XTickLabel',names(idx),'YTick',1:nA,'YTickLabel',names(idx),'FontSize',8);
xtickangle(90);
hold on; na_rot = sum(idx <= 15);
plot([na_rot+0.5 na_rot+0.5],[0.5 nA+0.5],'r-','LineWidth',1.5);
plot([0.5 nA+0.5],[na_rot+0.5 na_rot+0.5],'r-','LineWidth',1.5);
text(0.0, 1.06, 'Correlação absoluta dos estimadores em \Theta_0 (linha vermelha separa rotor de aerodinâmica)', ...
    'Units','normalized','FontWeight','bold','FontSize',9);

nexttile;                                    % CR% por velocidade
hb = bar(categorical(names(idx), names(idx)), CR, 'grouped');
set(gca,'YScale','log'); grid on; ylabel('CR% em \Theta_0 (log)');
yline(100,'r--','100% = sem informação','LabelHorizontalAlignment','left');
yline(20,'k--','20% = aceitável','LabelHorizontalAlignment','left');
legend(arrayfun(@(s) sprintf('V \\times %d', s), V_SCALES, 'UniformOutput',false), ...
    'Location','northwest','NumColumns',4);
text(0.0, 1.06, 'Cota de Cramér-Rao: efeito de voar mais rápido', 'Units','normalized','FontWeight','bold','FontSize',9);

nexttile;                                    % séries: p e os dois amortecimentos
tt = (0:N-1)'*dt;  hold on; grid on;
plot(tt, P0(13)*p_all, '-', 'Color',[0 0.45 0.7], 'LineWidth',1.1);
plot(tt, kA.Lp*P0(17)*V_all.*p_all/Jx, '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.1);
plot(tt, pd_all, 'k-', 'LineWidth',0.4);
ylabel('[rad/s²]'); xlabel('t concatenado de treino [s]'); xlim([0 tt(end)]);
legend({'rotor  c_p·p (\Theta_0)','estrutura  ¼\rhoSb²C_{lp}V·p/J_x (AVL)','ṗ medido'}, ...
    'Location','northoutside','Orientation','horizontal','FontSize',8);

nexttile;                                    % velocidade estimada no treino
hold on; grid on;
plot(tt, V_all, 'k-', 'LineWidth',1);
yline(prctile(V_all,95), '--', sprintf('p95 = %.2f m/s', prctile(V_all,95)), 'Color',[0.85 0.37 0.01]);
ylabel('V estimada [m/s]'); xlabel('t concatenado de treino [s]'); xlim([0 tt(end)]);
text(0.0, 1.06, 'velocidade que escala toda a aerodinâmica da estrutura', 'Units','normalized','FontWeight','bold','FontSize',9);

fn = fullfile(paths.images,'identifiability_check.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',140);
fprintf('\n  Figura: %s\n', fn);

%% ------------------------------------------------------------- utilidades
function e = trim7(e), e = e(1:end-7); end     % remove os 7 resíduos de regularização
function c = corr_(a, b)
    a = a(:) - mean(a(:));  b = b(:) - mean(b(:));
    c = (a.'*b) / max(sqrt((a.'*a)*(b.'*b)), eps);
end
function y = rms_(x), y = sqrt(mean(x(:).^2)); end
