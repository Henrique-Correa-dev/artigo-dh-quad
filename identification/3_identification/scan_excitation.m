% scan_excitation.m — Diagnóstico de EXCITAÇÃO para justificar a ENTRADA do OEM
% =========================================================================
% Visão geral dos quatro sinais da cadeia de identificação e quantificação de
% se as manobras escolhidas atendem aos critérios de uma boa entrada:
%
%   RCIN (piloto, sigma_P)  ->  RCOU/PWM (controlador->motores)  ->  MOMENTOS tau
%   (entrada exogena)            (saida do mixer)                    (entrada REAL do OEM)
%                                                                  ->  pqr (resposta)
%
% Os MOMENTOS sao calculados com a MESMA alocacao do estimador
% (vtol_dynamics -> forces_from_pwm_local), entao a manobra e avaliada na
% EXATA entrada que o OEM enxerga -- nao no PWM bruto nem no comando do piloto.
%
% Analises (Klein & Morelli / Jategaonkar):
%   1. Visao geral no tempo:  RCIN | PWM | momentos | pqr  (segmentos por eixo)
%   2. SNR e ISOLAMENTO (rho) por janela, por eixo
%   3. Conteudo de frequencia (PSD) do momento por eixo  (energia na banda, ~0 em DC)
%   4. Relatorio por segmento (yaw/roll/pitch)
%
% Nao modifica nada -- e so diagnostico.
% Uso:  >> scan_excitation

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---------------- CONFIG ----------------
LOG_FILE = 'logs_concat.mat';
dt   = 0.1;                 % 10 Hz (taxa do RCOU na fonte)
WIN  = 20;  STEP = 10;      % varredura: janela 20 s, passo 10 s (50% overlap)
SNR_MIN = 3;               % razao minima resposta/ruido (Jategaonkar)

% Janelas REAIS usadas no OEM (identify_v2.m) + validacao — referencia (nao plotadas)
% USED_TRAIN = [4 24; 42 62; 100 125];
% USED_VAL   = [605 625];

% Segmentos de manobra ISOLADA por eixo: {nome, [t0 t1], eixo}  (1=roll/p, 2=pitch/q, 3=yaw/r)
AXSEG = { 'Guinada (yaw)',   [  0  40], 3;
          'Rolagem (roll)',  [ 40  80], 1;
          'Arfagem (pitch)', [ 80 125], 2 };

axname = {'roll (p)','pitch (q)','yaw (r)'};
taulbl = {'\tau_x','\tau_y','\tau_z'};

%% ---------------- MODELO DE MOMENTO (idêntico ao identify_plant) ----------------
% Mesmo pré-processo do OEM: atraso de 1 amostra + lag de 1a ordem (tau_m) no PWM
% e modelo de motor FISICO (PWM->RPM->k*RPM^2). O PWM cru segue no painel de RCOU;
% o momento usa o PWM atrasado+lag (== o que o estimador enxerga).
DELAY_PWM = 1;                             % amostras (atraso puro, = identify_plant)
proj = parameters();
tau_lag = proj.motor.tau_m;               % s — lag 1a ordem (spin-up do motor)
[fT, fQ] = motor_models('physical');      % T,Q = k*(CR*PWM+w_b)^2  (= identify_plant)
dyn = vtol_dynamics('get_handles');        % dyn.forces: PWM -> [T, tx, ty, tz]
lmf = fullfile(paths.outputs, 'linear_model.mat');
if exist(lmf, 'file')
    S = load(lmf);  P = S.P;  Psrc = 'linear_model.mat (identificado)';
else
    P = proj.P0_J;            Psrc = 'parameters().P0_J (chute inicial)';
end

%% ---------------- CARREGAR + GRADE ----------------
L  = load_log_data(fullfile(paths.data, LOG_FILE));
t0 = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t1 = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t0:dt:t1)';
ip = @(tt, xx) interp1(tt, xx, tg, 'linear');

% RESPOSTA: taxas angulares (gyro)
p = ip(L.time_IMU, L.gyrX_raw);
q = ip(L.time_IMU, L.gyrY_raw);
r = ip(L.time_IMU, L.gyrZ_raw);
PQR = [p q r];

% ENTRADA NO SISTEMA: PWM dos 4 motores (RCOU = saida do controlador)
W = [ip(L.time_RCOU, L.pwm1_raw), ip(L.time_RCOU, L.pwm2_raw), ...
     ip(L.time_RCOU, L.pwm3_raw), ip(L.time_RCOU, L.pwm4_raw)];

% ENTRADA REAL DO OEM: PWM com atraso+lag (== identify_plant), depois alocacao
dl    = @(x) [repmat(x(1),DELAY_PWM,1); x(1:end-DELAY_PWM)];
a_lag = exp(-dt/tau_lag);
lp    = @(x) filter(1-a_lag, [1 -a_lag], x, a_lag*x(1));
W_in  = W;  for j = 1:4, W_in(:,j) = lp(dl(W(:,j))); end   % PWM cru fica em W (painel RCOU)
[Ttot, tx, ty, tz] = dyn.forces(W_in, P, fT, fQ);
TAU = [tx ty tz];

% COMANDO DO PILOTO (exogeno): RCIN
hasRC = isfield(L, 'time_RCIN') && ~isempty(L.time_RCIN);
if hasRC
    ipr = @(xx) interp1(L.time_RCIN, xx, tg, 'linear', 'extrap');
    RC = [ipr(L.rcin_roll), ipr(L.rcin_pitch), ipr(L.rcin_yaw)];
end

% GAPS entre logs concatenados (salto > 1 s no tempo bruto do IMU)
gap = false(size(tg));
for j = find(diff(L.time_IMU) > 1.0)'
    gap = gap | (tg > L.time_IMU(j) & tg < L.time_IMU(j+1));
end

%% ---------------- RUÍDO (robusto, por eixo) ----------------
% ruido de medicao ~ residuo de alta-frequencia (raw - Savitzky-Golay), sigma via MAD
hf   = @(x) x - sgolayfilt(x, 2, 7);
rstd = @(x) 1.4826 * median(abs(hf(x) - median(hf(x))));   % sigma robusto
noise = [rstd(p), rstd(q), rstd(r)];                       % rad/s

%% ---------------- VARREDURA POR JANELA ----------------
starts = (t0 : STEP : t1 - WIN)';
nW  = numel(starts);
Wc  = starts + WIN/2;
SNRw = nan(nW,3);  RHOw = nan(nW,3);  flyw = false(nW,1);  valw = false(nW,1);
for k = 1:nW
    m = tg >= starts(k) & tg < starts(k)+WIN;
    valw(k) = ~any(gap(m));
    flyw(k) = mean(mean(W(m,:), 2)) > 1400;          % acima de idle -> voando
    for ax = 1:3
        oth = setdiff(1:3, ax);
        SNRw(k,ax) = std(PQR(m,ax)) / noise(ax);
        RHOw(k,ax) = std(PQR(m,ax)) / sqrt(std(PQR(m,oth(1)))^2 + std(PQR(m,oth(2)))^2);
    end
end
ok = valw & flyw;

%% ---------------- RELATÓRIO ----------------
fprintf('\n==========================================================\n');
fprintf('  EXCITACAO -- %s  | P: %s\n', LOG_FILE, Psrc);
fprintf('  ruido gyro (sigma):  p=%.4f  q=%.4f  r=%.4f rad/s\n', noise);
fprintf('==========================================================\n');
fprintf('  SEGMENTOS DE MANOBRA ISOLADA (eixo-alvo):\n\n');
fprintf('  %-16s %-10s %8s %8s %9s %7s %9s\n', ...
        'segmento','alvo','std_tau','std_taxa','isol_rho','SNR','fpico(Hz)');
for s = 1:size(AXSEG,1)
    nm = AXSEG{s,1};  seg = AXSEG{s,2};  ax = AXSEG{s,3};
    m  = tg >= seg(1) & tg < seg(2) & ~gap;
    oth = setdiff(1:3, ax);
    s_tau  = std(TAU(m,ax));
    s_rate = std(PQR(m,ax));
    rho    = s_rate / sqrt(std(PQR(m,oth(1)))^2 + std(PQR(m,oth(2)))^2);
    snr    = s_rate / noise(ax);
    fpk    = peak_freq(TAU(m,ax), 1/dt);
    fprintf('  %-16s %-10s %8.4f %8.4f %9.2f %7.1f %9.2f%s\n', ...
            nm, axname{ax}, s_tau, s_rate, rho, snr, fpk, okmark(snr >= SNR_MIN));
end
fprintf('\n  Criterios: isol_rho >> 1 (eixo isolado) ; SNR >= %g (domina o ruido).\n', SNR_MIN);
fprintf('  Obs.: guinada (r) e o canal fraco (tau_z = so torque de reacao) -> menor SNR.\n');

%% ---------------- FIG 1: visão geral (RCIN | PWM | momentos | pqr) ----------------
fig1 = figure('Name','scan_excitation — visao geral','Position',[40 40 1320 920]);
set(fig1,'Color','w'); try, fig1.Theme = 'light'; catch, end
ax1 = subplot(4,1,1); hold on; grid on;
if hasRC
    plot(tg,RC(:,1),'b'); plot(tg,RC(:,2),'r'); plot(tg,RC(:,3),'Color',[0 .6 0]);
    shade_axsegs(AXSEG); shade_gaps(tg,gap);
    legend('roll','pitch','yaw','Location','eastoutside'); ylabel('RCIN (\mus)');
else
    text(0.5,0.5,'RCIN ausente no log','Units','normalized','HorizontalAlignment','center');
end
title('Comando do piloto — RCIN (entrada exogena \sigma_P)');
ax2 = subplot(4,1,2); hold on; grid on;
plot(tg,W(:,1)); plot(tg,W(:,2)); plot(tg,W(:,3)); plot(tg,W(:,4));
shade_axsegs(AXSEG); shade_gaps(tg,gap);
legend('M1','M2','M3','M4','Location','eastoutside'); ylabel('PWM (\mus)');
title('Saida do controlador — RCOU/PWM dos motores (entrada no sistema)');
ax3 = subplot(4,1,3); hold on; grid on;
plot(tg,tx,'b'); plot(tg,ty,'r'); plot(tg,tz,'Color',[0 .6 0]);
shade_axsegs(AXSEG); shade_gaps(tg,gap);
legend('\tau_x','\tau_y','\tau_z','Location','eastoutside'); ylabel('momento (N\cdotm)');
title('Momentos gerados (entrada REAL do OEM — alocacao identica ao vtol\_dynamics)');
ax4 = subplot(4,1,4); hold on; grid on;
plot(tg,p,'b'); plot(tg,q,'r'); plot(tg,r,'Color',[0 .6 0]);
shade_axsegs(AXSEG); shade_gaps(tg,gap);
legend('p','q','r','Location','eastoutside'); ylabel('taxas (rad/s)'); xlabel('t (s)');
title('Resposta — taxas angulares (saida)');
linkaxes([ax1 ax2 ax3 ax4],'x'); xlim([t0 t1]);
exportgraphics(fig1, fullfile(paths.images,'scan_excitation_overview.png'), 'BackgroundColor','white');
% versao recortada ate o fim do ultimo segmento (125 s) — remove o pouso posterior
maxseg = max(cellfun(@(s) s(2), AXSEG(:,2)));
xlim(ax4, [0 maxseg]);                                     % linkaxes propaga aos 4 paineis
exportgraphics(fig1, fullfile(paths.images,'scan_excitation_overview_tese.png'), 'BackgroundColor','white');
xlim(ax4, [t0 t1]);                                        % restaura visao completa

%% ---------------- FIG 2: SNR e isolamento por janela, por eixo ----------------
fig2 = figure('Name','scan_excitation — SNR e isolamento','Position',[60 60 1320 820]);
set(fig2,'Color','w'); try, fig2.Theme = 'light'; catch, end
order = [1 2 3];     % roll, pitch, yaw
for ii = 1:3
    ax = order(ii);
    subplot(3,1,ii); hold on; grid on;
    yyaxis left
    stem(Wc(ok), SNRw(ok,ax), 'filled', 'MarkerSize',4);
    yline(SNR_MIN,'--','SNR mín');
    ylabel(sprintf('SNR (%s)', axname{ax}));
    yyaxis right
    stem(Wc(ok), RHOw(ok,ax), 'MarkerSize',4);
    ylabel('isolamento \rho');
    seg = AXSEG{cellfun(@(a) a==ax, AXSEG(:,3)),2};
    yl = ylim; patch([seg(1) seg(2) seg(2) seg(1)],[yl(1) yl(1) yl(2) yl(2)], ...
        [.7 .9 .7],'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off'); ylim(yl);
    title(sprintf('Eixo %s — SNR e isolamento por janela (faixa verde = manobra deste eixo)', axname{ax}));
    if ii==3, xlabel('t (s)'); end
    xlim([t0 t1]);
end
exportgraphics(fig2, fullfile(paths.images,'scan_excitation_snr.png'), 'BackgroundColor','white');

%% ---------------- FIG 3: conteúdo de frequência (PSD do momento por eixo) ----------------
fig3 = figure('Name','scan_excitation — frequencia','Position',[80 80 1100 820]);
set(fig3,'Color','w'); try, fig3.Theme = 'light'; catch, end
fs = 1/dt;
for s = 1:size(AXSEG,1)
    seg = AXSEG{s,2};  ax = AXSEG{s,3};
    m = tg >= seg(1) & tg < seg(2) & ~gap;
    x = TAU(m,ax) - mean(TAU(m,ax));            % remove offset de trim
    nx = numel(x);
    [Pxx, ff] = pwelch(x, hann(min(nx, 128)), [], [], fs);
    subplot(3,1,s); hold on; grid on;
    plot(ff, 10*log10(Pxx + eps), 'LineWidth', 1.1);
    [~, ipk] = max(Pxx);
    xline(ff(ipk), 'r--', sprintf('pico %.2f Hz', ff(ipk)));
    ylabel('PSD (dB)');
    title(sprintf('%s — PSD de %s (energia na banda do modo, ~0 em DC)', AXSEG{s,1}, taulbl{ax}));
    if s==3, xlabel('frequencia (Hz)'); end
    xlim([0 fs/2]);
end
exportgraphics(fig3, fullfile(paths.images,'scan_excitation_psd.png'), 'BackgroundColor','white');

fprintf('\n  Figuras salvas em %s:\n', paths.images);
fprintf('    scan_excitation_overview.png  (RCIN | PWM | momentos | pqr)\n');
fprintf('    scan_excitation_snr.png       (SNR + isolamento por janela)\n');
fprintf('    scan_excitation_psd.png       (conteudo de frequencia dos momentos)\n\n');


%% ======================= HELPERS =======================
function s = okmark(tf)
    if tf, s = '  OK'; else, s = '  < min!'; end
end

function f = peak_freq(x, fs)
    x = x - mean(x);
    nx = numel(x);
    if nx < 8, f = NaN; return; end
    [Pxx, ff] = pwelch(x, hann(min(nx,128)), [], [], fs);
    [~, i] = max(Pxx);  f = ff(i);
end

function shade_axsegs(AX)
    yl = ylim;
    cols = [.2 .4 1; 1 .3 .3; .2 .7 .2];   % tom por eixo (1=roll,2=pitch,3=yaw)
    for i = 1:size(AX,1)
        seg = AX{i,2};  ax = AX{i,3};
        patch([seg(1) seg(2) seg(2) seg(1)], [yl(1) yl(1) yl(2) yl(2)], cols(ax,:), ...
              'FaceAlpha',0.06,'EdgeColor',[.5 .5 .5],'LineStyle',':','HandleVisibility','off');
    end
    ylim(yl);
end

function shade_gaps(tg, gap)
    yl = ylim;  d = diff([0; gap(:); 0]);
    s = tg(d(1:end-1) == 1);  e = tg(d(2:end) == -1);
    for i = 1:min(numel(s), numel(e))
        patch([s(i) e(i) e(i) s(i)], [yl(1) yl(1) yl(2) yl(2)], [.6 .6 .6], ...
              'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off');
    end
    ylim(yl);
end
