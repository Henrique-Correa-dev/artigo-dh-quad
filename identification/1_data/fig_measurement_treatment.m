% fig_measurement_treatment.m — Figuras da etapa de Measurements (Cap3 da tese)
% =========================================================================
% Gera TRÊS figuras, cada uma ilustrando um aspecto do tratamento das medidas:
%
%  (1) medidas_tratamento.png — ACELERÔMETRO (ax,ay,az): bruto NATIVO (~25 Hz)
%      vs tratado (reamostrado 10 Hz, sem bias, filtrado SG). Janelas com MOTORES
%      LIGADOS (vibração) → a supressão de ruído pela filtragem fica visível.
%      Esquerda: manobra; direita: zoom em voo pairado (mostra também 25 Hz vs 10 Hz).
%
%  (2) medidas_bias.png — calibração com o DRONE PARADO (motores desligados),
%      GIRO (esq.) e ACELERÔMETRO (dir.). Giro: bruto (no bias) vs tratado (em 0).
%      Acel.: plota a-g_proj (força específica menos a gravidade projetada pela
%      atitude) — o bruto fica no bias, o tratado em 0. Calibra os DOIS sensores
%      na mesma janela.
%
%  (3) medidas_espectro.png — GIRO (rolagem): taxa e aceleração angular, não-filtrado
%      vs filtrado, no tempo e no PSD (a diferenciação amplifica o ruído de alta freq.).
%
% Bias do giro = média das taxas na janela de calibração (esperado = 0, drone parado).
% Bias do acel. = média(acc) - (-g_body), com g_body projetado pela atitude do EKF.
%
% Uso:  >> fig_measurement_treatment

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();
set(groot, 'defaultAxesFontSize', 11);

%% ---------------- CONFIG ----------------
LOG_FILE = 'logs_concat.mat';
dt = 0.1;  fs = 1/dt;                   % grade comum (10 Hz)
SG_ORDER = 2;  SG_FRAME = 7;            % Savitzky-Golay (mesma da tese)
g  = 9.81;
W_MAN  = [605 625];                     % manobra (filtragem da vibração visível)
W_HOV  = [364 369];                     % voo pairado CALMO (log6) p/ o zoom + bias em voo
W_CAL  = [136 150];                     % drone parado p/ calibração (mede o bias)
W_CALZ = [139 143];                     % zoom da calibração
AX_B   = 1;  AXNM_B = 'p (rolagem)';    % eixo da fig. de espectro
W_SPEC = [42 62];                       % janela contígua p/ a PSD
W_TB   = [49 54];                       % zoom no tempo (um doublet) p/ a fig. de espectro
DEST = {paths.images, '/Users/graest/ita-master/artigo/artigo-dh-quad/dissertacao/Cap3'};

%% ---------------- CARREGAR + BIAS (giro e acelerômetro) ----------------
L = load_log_data(fullfile(paths.data, LOG_FILE));
tIMU = L.time_IMU;
fs_imu = 1/median(diff(tIMU));                       % taxa nativa do IMU (~25 Hz)
gyrN = {L.gyrX_raw, L.gyrY_raw, L.gyrZ_raw};
accN = {L.accX_raw, L.accY_raw, L.accZ_raw};

mcal = (tIMU >= W_CAL(1)) & (tIMU <= W_CAL(2));       % drone parado (IMU)
macA = (L.time_ATT >= W_CAL(1)) & (L.time_ATT <= W_CAL(2));
% bias do giro (esperado = 0)
biasG = cellfun(@(x) mean(x(mcal)), gyrN(:));
% gravidade projetada no corpo pela atitude do EKF (força específica esperada em repouso = -g_body)
phi_c = deg2rad(mean(L.roll_deg(macA)));  th_c = deg2rad(mean(L.pitch_deg(macA)));
exp_acc = [ g*sin(th_c); -g*sin(phi_c)*cos(th_c); -g*cos(phi_c)*cos(th_c) ];
biasA = cellfun(@(x) mean(x(mcal)), accN(:)) - exp_acc;   % bias do acelerômetro (solo)
% bias do acelerômetro EM VOO, calibrado num hover CALMO e nivelado: lá a_x,a_y
% devem ser ~0 e a_z ~ -g, então o offset medido é o bias (simples, sem projeção do EKF).
mhov = (tIMU>=W_HOV(1)) & (tIMU<=W_HOV(2));
biasA_fl = [mean(accN{1}(mhov)); mean(accN{2}(mhov)); mean(accN{3}(mhov))+g];
fprintf('IMU ~%.1f Hz | bias giro=[%+.5f %+.5f %+.5f] | bias acel solo=[%+.3f %+.3f %+.3f] voo=[%+.3f %+.3f %+.3f]\n', ...
    fs_imu, biasG, biasA, biasA_fl);

%% ---------------- CADEIA DE TRATAMENTO ----------------
t0 = max([tIMU(1), L.time_RCOU(1)]);  t1 = min([tIMU(end), L.time_RCOU(end)]);
tg = (t0:dt:t1)';
if exist('sgolayfilt','file') == 2, sgf = @(x) sgolayfilt(x, SG_ORDER, SG_FRAME);
else, warning('sgolayfilt ausente — SG local'); sgf = @(x) sg_local(x, SG_ORDER, SG_FRAME); end
resG = cell(3,1); filG = cell(3,1); resA = cell(3,1); filA = cell(3,1); filAdb = cell(3,1);
for a = 1:3
    resG{a} = interp1(tIMU, gyrN{a}, tg, 'linear');   % giro reamostrado (10 Hz)
    filG{a} = sgf(resG{a} - biasG(a));                 % giro: sem bias (estável) + filtrado
    resA{a} = interp1(tIMU, accN{a}, tg, 'linear');    % accel reamostrado (10 Hz)
    filA{a} = sgf(resA{a});                            % accel reamostrado + filtrado (bias removido por janela no plot)
    filAdb{a} = sgf(resA{a} - biasA(a));               % accel debiased — só p/ ilustrar a calibração no solo
end
lab_b = sprintf('bruto (%.0f Hz)', fs_imu);

%% ===== FIG 1: TRATAMENTO — taxas angulares (esq.) e acelerômetro (dir.) =====
% Mesma janela de manobra nas duas colunas: esquerda p,q,r e direita a_x,a_y,a_z,
% bruto (nativo ~25 Hz) vs tratado (reamostrado 10 Hz, sem bias, filtrado SG).
nmG = {'p [rad/s]','q [rad/s]','r [rad/s]'};
nmA = {'a_x [m/s^2]','a_y [m/s^2]','a_z [m/s^2]'};
colG = {'b','r',[0 .6 0]};  colA = {'b','r',[0 .6 0]};
figA = figure('Color','w','Position',[40 40 1450 850]);
mn1 = (tIMU>=W_MAN(1)) & (tIMU<=W_MAN(2));   % nativo (~25 Hz) na janela de manobra
m   = (tg>=W_MAN(1))   & (tg<=W_MAN(2));      % grade comum (10 Hz)
for a = 1:3
    % ---- ESQUERDA: taxas angulares p, q, r ----
    subplot(3,2,2*a-1);
    bLg = mean(filG{a}(m));   % bias do segmento (média) — centra a taxa em 0
    plot(tIMU(mn1), gyrN{a}(mn1),   '-', 'Color', [.6 .6 .6], 'LineWidth', 0.7); hold on;  % bruto (com bias)
    plot(tg(m),     filG{a}(m)-bLg, '-', 'Color', colG{a}, 'LineWidth', 1.4);              % tratado (sem bias + filtrado)
    grid on; ylabel(nmG{a});
    if a==1, title('Taxas angulares — bruto vs tratado');
        legend({lab_b,'tratado (10 Hz)'},'Location','best'); end
    if a==3, xlabel('t [s]'); end
    % ---- DIREITA: acelerômetro a_x, a_y, a_z ----
    subplot(3,2,2*a);
    bLa = mean(filA{a}(m)); if a==3, bLa = bLa + g; end   % bias do segmento — centra a_x,a_y em 0 e a_z em -g
    plot(tIMU(mn1), accN{a}(mn1),   '-', 'Color', [.6 .6 .6], 'LineWidth', 0.7); hold on;  % bruto (com bias)
    plot(tg(m),     filA{a}(m)-bLa, '-', 'Color', colA{a}, 'LineWidth', 1.4);              % tratado (sem bias + filtrado)
    grid on; ylabel(nmA{a});
    if a==1, title('Acelerações — bruto vs tratado');
        legend({lab_b,'tratado (10 Hz)'},'Location','best'); end
    if a==3, xlabel('t [s]'); end
end
% (sem título — a legenda do LaTeX descreve a figura)
save_to(figA, 'medidas_tratamento.png', DEST);

%% ===== FIG 2: CALIBRAÇÃO DE BIAS (giro + acelerômetro, drone parado) =====
nmG  = {'p [rad/s]','q [rad/s]','r [rad/s]'};
nmAr = {'a_x [m/s^2]','a_y [m/s^2]','a_z [m/s^2]'};
colG = {'b','r',[0 .6 0]};
figB = figure('Color','w','Position',[40 40 1450 850]);
mn = (tIMU>=W_CALZ(1)) & (tIMU<=W_CALZ(2));  mq = (tg>=W_CALZ(1)) & (tg<=W_CALZ(2));
for a = 1:3
    % esquerda: giro (esperado em repouso = 0)
    subplot(3,2,2*a-1);
    plot(tIMU(mn), gyrN{a}(mn), '.',  'Color', [.55 .55 .55], 'MarkerSize', 10); hold on;
    plot(tg(mq),   filG{a}(mq), '-o', 'Color', colG{a}, 'LineWidth', 1.3, 'MarkerSize', 4, 'MarkerFaceColor', colG{a});
    yline(0, ':k', 'HandleVisibility', 'off');
    grid on; ylabel(nmG{a});
    text(0.02, 0.86, sprintf('bias = %+.4f rad/s', biasG(a)), 'Units','normalized', 'FontSize', 9, 'BackgroundColor', [1 1 1 .6]);
    if a==1, title(sprintf('Giroscópio — bruto (%.0f Hz) vs tratado (10 Hz)', fs_imu));
        legend({lab_b,'tratado (sem bias)'},'Location','best'); end
    if a==3, xlabel('t [s]'); end
    % direita: acelerômetro, resíduo a - g_proj (esperado em repouso = 0)
    subplot(3,2,2*a);
    plot(tIMU(mn), accN{a}(mn)-exp_acc(a),   '.',  'Color', [.55 .55 .55], 'MarkerSize', 10); hold on;
    plot(tg(mq),   filAdb{a}(mq)-exp_acc(a), '-o', 'Color', colA{a}, 'LineWidth', 1.3, 'MarkerSize', 4, 'MarkerFaceColor', colA{a});
    yline(0, ':k', 'HandleVisibility', 'off');
    grid on; ylabel(nmAr{a});
    text(0.02, 0.86, sprintf('bias = %+.3f m/s^2', biasA(a)), 'Units','normalized', 'FontSize', 9, 'BackgroundColor', [1 1 1 .6]);
    if a==1, title(sprintf('Acelerômetro (a - g projetada) — bruto (%.0f Hz) vs tratado (10 Hz)', fs_imu));
        legend({lab_b,'tratado (sem bias)'},'Location','best'); end
    if a==3, xlabel('t [s]'); end
end
% (sem título — a legenda do LaTeX descreve a figura)
save_to(figB, 'medidas_bias.png', DEST);

%% ===== FIG 3: FILTRO + DIFERENCIAÇÃO (giro, eixo de rolagem) =====
a = AX_B;
rate_nf = resG{a};                 rate_f = sgf(resG{a});
acc_nf  = gradient(resG{a}, dt);   acc_f  = gradient(sgf(resG{a}), dt);
figC = figure('Color','w','Position',[80 80 1350 780]);
mt = (tg>=W_TB(1)) & (tg<=W_TB(2));
subplot(2,2,1);
plot(tg(mt), rate_nf(mt), '-', 'Color', [.6 .6 .6], 'LineWidth', 0.8); hold on;
plot(tg(mt), rate_f(mt),  'b', 'LineWidth', 1.5); grid on;
ylabel('taxa [rad/s]'); title('Taxa angular — tempo'); legend({'não filtrado','filtrado'},'Location','best');
subplot(2,2,3);
plot(tg(mt), acc_nf(mt), '-', 'Color', [.75 .75 .75], 'LineWidth', 0.8); hold on;
plot(tg(mt), acc_f(mt),  'b', 'LineWidth', 1.5); grid on;
ylabel('acel. ang. [rad/s^2]'); xlabel('t [s]');
title('Aceleração angular — a diferenciação amplifica o ruído');
legend({'derivada do não filtrado','derivada do filtrado'},'Location','best');
ms = (tg>=W_SPEC(1)) & (tg<=W_SPEC(2));
[fr, Pr] = welch_psd(rate_nf(ms), fs);  [~, Pf]  = welch_psd(rate_f(ms), fs);
[~, Par] = welch_psd(acc_nf(ms),  fs);  [~, Paf] = welch_psd(acc_f(ms),  fs);
subplot(2,2,2);
loglog(fr, Pr, '-', 'Color', [.6 .6 .6], 'LineWidth', 1.0); hold on;
loglog(fr, Pf, 'b', 'LineWidth', 1.4); grid on;
ylabel('PSD taxa'); title('Espectro da taxa angular'); legend({'não filtrado','filtrado'},'Location','southwest'); xlim([fr(2) fs/2]);
subplot(2,2,4);
loglog(fr, Par, '-', 'Color', [.75 .75 .75], 'LineWidth', 1.0); hold on;
loglog(fr, Paf, 'b', 'LineWidth', 1.4); grid on;
ylabel('PSD acel. ang.'); xlabel('f [Hz]'); title('Espectro da aceleração angular'); xlim([fr(2) fs/2]);
% (sem título — a legenda do LaTeX descreve a figura)
save_to(figC, 'medidas_espectro.png', DEST);

%% ===== FIG 4: ESPECTRO DO ACELERÔMETRO (vibração) =====
% Voo pairado, motores LIGADOS, janela contígua. O bruto NATIVO (25 Hz) alcança
% 12,5 Hz e revela a vibração; o tratado (10 Hz + SG) é limitado a 5 Hz (Nyquist).
WP = [353 369];   WT = [354 358];
nmA2 = {'a_x','a_y','a_z'};
% aviso se a janela cruzar gap entre logs (corromperia a PSD)
mnp = (tIMU>=WP(1)) & (tIMU<=WP(2));
if max(diff(tIMU(mnp))) > 0.2
    warning('janela da PSD [%g %g] cruza um gap — escolha outra (W P).', WP);
end
mqp = (tg>=WP(1))   & (tg<=WP(2));
mnt = (tIMU>=WT(1)) & (tIMU<=WT(2));   mqt = (tg>=WT(1)) & (tg<=WT(2));
figD = figure('Color','w','Position',[60 60 1350 850]);
for a = 1:3
    % tempo (zoom)
    subplot(3,2,2*a-1);
    plot(tIMU(mnt), accN{a}(mnt), '-', 'Color', [.6 .6 .6], 'LineWidth', 0.7); hold on;
    plot(tg(mqt),   filA{a}(mqt), '-', 'Color', colA{a}, 'LineWidth', 1.4);
    grid on; ylabel(sprintf('%s [m/s^2]', nmA2{a}));
    if a==1, title('Voo pairado — domínio do tempo'); legend({lab_b,'tratado (10 Hz)'},'Location','best'); end
    if a==3, xlabel('t [s]'); end
    % PSD (Welch): nativo 25 Hz vs tratado 10 Hz
    subplot(3,2,2*a);
    [fN, PN] = welch_psd(accN{a}(mnp), fs_imu);
    [fT, PT] = welch_psd(filA{a}(mqp), fs);
    loglog(fN, PN, '-', 'Color', [.6 .6 .6], 'LineWidth', 1.0); hold on;
    loglog(fT, PT, '-', 'Color', colA{a}, 'LineWidth', 1.5);
    xline(fs/2, 'k--', 'Nyquist 10 Hz', 'HandleVisibility','off', 'LabelVerticalAlignment','bottom', 'FontSize',8);
    grid on; ylabel(sprintf('PSD %s [(m/s^2)^2/Hz]', nmA2{a}));
    if a==1, title('Espectro (Welch)'); legend({lab_b,'tratado (10 Hz)'},'Location','southwest'); end
    if a==3, xlabel('f [Hz]'); end
    xlim([fN(2) fs_imu/2]);
end
% (sem título — a legenda do LaTeX descreve a figura)
save_to(figD, 'medidas_espectro_acc.png', DEST);

fprintf('Pronto.\n');
set(groot, 'defaultAxesFontSize', 'remove');

%% ---------------- HELPERS ----------------
function save_to(fig, name, dests)
    lighten(fig);
    for i = 1:numel(dests)
        if ~exist(dests{i}, 'dir'), mkdir(dests{i}); end
        p = fullfile(dests{i}, name);
        try, exportgraphics(fig, p, 'Resolution', 200); catch, saveas(fig, p); end
    end
    fprintf('  salvo: %s\n', name);
end

function lighten(fig)
    set(fig, 'Color', 'w');
    set(findall(fig,'Type','axes'), 'Color','w', 'XColor','k', 'YColor','k', ...
        'GridColor', [.15 .15 .15], 'MinorGridColor', [.3 .3 .3]);
    set(findall(fig,'Type','text'), 'Color', 'k');
    lg = findall(fig,'Type','legend'); for i=1:numel(lg), set(lg(i),'TextColor','k','Color','w'); end
    set(findall(fig,'Type','constantline'), 'Color', 'k');
end

function [f, P] = welch_psd(x, fs)
    x = x(:) - mean(x); N = numel(x);
    Lw = 2^floor(log2(max(N/4, 8)));  if Lw < 16, Lw = min(16, 2^floor(log2(N))); end
    ov = floor(Lw/2); step = max(Lw - ov, 1);
    win = 0.5*(1 - cos(2*pi*(0:Lw-1)'/(Lw-1)));  U = sum(win.^2);
    nseg = max(1, floor((N - ov)/step));  P = zeros(Lw, 1); cnt = 0;
    for i = 1:nseg
        idx = (i-1)*step + (1:Lw);  if idx(end) > N, break; end
        X = fft(x(idx).*win);  P = P + abs(X).^2;  cnt = cnt + 1;
    end
    P = P / (max(cnt,1) * fs * U);
    half = 1:floor(Lw/2)+1;  f = (0:floor(Lw/2))' * fs / Lw;  P = P(half);
    P(2:end-1) = 2*P(2:end-1);
end

function y = sg_local(x, order, frame)
    x = x(:); m = (frame-1)/2; t = (-m:m)'; A = t.^(0:order);
    Hc = A * pinv(A); c = Hc(m+1, :)'; y = conv(x, flipud(c), 'same');
end
