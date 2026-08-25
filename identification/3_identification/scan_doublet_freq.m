% scan_doublet_freq.m — Isola UM doublet por eixo e mostra o conteudo de frequencia.
% =========================================================================
% Demonstra o criterio de ENTRADA da manobra (conteudo de frequencia / Delta t,
% energia na banda): rolagem e arfagem atendem (energia ~1 Hz, na banda da
% dinamica de taxa); a guinada e o canal LENTO/fraco (tau_z so de reacao) ->
% so foi possivel um doublet lento, com energia em ~0.2 Hz.
%
% Para cada eixo:
%   - localiza o doublet de maior excursao no segmento (pico + pico oposto adjacente);
%   - mostra o doublet isolado no tempo: momento tau (entrada) + taxa (resposta);
%   - mostra a PSD do momento do doublet: pico de energia + banda ~1:3, ~0 em DC.
%
% Momentos calculados com a MESMA alocacao do OEM (vtol_dynamics). So diagnostico.
% Uso:  >> scan_doublet_freq

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---------------- CONFIG ----------------
LOG_FILE = 'logs_concat.mat';
dt = 0.1;  fs = 1/dt;

% {nome, [t0 t1] do segmento, eixo(1=roll,2=pitch,3=yaw), simbolo tau, simbolo taxa}
AXSEG = { 'Rolagem (roll)',  [ 40  80], 1, '\tau_x', 'p';
          'Arfagem (pitch)', [ 80 125], 2, '\tau_y', 'q';
          'Guinada (yaw)',   [  0  40], 3, '\tau_z', 'r' };

%% ---------------- MODELO DE MOMENTO (idêntico ao identify_plant) ----------------
% Replica o pré-processo do OEM: (1) atraso de transporte de 1 amostra no PWM,
% (2) lag de 1a ordem (spin-up do motor, tau_m) no PWM == lag na RPM, e (3) modelo
% de motor FISICO (PWM->RPM->k*RPM^2). Sem isso o momento sai ~150 ms adiantado da
% resposta e ainda aplica os k_T/k_Q (calibrados no fisico) sobre a curva spline.
DELAY_PWM = 1;                               % amostras (atraso puro, = identify_plant)
proj = parameters();
tau_lag = proj.motor.tau_m;                  % s — lag 1a ordem (USE_MOTOR_MODEL=true no OEM)
[fT, fQ] = motor_models('physical');         % T,Q = k*(CR*PWM+w_b)^2  (= identify_plant)
dyn = vtol_dynamics('get_handles');
lmf = fullfile(paths.outputs, 'linear_model.mat');
if exist(lmf, 'file'), S = load(lmf);  P = S.P;  else, P = proj.P0_J;  end

%% ---------------- CARREGAR + GRADE ----------------
L  = load_log_data(fullfile(paths.data, LOG_FILE));
t0 = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t1 = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t0:dt:t1)';
ip = @(tt, xx) interp1(tt, xx, tg, 'linear');
PQR = [ip(L.time_IMU, L.gyrX_raw), ip(L.time_IMU, L.gyrY_raw), ip(L.time_IMU, L.gyrZ_raw)];
W   = [ip(L.time_RCOU, L.pwm1_raw), ip(L.time_RCOU, L.pwm2_raw), ...
       ip(L.time_RCOU, L.pwm3_raw), ip(L.time_RCOU, L.pwm4_raw)];
% atraso (1 amostra) + lag (tau_m) no PWM, idêntico ao identify_plant (por coluna)
dl    = @(x) [repmat(x(1),DELAY_PWM,1); x(1:end-DELAY_PWM)];
a_lag = exp(-dt/tau_lag);
lp    = @(x) filter(1-a_lag, [1 -a_lag], x, a_lag*x(1));
for j = 1:4, W(:,j) = lp(dl(W(:,j))); end
[~, tx, ty, tz] = dyn.forces(W, P, fT, fQ);
TAU = [tx ty tz];

%% ---------------- FIGURA: doublet isolado por eixo (tempo | frequencia) ----------------
fig = figure('Name','scan_doublet_freq', 'Position',[50 50 1200 880]);
set(fig, 'Color', 'w');  try, fig.Theme = 'light'; catch, end

fprintf('\n  eixo              fpico(Hz)   Dt(s)    std_tau   std_taxa\n');
fprintf('  --------------------------------------------------------------\n');
for s = 1:3
    seg = AXSEG{s,2};  ax = AXSEG{s,3};
    in  = tg >= seg(1) & tg < seg(2);
    tt  = tg(in);
    xx  = TAU(in,ax) - mean(TAU(in,ax));     % momento (entrada), sem trim
    yy  = PQR(in,ax) - mean(PQR(in,ax));     % taxa (resposta), sem trim

    % escala de tempo do canal: freq dominante do segmento (PSD grosseira)
    [Ps, ffs] = pwelch(xx, hann(min(numel(xx),128)), [], [], fs);
    [~, isd]  = max(Ps);  fseg = max(ffs(isd), 0.05);
    Tseg = 1/fseg;  Wd = 3*Tseg;  rad = Tseg;   % janela ~3 periodos; busca do par dentro de 1 periodo

    % localizar doublet: maior pico e o maior pico de sinal OPOSTO dentro de 1 periodo
    [~, i1] = max(abs(xx));  sgn = sign(xx(i1));
    opp = (abs(tt - tt(i1)) <= rad) & (sign(xx) ~= sgn);
    if any(opp)
        tmp = xx;  tmp(~opp) = 0;  [~, i2] = max(abs(tmp));
    else
        i2 = i1;
    end
    tc = (tt(i1) + tt(i2)) / 2;               % centro entre os dois pulsos
    dt_meas = abs(tt(i2) - tt(i1));           % Delta t medido (largura do pulso)

    sel = tt >= tc - Wd/2 & tt <= tc + Wd/2;
    td  = tt(sel) - tc;  xd = xx(sel) - mean(xx(sel));  yd = yy(sel) - mean(yy(sel));

    % PSD do doublet isolado
    [Pxx, ff] = pwelch(xd, hann(numel(xd)), [], 512, fs);
    [~, im]   = max(Pxx);  fpk = ff(im);
    band = [fpk/sqrt(3), fpk*sqrt(3)];        % banda util ~1:3 em torno do pico

    fprintf('  %-16s %8.2f %8.2f %9.4f %9.4f\n', AXSEG{s,1}, fpk, dt_meas, std(xd), std(yd));

    % --- TEMPO: momento (entrada) + taxa (resposta) ---
    subplot(3,2,2*s-1); hold on; grid on;
    yyaxis left;  plot(td, xd, '-',  'LineWidth', 1.4);  ylabel(sprintf('%s (N\\cdotm)', AXSEG{s,4}));
    yyaxis right; plot(td, yd, '--', 'LineWidth', 1.0);  ylabel(sprintf('%s (rad/s)', AXSEG{s,5}));
    xlabel('t (s)');  xlim([td(1) td(end)]);
    title(sprintf('%s — doublet isolado (\\Deltat \\approx %.2f s)', AXSEG{s,1}, dt_meas));

    % --- FREQUENCIA: PSD do momento ---
    subplot(3,2,2*s); hold on; grid on;
    yl = [min(10*log10(Pxx+eps)), max(10*log10(Pxx+eps))+3];
    patch([band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [.7 .9 .7], 'FaceAlpha',0.25, 'EdgeColor','none', 'HandleVisibility','off');
    plot(ff, 10*log10(Pxx+eps), '-', 'LineWidth', 1.3);
    xline(fpk, 'r--', sprintf('%.2f Hz', fpk));
    xlim([0 3]);  ylim(yl);  ylabel('PSD (dB)');  xlabel('frequência (Hz)');
    title(sprintf('%s — energia em %.2f Hz (banda \\approx 1:3)', AXSEG{s,1}, fpk));
end

exportgraphics(fig, fullfile(paths.images,'scan_doublet_freq.png'), 'BackgroundColor','white');
fprintf('\n  Figura salva: %s\n', fullfile(paths.images,'scan_doublet_freq.png'));
