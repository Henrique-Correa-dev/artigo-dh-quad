% lag_sweep.m — Testa a hipótese de ERRO DE FASE (lag de motor) na malha de taxa
% =========================================================================
% A ACF do resíduo de p,q OSCILA (~0.6-0.7 Hz) → assinatura clássica de ERRO DE
% FASE: o modelo responde com fase ligeiramente diferente do real, e a diferença
% de dois sinais quase-iguais defasados é uma senoide na freq de excitação.
%
% Este script aplica um lag de 1ª ordem 1/(τs+1) ao COMANDO dos motores (PWM,
% antes da curva de empuxo → portanto no modelo de FORÇA) para uma grade de τ,
% simula em modo HYBRID com o MESMO P_final, e mede:
%   • % de brancura do resíduo (ACF dentro das bandas) por eixo
%   • TIC por eixo
% Se a brancura SOBE e o TIC CAI até um τ* e depois piora → erro de fase
% CONFIRMADO, e τ* é a constante de tempo do motor. Próximo passo: assar τ* no
% vtol_dynamics (saída T/Q) e RE-IDENTIFICAR (aqui P fica fixo: é só teste de fase).
%
% Também plota a PSD do resíduo (τ=0) p/ confirmar a frequência da oscilação.
%
% Uso:  >> lag_sweep

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---- CONFIG ----
LOG_FILE = 'logs_concat.mat';
t_val    = [605, 625];
dt       = 0.1;
SG_ORDER = 2;  SG_FRAME = 7;
TAU_GRID = [0, 0.05, 0.10, 0.15, 0.20, 0.30];   % s — 0 = sem lag (baseline)
MAXLAG   = 30;
P_FILE   = 'P_identified.mat';                   % usa P_final do full mode

%% ---- carregar P + log (MESMA filtragem do identify_v2 / lag_diagnostic) ----
proj = parameters();  constants = struct('m',proj.m,'g',proj.g);
d = load(fullfile(paths.outputs, P_FILE));  P = d.P_final;

L = load_log_data(fullfile(paths.data, LOG_FILE));
t0 = max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]);
t1 = min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]);
tg = (t0:dt:t1)';
ip = @(t,x) interp1(t,x,tg,'linear');
sg = @(x) sgolayfilt(x,SG_ORDER,SG_FRAME);
Pm=sg(ip(L.time_IMU,L.gyrX_raw)); Qm=sg(ip(L.time_IMU,L.gyrY_raw)); Rm=sg(ip(L.time_IMU,L.gyrZ_raw));
w1=sg(ip(L.time_RCOU,L.pwm1_raw)); w2=sg(ip(L.time_RCOU,L.pwm2_raw));
w3=sg(ip(L.time_RCOU,L.pwm3_raw)); w4=sg(ip(L.time_RCOU,L.pwm4_raw));
roll=ip(L.time_ATT,L.roll_deg); pitch=ip(L.time_ATT,L.pitch_deg); yaw=ip(L.time_ATT,L.yaw_deg);

m=(tg>=t_val(1))&(tg<=t_val(2));
time=tg(m); pqr=[Pm(m),Qm(m),Rm(m)]; att=[roll(m),pitch(m),yaw(m)];
W = [w1(m),w2(m),w3(m),w4(m)];

%% ---- SWEEP em τ (P fixo = teste de FASE) ----
fprintf('========================================================================\n');
fprintf('  SWEEP de lag de motor τ  [%g-%gs] — P FIXO (teste de fase)\n', t_val(1), t_val(2));
fprintf('========================================================================\n');
fprintf('  %6s | %16s | %16s | %16s\n', 'tau(s)', 'p: %branco/TIC', 'q: %branco/TIC', 'r: %branco/TIC');
fprintf('  %s\n', repmat('-',1,72));

nT = numel(TAU_GRID);
white = zeros(nT,3);  tic_v = zeros(nT,3);
for it = 1:nT
    tau = TAU_GRID(it);
    if tau > 0
        a = exp(-dt/tau);  lp = @(x) filter(1-a,[1 -a],x);   % 1ª ordem, DC=1
        Wl = [lp(W(:,1)), lp(W(:,2)), lp(W(:,3)), lp(W(:,4))];
    else
        Wl = W;
    end
    res = sim_window('hybrid', P, time, Wl, pqr, att, constants);
    sim = [res.p, res.q, res.r];
    for c = 1:3
        e = pqr(:,c) - sim(:,c);
        acf = acf_(e, MAXLAG);  N = numel(e);  b = 1.96/sqrt(N);
        white(it,c) = mean(abs(acf(2:end)) <= b) * 100;
        fm = fit_metrics(pqr(:,c), sim(:,c));
        tic_v(it,c) = fm.TIC;
    end
    fprintf('  %6.2f | %9.0f%% /%5.3f | %9.0f%% /%5.3f | %9.0f%% /%5.3f\n', ...
        tau, white(it,1),tic_v(it,1), white(it,2),tic_v(it,2), white(it,3),tic_v(it,3));
end

[~,bp]=max(white(:,1)); [~,bq]=max(white(:,2));
fprintf('\n  Melhor τ p/ p: %.2fs (%.0f%% branco) | q: %.2fs (%.0f%% branco)\n', ...
    TAU_GRID(bp), white(bp,1), TAU_GRID(bq), white(bq,2));
fprintf('  (brancura sobe + TIC cai com τ → erro de FASE confirmado; re-identifique com esse τ)\n');

%% ---- plot: brancura e TIC vs τ ----
fig=figure('Name','lag_sweep','Position',[60 60 1100 760]);
subplot(2,1,1); hold on; grid on;
plot(TAU_GRID*1000, white(:,1),'-o','LineWidth',1.5,'DisplayName','p');
plot(TAU_GRID*1000, white(:,2),'-s','LineWidth',1.5,'DisplayName','q');
plot(TAU_GRID*1000, white(:,3),'-^','LineWidth',1.5,'DisplayName','r');
yline(90,'k--','90% = branco'); ylabel('% lags dentro da banda'); legend('Location','best');
title('Brancura do resíduo vs lag de motor τ — pico = fase certa');
subplot(2,1,2); hold on; grid on;
plot(TAU_GRID*1000, tic_v(:,1),'-o','LineWidth',1.5,'DisplayName','p');
plot(TAU_GRID*1000, tic_v(:,2),'-s','LineWidth',1.5,'DisplayName','q');
plot(TAU_GRID*1000, tic_v(:,3),'-^','LineWidth',1.5,'DisplayName','r');
yline(0.30,'k--','0.30'); ylabel('TIC'); xlabel('τ (ms)'); legend('Location','best');
saveas(fig, fullfile(paths.images,'lag_sweep.png'));

%% ---- PSD do resíduo (τ=0) p/ confirmar a frequência da oscilação ----
res0 = sim_window('hybrid', P, time, W, pqr, att, constants);
sim0 = [res0.p, res0.q, res0.r];  names = {'p','q','r'};
fig2 = figure('Name','lag_sweep PSD','Position',[60 60 1100 760]);
for c = 1:3
    e = pqr(:,c) - sim0(:,c);  e = e - mean(e);  N = numel(e);
    w = 0.5*(1 - cos(2*pi*(0:N-1)'/(N-1)));      % janela de Hann (sem toolbox)
    Y = fft(e .* w);  P2 = abs(Y/N).^2;  P1 = P2(1:floor(N/2)+1);  P1(2:end-1) = 2*P1(2:end-1);
    f = (0:floor(N/2))' / (N*dt);
    [~,ipk] = max(P1(2:end));  fpk = f(ipk+1);
    subplot(3,1,c); plot(f, P1, 'LineWidth', 1.3); grid on; hold on;
    xline(fpk,'r--'); ylabel(['PSD ' names{c}]); xlim([0 2]);
    text(fpk+0.05, 0.8*max(P1), sprintf('%.2f Hz', fpk), 'Color','r');
    if c==1, title('PSD do resíduo (τ=0) — pico = freq da dinâmica não modelada'); end
    if c==3, xlabel('frequência (Hz)'); end
end
saveas(fig2, fullfile(paths.images,'lag_sweep_psd.png'));
fprintf('  Figuras: lag_sweep.png e lag_sweep_psd.png em %s\n', paths.images);

%% ---- helper: ACF normalizada (sem toolbox) ----
function acf = acf_(e, maxlag)
    e = e(:) - mean(e);  N = numel(e);  c0 = sum(e.^2) + 1e-12;
    acf = zeros(maxlag+1,1);
    for k = 0:maxlag
        acf(k+1) = sum(e(1:N-k) .* e(1+k:N)) / c0;
    end
end
