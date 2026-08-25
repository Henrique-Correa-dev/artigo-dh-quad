% excitation_quality.m — Minhas manobras EXCITAM o sistema? (roll/pitch)
% =========================================================================
% Critério clássico (Jategaonkar cap. input design / CIFER / Morelli):
%   1. PSD da ENTRADA espalhada pela banda de interesse (não só DC/1 pico)
%      → assinatura de bom multistep / 3-2-1-1 / sweep.
%   2. COERÊNCIA entrada↔saída γ²(f) ≥ 0.6 na banda (regra CIFER): abaixo disso
%      a saída não é explicada pela entrada (ruído, não-linearidade, ou pouca
%      excitação naquela frequência).
%   3. Banda excitada COBRE a dinâmica (polos de roll/pitch do modelo).
%
% Entrada = PWM DIFERENCIAL (independe do modelo — é a excitação real):
%   roll:  u_p = (pwm1+pwm4) - (pwm2+pwm3)   [direita - esquerda]
%   pitch: u_q = (pwm1+pwm3) - (pwm2+pwm4)   [frente  - trás]   (QuadX: 1FR 2RL 3FL 4RR)
% Saída  = p, q (gyro).
%
% ⚠️ Dados de MALHA FECHADA (Pixhawk estabilizando) → o PWM é em parte resposta
%    do controlador. A coerência aqui é um LIMITE INFERIOR da excitação real
%    (o feedback reduz coerência aparente). CIFER trata isso; EEM/OEM clássico
%    assume entrada de malha aberta.
%
% Uso:  >> excitation_quality

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---- CONFIG ----
LOG_FILE = 'logs_concat.mat';
WINDOWS  = {[4,24];[42,62];[100,125];[605,625]};  % treino (3) + validação (1)
dt       = 0.1;  fs = 1/dt;
SG_ORDER = 2;  SG_FRAME = 7;
COH_MIN  = 0.6;                 % limiar CIFER
BAND     = [0.2, 2.0];          % banda de interesse (Hz) p/ malha de taxa

%% ---- carregar + filtrar ----
d = load(fullfile(paths.data, LOG_FILE));  L = d.L;
t0 = max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]);
t1 = min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]);
tg = (t0:dt:t1)';
ip = @(t,x) interp1(t,x,tg,'linear');
sg = @(x) sgolayfilt(x,SG_ORDER,SG_FRAME);
P  = sg(ip(L.time_IMU,L.gyrX_raw));  Q = sg(ip(L.time_IMU,L.gyrY_raw));
R  = sg(ip(L.time_IMU,L.gyrZ_raw));
w1 = sg(ip(L.time_RCOU,L.pwm1_raw)); w2 = sg(ip(L.time_RCOU,L.pwm2_raw));
w3 = sg(ip(L.time_RCOU,L.pwm3_raw)); w4 = sg(ip(L.time_RCOU,L.pwm4_raw));
u_roll  = (w1+w4) - (w2+w3);    % direita - esquerda
u_pitch = (w1+w3) - (w2+w4);    % frente  - trás
u_yaw   = (w1+w2) - (w3+w4);    % CCW(1,2) - CW(3,4)  (= padrão do Mz do modelo)

%% ---- por janela ----
fprintf('========================================================================\n');
fprintf('  QUALIDADE DE EXCITAÇÃO (roll/pitch) — coerência CIFER γ²>=%.1f, banda %.1f-%.1f Hz\n', ...
    COH_MIN, BAND(1), BAND(2));
fprintf('========================================================================\n');
fprintf('  %-11s %6s | %-15s | %-15s | %-15s\n', ...
    'janela','tipo','ROLL (u_p→p)','PITCH (u_q→q)','YAW (u_r→r)');
fprintf('  %-11s %6s | %7s %7s | %7s %7s | %7s %7s\n', ...
    '','','γ²méd','%band>','γ²méd','%band>','γ²méd','%band>');
fprintf('  %s\n', repmat('-',1,74));

nW = numel(WINDOWS);
fig = figure('Name','excitation_quality','Position',[40 40 1500 820]);
axes_in  = {u_roll, u_pitch, u_yaw};
axes_out = {P, Q, R};
axes_ttl = {'ROLL  u_p \rightarrow p', 'PITCH  u_q \rightarrow q', 'YAW  u_r \rightarrow r'};
axes_col = {'b','m',[0 0.6 0]};
axes_lbl = {'roll','pitch','yaw'};
for iw = 1:nW
    tw = WINDOWS{iw};  m = (tg>=tw(1))&(tg<=tw(2));
    tipo = 'treino'; if iw==nW, tipo='valid'; end
    g = zeros(1,3); pct = zeros(1,3);
    for ax = 1:3
        [f, coh] = mscoh_(axes_in{ax}(m), axes_out{ax}(m), fs);
        bandm = (f>=BAND(1)) & (f<=BAND(2));
        g(ax) = mean(coh(bandm));  pct(ax) = 100*mean(coh(bandm)>=COH_MIN);
        subplot(nW,3,3*(iw-1)+ax); hold on; grid on;
        plot(f, coh, 'Color', axes_col{ax}, 'LineWidth',1.2);
        yline(COH_MIN,'r--'); xline(BAND(1),'k:'); xline(BAND(2),'k:');
        xlim([0 fs/2]); ylim([0 1]);
        if ax==1, ylabel(sprintf('[%g-%gs]\n\\gamma^2',tw(1),tw(2))); end
        if iw==1, title(['Coerência ' axes_ttl{ax}]); end
        if iw==nW, xlabel('f (Hz)'); end
    end
    fprintf('  [%4.0f,%4.0f]s %6s | %6.2f %6.0f%% | %6.2f %6.0f%% | %6.2f %6.0f%%\n', ...
        tw(1),tw(2),tipo, g(1),pct(1), g(2),pct(2), g(3),pct(3));
end
saveas(fig, fullfile(paths.images,'excitation_quality.png'));

fprintf('\n  Como ler:\n');
fprintf('   • γ²méd >= %.1f e %%band> alto → entrada EXCITA bem aquele eixo (CIFER ok).\n', COH_MIN);
fprintf('   • γ² baixo → pouca excitação / ruído / não-linearidade / efeito do feedback.\n');
fprintf('   • revers = nº de reversões de sinal da entrada (proxy de riqueza tipo 3-2-1-1).\n');
fprintf('   Figura: %s\n', fullfile(paths.images,'excitation_quality.png'));

%% ---- helper: PSD + coerência de Welch (Hann, 50%% overlap) sem toolbox ----
function [f, coh, Pxx_out] = mscoh_(x, y, fs)
    x = x(:)-mean(x);  y = y(:)-mean(y);  N = numel(x);
    L = 2^floor(log2(N/4));  if L < 16, L = 2^max(4,floor(log2(N/2))); end
    ov = floor(L/2);  step = L - ov;  win = 0.5*(1-cos(2*pi*(0:L-1)'/(L-1)));
    nseg = max(1, floor((N - ov)/step));
    Pxx=zeros(L,1); Pyy=zeros(L,1); Pxy=zeros(L,1);
    for i = 1:nseg
        idx = (i-1)*step + (1:L);
        X = fft(x(idx).*win);  Y = fft(y(idx).*win);
        Pxx = Pxx + X.*conj(X);  Pyy = Pyy + Y.*conj(Y);  Pxy = Pxy + X.*conj(Y);
    end
    coh = abs(Pxy).^2 ./ (Pxx.*Pyy + eps);
    half = 1:floor(L/2)+1;
    f = (0:floor(L/2))' * fs / L;
    coh = coh(half);  Pxx_out = Pxx(half)/nseg;
end
