% lag_diagnostic.m — ATRASO (delay) entre modelo e medido por cross-correlação
% =========================================================================
% Você notou um delay sim×medido em p,q. Este script mede o atraso por
% cross-correlação e ajuda a decidir a causa:
%   • ~0.5 amostra (~50 ms), com arredondamento  → LAG DE MOTOR/ROTOR (1ª ordem)
%   • 1–2 amostras IGUAIS nos 3 eixos, shift puro → DESSINCRONIA RCOU↔IMU (timestamp)
%
% Convenção: lag > 0  ⇒  MEDIDO chega DEPOIS do SIM (modelo responde antes).
%
% Uso:  >> lag_diagnostic

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---- CONFIG ----
LOG_FILE = 'logs_concat.mat';
t_val    = [605, 625];
dt       = 0.1;
SG_ORDER = 2;  SG_FRAME = 7;
MAXLAG   = 15;                 % amostras (= 1.5 s a 10 Hz)
P_FILE   = 'P_identified.mat'; % usa P_final do full mode

%% ---- carregar P + log (mesma filtragem do identify_v2) ----
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
time=tg(m);  pqr=[Pm(m),Qm(m),Rm(m)];  pwm=[w1(m),w2(m),w3(m),w4(m)];  att=[roll(m),pitch(m),yaw(m)];

%% ---- simular (modo hybrid: pqr integrado, att medida) ----
res = sim_window('hybrid', P, time, pwm, pqr, att, constants);
simv = [res.p, res.q, res.r];

%% ---- cross-correlação por eixo ----
names = {'p','q','r'};
fprintf('========================================================\n');
fprintf('  ATRASO modelo×medido [%g-%gs] | dt=%.0f ms\n', t_val(1),t_val(2),dt*1000);
fprintf('========================================================\n');
fprintf('  %-4s %12s %10s %9s %9s\n','eixo','lag(amostr)','lag(ms)','xc@lag','xc@0');
fig=figure('Name','lag_diagnostic','Position',[60 60 1100 720]);
for i=1:3
    [xc,lags] = xcorr_norm(pqr(:,i), simv(:,i), MAXLAG);
    [xcmax,idx]=max(xc);  klag=lags(idx);  xc0=xc(lags==0);
    fprintf('  %-4s %12d %10.0f %9.3f %9.3f\n', names{i}, klag, klag*dt*1000, xcmax, xc0);
    subplot(3,1,i); hold on; grid on;
    stem(lags,xc,'filled','MarkerSize',3,'Color',[.3 .6 1]);
    plot(klag,xcmax,'ro','MarkerFaceColor','r');
    xline(0,'--','Color',[.6 .6 .6]);
    ylabel(['xcorr ' names{i}]); ylim([-0.2 1]);
    if i==1, title(sprintf('Cross-correlação medido×sim [%g-%gs] — pico fora de 0 = atraso', t_val(1),t_val(2))); end
    if i==3, xlabel('lag (amostras)  —  >0: MEDIDO atrasado em relação ao SIM'); end
end
saveas(fig, fullfile(paths.images,'lag_diagnostic.png'));

fprintf('\n  Como ler:\n');
fprintf('   • lag ~0.5 amostra (~50 ms) e parecido nos eixos excitados → LAG DE MOTOR.\n');
fprintf('   • lag de 1-2 amostras IGUAL nos 3 eixos → DESSINCRONIA RCOU↔IMU (corrige com shift).\n');
fprintf('   • xc@lag >> xc@0 → o atraso é real e significativo.\n');
fprintf('   Figura: %s\n', fullfile(paths.images,'lag_diagnostic.png'));


%% ---- helper: cross-correlação normalizada (sem toolbox) ----
function [xc,lags] = xcorr_norm(x,y,maxlag)
% xc(k) = corr( x(t), y(t-k) ).  k>0 ⇒ x casa com y atrasado ⇒ x ATRASADO vs y.
    x=x(:)-mean(x);  y=y(:)-mean(y);  N=numel(x);
    den=max(sqrt(sum(x.^2))*sqrt(sum(y.^2)),1e-12);
    lags=(-maxlag:maxlag)';  xc=zeros(numel(lags),1);
    for j=1:numel(lags)
        k=lags(j);
        if k>=0, a=x(1+k:N); b=y(1:N-k);
        else,    a=x(1:N+k); b=y(1-k:N);
        end
        xc(j)=sum(a.*b)/den;
    end
end
