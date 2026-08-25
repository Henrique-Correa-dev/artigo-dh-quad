% scan_excitation_rc.m — Excitação de roll/pitch/yaw via RCIN e RCOU
% =========================================================================
% Justifica que as manobras (doublets) excitam bem os 3 eixos em MALHA FECHADA,
% avaliando as DUAS entradas relevantes:
%   RCIN  = comando exógeno do piloto (rádio)  -> o "d_p" do cap.9 (independe dos estados)
%   RCOU  = PWM enviado aos motores            -> a entrada EFETIVA da planta "u"
%
% Critérios (Jategaonkar cap.2/cap.9, Klein & Morelli, CIFER):
%   (1) Excitação exógena POR EIXO, um eixo de cada vez  -> reduz colinearidade (cap.9)
%   (2) RCOU diferencial acompanha o comando             -> a planta de fato recebe o momento
%   (3) Resposta p/q/r bem acima do piso de hover, SNR>=3 (regra do cap.9)
%   (4) Coerência entrada->saída gamma^2 >= 0.6 na banda  (regra CIFER)
%
% RCOU diferencial (QuadX: M1=FR, M2=RL, M3=FL, M4=RR):
%   roll  = (1+4)-(2+3)   pitch = (1+3)-(2+4)   yaw = (1+2)-(3+4)
%
% Uso:  >> scan_excitation_rc

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ---------------- CONFIG ----------------
LOG_FILE = 'logs_concat.mat';
dt   = 0.1;  fs = 1/dt;
BAND = [0.2 2.0];        % banda da malha de taxa (Hz)
COH_MIN = 0.6;           % limiar CIFER
% segmentos de treino (5) + validação (1), com o eixo-alvo de cada um
SEG = [4 24; 25 41; 42 62; 63 99; 100 125; 605 625];
TGT = {'yaw','yaw','roll','roll+pitch','pitch','VALID'};

%% ---------------- CARREGAR + GRADE ----------------
L = load_log_data(fullfile(paths.data, LOG_FILE));
t0 = max([L.time_IMU(1), L.time_RCOU(1), L.time_RCIN(1)]);
t1 = min([L.time_IMU(end), L.time_RCOU(end), L.time_RCIN(end)]);
tg = (t0:dt:t1)';
ip = @(t,x) interp1(t,x,tg,'linear','extrap');

% saída: taxas angulares (gyro)
p = ip(L.time_IMU,L.gyrX_raw); q = ip(L.time_IMU,L.gyrY_raw); r = ip(L.time_IMU,L.gyrZ_raw);

% RCOU -> entrada efetiva (diferencial de PWM)
w1=ip(L.time_RCOU,L.pwm1_raw); w2=ip(L.time_RCOU,L.pwm2_raw);
w3=ip(L.time_RCOU,L.pwm3_raw); w4=ip(L.time_RCOU,L.pwm4_raw);
wmean=(w1+w2+w3+w4)/4;
o_roll=(w1+w4)-(w2+w3); o_pitch=(w1+w3)-(w2+w4); o_yaw=(w1+w2)-(w3+w4);

% RCIN -> comando exógeno (desvio do neutro de cada canal)
i_roll = ip(L.time_RCIN,L.rcin_roll)  - median(L.rcin_roll);
i_pitch= ip(L.time_RCIN,L.rcin_pitch) - median(L.rcin_pitch);
i_yaw  = ip(L.time_RCIN,L.rcin_yaw)   - median(L.rcin_yaw);

% gaps entre logs concatenados -> invalida
gap = false(size(tg)); jp=find(diff(L.time_IMU)>1.0);
for j=jp', gap = gap | (tg>L.time_IMU(j) & tg<L.time_IMU(j+1)); end

%% ---------------- PISO DE HOVER (ruído) ----------------
% 10º percentil do std do gyro em janelas de 5 s VOANDO e válidas = hover quieto
wl=5; ws=(tg(1):wl:tg(end)-wl)';
fl=@(x) prctile(arrayfun(@(s) localstd(x,tg,s,wl,wmean,gap), ws), 10);
fp=fl(p); fq=fl(q); fr=fl(r); FL=[fp fq fr];

%% ---------------- RELATÓRIO ----------------
fprintf('================================================================================\n');
fprintf('  EXCITAÇÃO RCIN/RCOU por eixo — banda %.1f-%.1f Hz, CIFER g2>=%.1f, log %s\n',BAND,COH_MIN,LOG_FILE);
fprintf('================================================================================\n');
fprintf('  Piso de hover (ruído): p=%.3f  q=%.3f  r=%.3f rad/s\n\n', fp,fq,fr);

GY={p,q,r}; IN={i_roll,i_pitch,i_yaw}; OU={o_roll,o_pitch,o_yaw};

% --- Tabela 1: a excitação ESTÁ presente nas entradas? (std de RCIN e RCOU) ---
fprintf('  [1] AMPLITUDE DE EXCITAÇÃO  (std do desvio; us)\n');
fprintf('  %-12s %-11s |  RCIN  roll/pitch/yaw   |  RCOU  roll/pitch/yaw\n','seg','alvo');
fprintf('  %s\n',repmat('-',1,72));
for s=1:size(SEG,1)
  m=(tg>=SEG(s,1))&(tg<SEG(s,2));
  fprintf('  [%3d %3d]   %-11s |  %6.1f %6.1f %6.1f    |  %6.1f %6.1f %6.1f\n', ...
    SEG(s,1),SEG(s,2),TGT{s}, std(i_roll(m)),std(i_pitch(m)),std(i_yaw(m)), ...
    std(o_roll(m)),std(o_pitch(m)),std(o_yaw(m)));
end

% --- Tabela 2: a excitação é EFETIVA? (SNR da resposta + coerência) ---
fprintf('\n  [2] EFETIVIDADE  (SNR = std_gyro/piso ; g2 = coerência RCOU->gyro na banda)\n');
fprintf('  %-12s %-11s |   ROLL p     |   PITCH q    |   YAW r\n','seg','alvo');
fprintf('  %-12s %-11s | SNR   g2  %%b | SNR   g2  %%b | SNR   g2  %%b\n','','');
fprintf('  %s\n',repmat('-',1,74));
for s=1:size(SEG,1)
  m=(tg>=SEG(s,1))&(tg<SEG(s,2));
  vals=zeros(3,3);
  for ax=1:3
    gi=GY{ax}(m); oo=OU{ax}(m);
    snr=std(gi)/FL(ax);
    [f,c]=mscoh_(oo,gi,fs); bm=(f>=BAND(1))&(f<=BAND(2));
    vals(ax,:)=[snr mean(c(bm)) 100*mean(c(bm)>=COH_MIN)];
  end
  fprintf('  [%3d %3d]   %-11s | %4.1f %4.2f %3.0f | %4.1f %4.2f %3.0f | %4.1f %4.2f %3.0f\n', ...
    SEG(s,1),SEG(s,2),TGT{s}, vals(1,:),vals(2,:),vals(3,:));
end
fprintf('\n  Leitura: SNR>=3 e g2>=%.1f na banda => o eixo está bem excitado por aquela manobra.\n',COH_MIN);

%% ---------------- FIGURA: comando -> entrada -> resposta ----------------
% segmento dedicado de cada eixo (excitação mais limpa)
DED=[42 62; 100 125; 4 24]; nm={'ROLL  (p)','PITCH  (q)','YAW  (r)'};
DI={i_roll,i_pitch,i_yaw}; DO={o_roll,o_pitch,o_yaw}; DG={p,q,r}; col={'b','r',[0 .6 0]};
fig=figure('Color','w','Name','scan_excitation_rc','Position',[40 40 1400 780]);
for k=1:3
  a=DED(k,1); b=DED(k,2); m=(tg>=a)&(tg<b);
  subplot(3,1,k);
  yyaxis left
  plot(tg(m),DI{k}(m),'-','Color',[.6 .6 .6],'LineWidth',1.0); hold on;
  plot(tg(m),DO{k}(m),'-','Color',col{k},'LineWidth',1.1);
  ylabel('PWM dif. [\mus]');
  yyaxis right
  plot(tg(m),DG{k}(m),'k-','LineWidth',1.3);
  ylabel('taxa [rad/s]');
  grid on; title(sprintf('%s — segmento dedicado [%d, %d] s', nm{k}, a, b));
  if k==1, legend({'RCIN (comando exógeno)','RCOU diferencial (entrada da planta)','gyro (resposta)'}, ...
      'Location','northeast'); end
  if k==3, xlabel('t [s]'); end
end
saveas(fig, fullfile(paths.images,'scan_excitation_rc.png'));
fprintf('\n  Figura: %s\n', fullfile(paths.images,'scan_excitation_rc.png'));

%% ---------------- HELPERS ----------------
function v=localstd(x,tg,s,wl,wmean,gap)
  m=(tg>=s)&(tg<s+wl);
  if ~any(m) || any(gap(m)) || mean(wmean(m))<1400, v=NaN; return; end
  v=std(x(m));
end
function [f, coh] = mscoh_(x, y, fs)
% PSD/coerência de Welch (Hann, 50%% overlap) — sem Signal Toolbox
  x=x(:)-mean(x); y=y(:)-mean(y); N=numel(x);
  Ln=2^floor(log2(N/4)); if Ln<16, Ln=2^max(4,floor(log2(N/2))); end
  ov=floor(Ln/2); step=Ln-ov; win=0.5*(1-cos(2*pi*(0:Ln-1)'/(Ln-1)));
  nseg=max(1,floor((N-ov)/step)); Pxx=zeros(Ln,1);Pyy=zeros(Ln,1);Pxy=zeros(Ln,1);
  for i=1:nseg
    idx=(i-1)*step+(1:Ln); X=fft(x(idx).*win); Y=fft(y(idx).*win);
    Pxx=Pxx+X.*conj(X); Pyy=Pyy+Y.*conj(Y); Pxy=Pxy+X.*conj(Y);
  end
  coh=abs(Pxy).^2./(Pxx.*Pyy+eps); half=1:floor(Ln/2)+1;
  f=(0:floor(Ln/2))'*fs/Ln; coh=coh(half);
end
