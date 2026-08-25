% envelope_velocity.m — Envelope de velocidade da campanha e ordem de grandeza das
%                       forças aerodinâmicas (o que a hipótese "baixa velocidade" vale)
% =========================================================================
% O modelo despreza as forças aerodinâmicas da estrutura (asa, fuselagem) e o
% arrasto dos rotores no modo multirrotor. Para sustentar isso é preciso (1) mostrar
% em que velocidades a aeronave voou e (2) quanto essas forças valeriam lá.
%
% FONTES DE VELOCIDADE
%   • Voo 4 (log 1 do concat, 0–268 s): GPS com fix (5 Hz) → velocidade sobre o
%     solo MEDIDA. É onde estão as janelas de TREINO.
%   • Voo 7 (517–749 s, janelas de VALIDAÇÃO): GPS sem fix (0 satélites) e EKF sem
%     auxílio de velocidade (XKF4.SV = 0) → NÃO há medida de velocidade. Usa-se
%     como proxy a atitude (|φ|,|θ|) e a atividade de PWM, que são da mesma classe
%     das janelas de treino (mesmo local, mesmo piloto, mesmas manobras).
%
% FORÇAS AVALIADAS (em % do peso, mg = 19,5 N)
%   placa plana (pior caso, C=2):   F = q̄·S·2
%   sustentação da asa a α=0:       F = q̄·S·CL0   (USA-35B, CL0 ≈ 0,36)
%   arrasto de rotor (força H):     F = 4·k_h·V   (prior_damping.m, teoria de pá)
%   arrasto induzido (Beard, C_d=0,1, limite alto da literatura): F = m·g·C_d·V
%
% Uso:  >> envelope_velocity
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();  g = proj.g;  m = proj.m;  mg = m*g;
S = proj.wing.S;  rho = 1.225;  CL0 = 0.36;
try, pr = load(fullfile(paths.outputs,'prior_damping.mat')).prior; k_h = pr.k_h; catch, k_h = 0.0156; end

L  = load_log_data(fullfile(paths.data,'logs_concat.mat'));
dt = 0.1;
T_TRAIN = [4 24; 25 41; 42 62; 63 99; 100 125];
T_VAL   = [610 630];
T_VAL2  = [550 600];

%% ---------- Voo 4: GPS ----------
tG = L.time_GPS(:);  spd = L.gps_spd(:);  vd = L.gps_vd(:);
inTrain = false(size(tG)); for i=1:size(T_TRAIN,1), inTrain = inTrain | (tG>=T_TRAIN(i,1) & tG<=T_TRAIN(i,2)); end
inFlight = tG>=4 & tG<=268 & spd>0;
st = @(x) [median(x), prctile(x,95), max(x)];
sT = st(spd(inTrain));  sF = st(spd(inFlight));  vT = st(abs(vd(inTrain)));
fprintf('\n=== Voo 4 (GPS): velocidade sobre o solo ===\n');
fprintf('  janelas de TREINO (%d s): mediana %.2f | p95 %.2f | máx %.2f m/s   (vertical |Vd| p95 %.2f)\n', ...
    sum(diff(T_TRAIN,1,2)), sT, vT(2));
fprintf('  voo inteiro (4–268 s):     mediana %.2f | p95 %.2f | máx %.2f m/s\n', sF);

%% ---------- Atitude e PWM: treino × validação (mesma classe de manobra?) ----------
tg = (max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]):dt:min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]))';
ip = @(tt,xx) interp1(tt,xx,tg,'linear');
tilt = hypot(ip(L.time_ATT,L.roll_deg), ip(L.time_ATT,L.pitch_deg));         % inclinação total [°]
W = [ip(L.time_RCOU,L.pwm1_raw) ip(L.time_RCOU,L.pwm2_raw) ip(L.time_RCOU,L.pwm3_raw) ip(L.time_RCOU,L.pwm4_raw)];
pwm_spread = max(W,[],2) - min(W,[],2);                                       % atividade de comando [µs]
mTr = false(size(tg)); for i=1:size(T_TRAIN,1), mTr = mTr | (tg>=T_TRAIN(i,1) & tg<=T_TRAIN(i,2)); end
mVa = tg>=T_VAL(1) & tg<=T_VAL(2);   mVb = tg>=T_VAL2(1) & tg<=T_VAL2(2);
fprintf('\n=== Proxies de envelope (atitude e comando) ===\n');
fprintf('  %-22s | tilt p50 %5s p95 %5s máx %5s [°] | ΔPWM p95 %6s [µs]\n','trecho','','','','');
prx = @(nm,mm) fprintf('  %-22s |      %5.1f     %5.1f     %5.1f     |         %6.0f\n', nm, median(tilt(mm)), prctile(tilt(mm),95), max(tilt(mm)), prctile(pwm_spread(mm),95));
prx('treino (voo 4)', mTr);  prx('validação 610–630 (v7)', mVa);  prx('validação 550–600 (v7)', mVb);
% aceleração horizontal quase-estática ≈ g·tan(tilt): dá ideia da velocidade que a manobra gera
fprintf('  aceleração horizontal quase-estática g·tan(tilt p95): treino %.2f | val %.2f m/s²\n', g*tand(prctile(tilt(mTr),95)), g*tand(prctile(tilt(mVa),95)));

%% ---------- Ordem de grandeza das forças ----------
V = linspace(0,4,401);  qbar = 0.5*rho*V.^2;
F_plate = qbar*S*2;         F_lift = qbar*S*CL0;
F_H     = 4*k_h*V;          F_beard = mg*0.1*V;
pct = @(F) 100*F/mg;
Vref = [sT(2), sT(3), sF(3), 4];  Vlab = {'p95 treino','máx treino','máx voo','4 m/s'};
fprintf('\n=== Forças aerodinâmicas em %% do peso (mg = %.1f N) ===\n', mg);
fprintf('  %-12s %6s | %10s %10s %10s %12s\n','ponto','V','placa(C=2)','asa CL0','rotor H','Beard Cd=0,1');
for i=1:numel(Vref)
    v = Vref(i); qb = 0.5*rho*v^2;
    fprintf('  %-12s %6.2f | %9.1f%% %9.1f%% %9.1f%% %11.1f%%\n', Vlab{i}, v, pct(qb*S*2), pct(qb*S*CL0), pct(4*k_h*v), pct(mg*0.1*v));
end
fprintf('  referência de comando: momento de controle ≈ 1 N·m ↔ força ≈ 4 N (20%% de mg) no braço; ΔT de manobra ≈ 2 N (10%%)\n');

%% ---------- FIGURA ----------
f = figure('Position',[40 40 1150 900],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% (a) GPS voo 4
nexttile; hold on; grid on;
yl = [0 3.2];
for i=1:size(T_TRAIN,1), patch([T_TRAIN(i,1) T_TRAIN(i,2) T_TRAIN(i,2) T_TRAIN(i,1)],[yl(1) yl(1) yl(2) yl(2)],[0.85 0.93 0.85],'EdgeColor','none'); end
plot(tG(inFlight), spd(inFlight), 'k-', 'LineWidth',1.2);
yline(sT(2),'--','Color',[0.85 0.37 0.01],'LineWidth',1.4,'Label',sprintf('p95 treino = %.2f m/s',sT(2)),'LabelHorizontalAlignment','left');
yline(sT(3),':','Color',[0.85 0.37 0.01],'LineWidth',1.2,'Label',sprintf('máx treino = %.2f m/s',sT(3)),'LabelHorizontalAlignment','left');
ylim(yl); xlim([0 268]); ylabel('|V| GPS [m/s]'); text(0.005,1.06,'(a) Voo 4 (treino): velocidade sobre o solo medida pelo GPS; janelas de treino em verde','Units','normalized','FontWeight','bold');
% (b) proxies: tilt treino × validação
nexttile; hold on; grid on;
edges = 0:1:30;
histogram(tilt(mTr), edges, 'Normalization','probability', 'FaceColor',[0.3 0.6 0.3], 'FaceAlpha',0.6);
histogram(tilt(mVa), edges, 'Normalization','probability', 'FaceColor',[0 0.45 0.7], 'FaceAlpha',0.6);
histogram(tilt(mVb), edges, 'Normalization','probability', 'FaceColor',[0.85 0.37 0.01], 'FaceAlpha',0.4);
xlabel('inclinação \surd(\phi^2+\theta^2) [°]'); ylabel('fração do tempo');
legend({'treino (voo 4, GPS)','validação 610–630 (voo 7)','validação 550–600 (voo 7)'},'Location','northeast');
text(0.005,1.06,'(b) Voo 7 não tem velocidade medida: a atitude mostra que as manobras são da mesma classe do treino','Units','normalized','FontWeight','bold');
% (c) forças
nexttile; hold on; grid on;
plot(V, pct(F_plate), 'k-', 'LineWidth',1.6); plot(V, pct(F_lift), 'k--', 'LineWidth',1.4);
plot(V, pct(F_beard), '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.4); plot(V, pct(F_H), '-', 'Color',[0 0.45 0.7], 'LineWidth',1.4);
yline(5,':','Color',[0.4 0.4 0.4],'Label','5% de mg'); yline(1,':','Color',[0.4 0.4 0.4],'Label','1% de mg');
xline(sT(2),'--','Color',[0.3 0.6 0.3],'LineWidth',1.4,'Label','p95 treino','LabelVerticalAlignment','bottom');
xline(sT(3),':','Color',[0.3 0.6 0.3],'LineWidth',1.2,'Label','máx treino','LabelVerticalAlignment','bottom');
ylim([0 25]); xlabel('velocidade [m/s]'); ylabel('força / mg [%]');
legend({'placa plana C=2 (pior caso)','sustentação da asa a \alpha=0','arrasto induzido Beard, C_d=0,1','força H dos rotores (teoria)'},'Location','northwest');
text(0.005,1.06,'(c) Forças aerodinâmicas não modeladas em função da velocidade (% do peso)','Units','normalized','FontWeight','bold');
exportgraphics(f, fullfile(paths.images,'envelope_velocity.png'), 'BackgroundColor','white','Resolution',150);
fprintf('\n  Figura: %s\n', fullfile(paths.images,'envelope_velocity.png'));
