% envelope_defesa.m — fundamentação do envelope de voo e do peso da aerodinâmica
% =========================================================================
% Figura única para responder à observação da banca sobre desprezar os efeitos
% aerodinâmicos. Quatro painéis:
%
%   (a) COMO a velocidade foi obtida, e sua aferição contra o GPS no trecho em
%       que o GPS ainda tinha solução. O método não usa GPS:
%         vertical    w  = V_D do EKF, auxiliado pelo barômetro
%         horizontal     = integração de a_NED = R·[0;0;−T/m] + g·e3, com atitude
%                          medida e empuxo da cadeia de motor, ancorada em v = 0
%                          nas bordas de sub-janelas de 10 s (tira a deriva)
%       A ancoragem faz o método NÃO servir como medida instantânea, e sim como
%       medida de ENVELOPE: o que ele reproduz bem são percentis e máximos.
%
%   (b) ENVELOPE propriamente dito: distribuição acumulada de V no treino e na
%       validação, com p95 e máximo marcados.
%
%   (c) FORÇA aerodinâmica contra o empuxo. Limite SUPERIOR conservador, tomando
%       o coeficiente de força total igual a 1, que é teto para esta fuselagem
%       em qualquer incidência:  F ≤ ½ρV²S·1.
%
%   (d) MOMENTO de arfagem: o de controle, o do amortecimento do modelo
%       identificado, e o da aerodinâmica da estrutura com C_mq MEDIDO no voo de
%       asa fixa de dez/2024 (não com o valor de painéis da AVL).
%
% Uso:  >> envelope_defesa
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();
rho = p.rho;  S = p.wing.S;  cbar = p.wing.c;  b = p.wing.b;  m = p.m;  g = p.g;

T_TRAINS = [4 24; 25 41; 42 62; 63 99; 100 125];
T_VALS   = [605 625; 610 630];
CMQ_VOO  = -5.540;      % medido em 4 doublets de profundor (11_fixed_wing)

L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU) min(L.time_ATT) min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU) max(L.time_ATT) max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';
VE = estimate_velocity(L, tg);
V  = VE.V;  V(~isfinite(V)) = 0;
ip = @(tt,xx) interp1(tt, xx, tg, 'linear');
W4 = [ip(L.time_RCOU,L.pwm1_raw) ip(L.time_RCOU,L.pwm2_raw) ip(L.time_RCOU,L.pwm3_raw) ip(L.time_RCOU,L.pwm4_raw)];
[rpm, fT] = motor_chain(tg, W4);
Ttot = sum(fT(rpm),2) * 1.06;                 % k_T identificado
q_g  = ip(L.time_IMU, L.gyrY_raw);

msk = @(WIN) any(cell2mat(arrayfun(@(k) tg>=WIN(k,1) & tg<=WIN(k,2), (1:size(WIN,1))', 'UniformOutput',false)'),2);
iT = msk(T_TRAINS);  iV = msk(T_VALS);

%% ---------------- aferição contra o GPS ----------------
gps_t = []; gps_v = [];
try
    raw = L.log_names{1};  if ~exist(raw,'file'), raw = fullfile(paths.data, raw); end
    G = load(raw, 'GPS_0').GPS_0;
    t0_log = min([min(L.time_IMU) min(L.time_ATT) min(L.time_RCOU)]);
    Lr = load_log_data(raw);
    shift = L.log_starts(1) - min([Lr.time_IMU(1) Lr.time_ATT(1) Lr.time_RCOU(1)]);
    ns = G(:,7);  ok = ns >= 5;
    gps_t = G(ok,2)/1e6 + shift;  gps_v = G(ok,12);      % col 12 = Spd [m/s]
catch ME
    fprintf('  (GPS indisponível: %s)\n', ME.message);
end

fprintf('\n  ================= ENVELOPE DE VOO =================\n');
fprintf('  %-12s %8s %8s %8s %8s\n','condição','média','p95','máx','amostras');
fprintf('  %-12s %8.2f %8.2f %8.2f %8d\n','treino',    mean(V(iT)), prctile(V(iT),95), max(V(iT)), sum(iT));
fprintf('  %-12s %8.2f %8.2f %8.2f %8d\n','validação', mean(V(iV)), prctile(V(iV),95), max(V(iV)), sum(iV));
if ~isempty(gps_t)
    jj = gps_t >= tg(1) & gps_t <= tg(end);
    Vg = gps_v(jj);  tgs = gps_t(jj);
    Ve_at = interp1(tg, V, tgs, 'linear');
    kk = isfinite(Ve_at) & tgs >= T_TRAINS(1,1) & tgs <= T_TRAINS(end,2);
    fprintf('\n  Aferição contra GPS (%d amostras, %.0f a %.0f s):\n', sum(kk), min(tgs(kk)), max(tgs(kk)));
    fprintf('    GPS       : média %.2f | p95 %.2f | máx %.2f m/s\n', mean(Vg(kk)), prctile(Vg(kk),95), max(Vg(kk)));
    fprintf('    estimada  : média %.2f | p95 %.2f | máx %.2f m/s\n', mean(Ve_at(kk)), prctile(Ve_at(kk),95), max(Ve_at(kk)));
    cc = corr_(Vg(kk), Ve_at(kk));
    fprintf('    correlação instantânea: %.2f   (o método é de ENVELOPE, não instantâneo)\n', cc);
end

%% ---------------- força e momento ----------------
F_max = 0.5*rho*V.^2*S*1.0;                 % teto conservador: C_F = 1
frac_F = 100*F_max./max(Ttot,1e-6);
kMq = 0.25*rho*S*cbar^2;
M_aero = kMq*abs(CMQ_VOO)*V.*abs(q_g);      % amortecimento da estrutura, C_mq MEDIDO
Mq_mod = 0.3372*abs(q_g);                   % amortecimento do modelo (M_q identificado)
fprintf('\n  ================= PESO DA AERODINÂMICA =================\n');
fprintf('  força / empuxo (teto C_F = 1):   treino p95 %.1f%% máx %.1f%%   |   validação p95 %.1f%% máx %.1f%%\n', ...
    prctile(frac_F(iT),95), max(frac_F(iT)), prctile(frac_F(iV),95), max(frac_F(iV)));
r_M = 100*M_aero./max(Mq_mod,1e-9);
fprintf('  amortecimento estrutura / modelo: treino mediana %.1f%%   |   validação mediana %.1f%%\n', ...
    median(r_M(iT & Mq_mod>0.01)), median(r_M(iV & Mq_mod>0.01)));

%% ================= FIGURA =================
f = figure('Position',[30 30 1250 940],'Color','w'); try, f.Theme='light'; catch, end
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% (a) aferição
nexttile; hold on; grid on;
if ~isempty(gps_t)
    jj = gps_t >= 63 & gps_t <= 99;
    plot(gps_t(jj), gps_v(jj), 'k-', 'LineWidth',1.5);
    ii = tg >= 63 & tg <= 99;
    plot(tg(ii), V(ii), '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.4);
    legend({'GPS (velocidade sobre o solo)','estimada (atitude + empuxo + barômetro)'}, ...
        'Location','northwest','FontSize',8);
    xlim([63 99]);
end
xlabel('t [s]'); ylabel('|V| [m/s]');
text(0, 1.05, '(a) aferição da estimativa no trecho com GPS (segmento de treino 63–99 s)', ...
    'Units','normalized','FontWeight','bold','FontSize',9);

% (b) envelope
nexttile; hold on; grid on;
[fT_,xT] = ecdf(V(iT));  [fV_,xV] = ecdf(V(iV));
plot(xT, 100*fT_, '-', 'Color',[0 0.45 0.7], 'LineWidth',2);
plot(xV, 100*fV_, '-', 'Color',[0.85 0.37 0.01], 'LineWidth',2);
yline(95,'k--','p95');
plot(prctile(V(iT),95), 95, 'o', 'MarkerFaceColor',[0 0.45 0.7], 'MarkerEdgeColor','k','MarkerSize',8);
plot(prctile(V(iV),95), 95, 'o', 'MarkerFaceColor',[0.85 0.37 0.01], 'MarkerEdgeColor','k','MarkerSize',8);
xlabel('|V| [m/s]'); ylabel('distribuição acumulada [%]'); xlim([0 3]);
legend({sprintf('treino (p95 %.2f, máx %.2f)', prctile(V(iT),95), max(V(iT))), ...
        sprintf('validação (p95 %.2f, máx %.2f)', prctile(V(iV),95), max(V(iV)))}, ...
        'Location','southeast','FontSize',8);
text(0, 1.05, '(b) envelope de velocidade dos ensaios', 'Units','normalized','FontWeight','bold','FontSize',9);

% (c) força: só as curvas de teto, com o envelope marcado
nexttile; hold on; grid on;
Vg_ = linspace(0, 3, 200)';
W_ = m*g;
plot(Vg_, 100*(0.5*rho*Vg_.^2*S*2.0)/W_, '-',  'Color',[0.20 0.20 0.20], 'LineWidth',2.2);
plot(Vg_, 100*(0.5*rho*Vg_.^2*S*1.0)/W_, '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1.8);
mk = [prctile(V(iT),95) 0 0.45 0.7; prctile(V(iV),95) 0.85 0.37 0.01; max(V) 0.4 0.4 0.4];
for k = 1:3
    vk = mk(k,1);  ck = mk(k,2:4);
    yk = 100*(0.5*rho*vk^2*S*2.0)/W_;
    plot([vk vk],[0 yk], '--', 'Color', ck, 'LineWidth',1.2);
    plot(vk, yk, 'o', 'MarkerFaceColor',ck, 'MarkerEdgeColor','k', 'MarkerSize',8);
    text(vk+0.04, yk, sprintf('%.1f%%', yk), 'FontSize',9, 'Color',ck, 'FontWeight','bold');
end
xlabel('|V| [m/s]'); ylabel('força aerodinâmica / peso [%]'); xlim([0 3]); ylim([0 16]);
legend({'teto estrito: C_F = 2 (placa plana de través)','referência: C_F = 1', ...
        sprintf('p95 treino %.2f m/s', mk(1,1)), '', ...
        sprintf('p95 validação %.2f m/s', mk(2,1)), '', ...
        sprintf('máximo %.2f m/s', mk(3,1))}, 'Location','northwest','FontSize',8);
text(0, 1.05, '(c) LIMITE SUPERIOR da força aerodinâmica, para qualquer incidência', ...
    'Units','normalized','FontWeight','bold','FontSize',9);

% (d) fração do amortecimento que a estrutura explica, em cada eixo
nexttile; hold on; grid on;
kLp = 0.25*rho*S*b^2;  kNr = kLp;
Lp_id = 0.2329;  Mq_id = 0.3372;  Nr_id = 0.1122;      % modelo oficial [N·m·s]
CLP = -0.0769;  CNR = -0.2037;                          % medidos em voo de asa fixa
sl = [kLp*abs(CLP)/Lp_id, kMq*abs(CMQ_VOO)/Mq_id, kNr*abs(CNR)/Nr_id];  % por m/s
cor = [0 0.45 0.7; 0.85 0.37 0.01; 0.45 0.70 0.35];
nomes = {'rolagem', 'arfagem', 'guinada'};
for k = 1:3
    plot(Vg_, 100*sl(k)*Vg_, '-', 'Color', cor(k,:), 'LineWidth', 2.2);
end
xregion(0, prctile(V(iT),95), 'FaceColor',[0 0.45 0.7], 'FaceAlpha',0.07);
xregion(prctile(V(iT),95), prctile(V(iV),95), 'FaceColor',[0.85 0.37 0.01], 'FaceAlpha',0.07);
for k = 1:3
    for vk = [prctile(V(iT),95), prctile(V(iV),95)]
        plot(vk, 100*sl(k)*vk, 'o', 'MarkerFaceColor',cor(k,:), 'MarkerEdgeColor','k','MarkerSize',6);
    end
    text(2.85, 100*sl(k)*2.85, sprintf(' %.1f%%/(m/s)', 100*sl(k)), 'Color',cor(k,:), 'FontSize',8.5, 'FontWeight','bold');
end
xlabel('|V| [m/s]'); ylabel('estrutura / amortecimento do modelo [%]'); xlim([0 3.3]); ylim([0 75]);
legend(nomes, 'Location','northwest','FontSize',9);
text(0, 1.05, '(d) fração do amortecimento explicada pela estrutura (C_{lp}, C_{mq}, C_{nr} medidos em voo)', ...
    'Units','normalized','FontWeight','bold','FontSize',9);

fn = fullfile(paths.images,'envelope_defesa.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',150);
fprintf('\n  Figura: %s\n', fn);

function c = corr_(a,b)
    a = a(:)-mean(a(:)); b = b(:)-mean(b(:));
    c = (a.'*b)/max(sqrt((a.'*a)*(b.'*b)),eps);
end
