%% analise_damping.m
%  Fecha quantitativamente a pergunta: "Dp, Dq, Dr são física ou fudge?"
%  (1) ABLAÇÃO: roda a sim com D e com D=0 -> os D são "load-bearing"?
%  (2) IDENTIFICABILIDADE: correlação das sensibilidades + CRB -> Dp troca
%      com Jx/kT (=fudge) ou tem assinatura própria (=física)?
%
%  Usa o modelo do usuário (sim_window/vtol_dynamics) SEM editá-lo.
%  Janela [62,98]s (a mesma do plot validate_params 'full').
%  Rodar:  >> analise_damping

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();

% ---- P_MANUAL REAL do usuário ----
P = [0.055; 0.095; 0.253; 0.000933; ...   % Jx Jy Jz Jxz
     0.76; 0.76; 0.63; 0.64; ...           % kT1..4
     1.16; 1.15; 1.30; 1.28; ...           % kQ1..4
     3.57; 1.83; 0.56];                    % Dp Dq Dr
pnames = {'Jx','Jy','Jz','Jxz','kT1','kT2','kT3','kT4','kQ1','kQ2','kQ3','kQ4','Dp','Dq','Dr'};

W = [1 30];

%% ---- carrega log + grade 0.1s (igual validate_params) ----
proj = parameters();
constants = struct('m', proj.m, 'g', proj.g);
L = load_log_data(fullfile(setup_paths().data, 'logs_concat.mat'));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)'; sel = tg>=W(1) & tg<=W(2); time = tg(sel); N = numel(time);
pwm = [interp1(L.time_RCOU,L.pwm1_raw,time), interp1(L.time_RCOU,L.pwm2_raw,time), ...
       interp1(L.time_RCOU,L.pwm3_raw,time), interp1(L.time_RCOU,L.pwm4_raw,time)];
pqr_meas = [interp1(L.time_IMU,L.gyrX_raw,time), interp1(L.time_IMU,L.gyrY_raw,time), ...
            interp1(L.time_IMU,L.gyrZ_raw,time)];
att_meas = [interp1(L.time_ATT,L.roll_deg,time), interp1(L.time_ATT,L.pitch_deg,time), ...
            interp1(L.time_ATT,L.yaw_deg,time)];
R2 = @(ye,ys) 1 - sum((ye-ys).^2)/max(sum((ye-mean(ye)).^2),1e-12);
fprintf('Janela [%g,%g]s | %d amostras\n', W(1),W(2),N);

% rates independem da atitude no modelo -> 'hybrid' dá p,q,r idênticos ao 'full'
sim_fn = @(Pv) sim_window('hybrid', Pv, time, pwm, pqr_meas, att_meas, constants);

%% ===== 1. ABLAÇÃO D=0 =====
r_nom = sim_fn(P);
P0 = P; P0(13:15) = 0;
r_abl = sim_fn(P0);
fprintf('\n===== (1) ABLAÇÃO: com D  vs  D=0 =====\n');
fprintf('  canal |  R²(com D) |  R²(D=0)  | RMS(com D) | RMS(D=0)  [rad/s]\n');
ch = {'p','q','r'};
for i = 1:3
    yn = r_nom.(ch{i}); ya = r_abl.(ch{i}); ym = pqr_meas(:,i);
    fprintf('    %s   |  %+8.3f  | %+8.3f  |  %7.3f   |  %7.3f\n', ...
        ch{i}, R2(ym,yn), R2(ym,ya), rms(ym-yn), rms(ym-ya));
end

%% ===== 2. IDENTIFICABILIDADE (Jacobiano por diferenças centrais) =====
wstd = std(pqr_meas,0,1);                          % escala por canal
stack = @(r) [r.p/wstd(1); r.q/wstd(2); r.r/wstd(3)];
nP = numel(P);  Jac = zeros(3*N, nP);  e = 1e-3;
for i = 1:nP
    d = e*max(abs(P(i)), 1e-6);
    Pp = P; Pp(i) = P(i)+d;  Pm = P; Pm(i) = P(i)-d;
    Jac(:,i) = (stack(sim_fn(Pp)) - stack(sim_fn(Pm)))/(2*d) * P(i);  % ∂/∂lnP
end
nrm  = sqrt(sum(Jac.^2,1));                          % observabilidade de cada param
Rs   = (Jac'*Jac) ./ max(nrm'*nrm, 1e-30);           % correlação das sensibilidades
M    = Jac'*Jac;

fprintf('\n===== (2) IDENTIFICABILIDADE =====\n');
fprintf('cond(M) = %.2e\n', cond(M));
fprintf('\nObservabilidade |∂y/∂lnP| (maior = mais "visível" nos dados):\n');
fprintf('  Dp=%.2f  Dq=%.2f  Dr=%.2f  |  Jx=%.2f Jy=%.2f Jz=%.2f Jxz=%.4f\n', ...
    nrm(13),nrm(14),nrm(15), nrm(1),nrm(2),nrm(3),nrm(4));

fprintf('\nCorrelação de sensibilidade dos D (|.|~1 => confundido = fudge):\n');
fprintf('        vs Jx   vs Jy   vs Jz   vs Jxz  vs kT(max) vs kQ(max)\n');
for d = [13 14 15]
    fprintf('  %-3s  %+5.2f   %+5.2f   %+5.2f   %+5.2f    %5.2f      %5.2f\n', ...
        pnames{d}, Rs(d,1),Rs(d,2),Rs(d,3),Rs(d,4), ...
        max(abs(Rs(d,5:8))), max(abs(Rs(d,9:12))));
end

% CRB (incerteza relativa) — só se M razoavelmente condicionada
if cond(M) < 1e12
    res = stack(r_nom) - [pqr_meas(:,1)/wstd(1); pqr_meas(:,2)/wstd(2); pqr_meas(:,3)/wstd(3)];
    s2  = (res'*res)/(3*N - nP);
    C   = s2 * pinv(M);
    crb = sqrt(diag(C));                              % incerteza relativa (∂lnP)
    fprintf('\nCRB (incerteza relativa 1σ; <~0.3 = bem identificado):\n');
    fprintf('  Dp=%.2f (%.0f%%)  Dq=%.2f (%.0f%%)  Dr=%.2f (%.0f%%)\n', ...
        crb(13),100*crb(13), crb(14),100*crb(14), crb(15),100*crb(15));
end

%% ---- figura: ablação visual ----
fig = figure('Color','w','Position',[60 60 1300 380]);
for i = 1:3
    subplot(1,3,i); hold on; grid on;
    plot(time, pqr_meas(:,i), 'b-', 'LineWidth',1.2, 'DisplayName','medido');
    plot(time, r_nom.(ch{i}), 'r--','LineWidth',1.2, 'DisplayName','sim com D');
    plot(time, r_abl.(ch{i}), 'Color',[1 .6 0], 'LineStyle',':', 'LineWidth',1.4, 'DisplayName','sim D=0');
    title(ch{i}); xlabel('t [s]'); ylabel('rad/s');
    if i==1, legend('Location','best'); end
end
sgtitle('Ablação do amortecimento: medido vs sim(com D) vs sim(D=0)');
saveas(fig, fullfile(fileparts(mfilename('fullpath')), 'analise_damping.png'));
fprintf('\nFigura: analise_damping.png\n');
