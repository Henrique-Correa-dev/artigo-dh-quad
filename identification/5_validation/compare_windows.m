% compare_windows.m — Θ_OEM da dissertação × Θ_OEM do programa em VÁRIAS janelas
% =========================================================================
% Complemento do compare_dissertacao.m: só simula (não otimiza) os dois vetores
% OEM em uma lista de janelas de validação e nos dois modos de simulação
% (hybrid = atitude medida, isola o modelo rotacional; full = integra tudo,
% como na Tabela 5.5). Objetivo: ver se a conclusão "os dois ajustam igual"
% depende da janela.
%
% Uso:  >> compare_windows
% =========================================================================
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

LOG_FILE  = 'logs_concat.mat';
T_WINDOWS = [605 625; 610 630];       % linhas = janelas [t0 t1]
MODES     = {'hybrid','full'};
DELAY_PWM = 0;   % atraso já dentro de sim_window (motor_chain)

kT_ref = 1.764e-7;  kQ_ref = 5.04e-9;
sc = @(cT, cM) [cT(:)*1e-7/kT_ref; cM(:)*1e-9/kQ_ref];
proj = parameters();
prog = load(fullfile(paths.outputs,'P_identified.mat'));
S(1).tag = 'diss_OEM';  S(1).P = [0.0469; 0.0984; 0.1453; 0.00261; sc([1.87 1.86 1.86 1.89],[3.38 3.16 3.92 3.47]); 6.85; 4.68; 0.81];
S(2).tag = 'prog_OEM';  S(2).P = prog.P_final(:);

L = load_log_data(fullfile(paths.data, LOG_FILE));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg   = (t_lo:0.1:t_hi)';
constants = struct('m', proj.m, 'g', proj.g);
sig = {'p','q','r','aX','aY','aZ'};

for w = 1:size(T_WINDOWS,1)
    tw = T_WINDOWS(w,:);
    idx  = (tg >= tw(1)) & (tg <= tw(2));  time = tg(idx);
    ip   = @(tt,xx) interp1(tt, xx, time, 'linear');
    pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
    if DELAY_PWM > 0
        dl = @(x) [repmat(x(1),DELAY_PWM,1); x(1:end-DELAY_PWM)];
        for ci = 1:4, pwm(:,ci) = dl(pwm(:,ci)); end
    end
    pqr_meas = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
    acc_meas = [ip(L.time_IMU,L.accX_raw), ip(L.time_IMU,L.accY_raw), ip(L.time_IMU,L.accZ_raw)];
    att_meas = [ip(L.time_ATT,L.roll_deg),  ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
    ms = [pqr_meas, acc_meas];

    for mi = 1:numel(MODES)
        mode = MODES{mi};
        fprintf('\n=== Janela %g–%g s | modo %s | %d amostras ===\n', tw, mode, numel(time));
        fprintf('%-9s | %6s %6s %6s %6s %6s %6s | %6s %6s %6s %6s %6s %6s\n', 'conjunto', ...
            'R2 p','R2 q','R2 r','R2 aX','R2 aY','R2 aZ','TIC p','TIC q','TIC r','TIC aX','TIC aY','TIC aZ');
        for k = 1:numel(S)
            r = sim_window(mode, S(k).P, time, pwm, pqr_meas, att_meas, constants);
            sm = [r.p, r.q, r.r, r.accX, r.accY, r.accZ];
            R2 = zeros(1,6); TIC = zeros(1,6);
            for i = 1:6, m = fit_metrics(ms(:,i), sm(:,i)); R2(i) = m.R2; TIC(i) = m.TIC; end
            fprintf('%-9s | %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n', S(k).tag, R2, TIC);
        end
    end
end
