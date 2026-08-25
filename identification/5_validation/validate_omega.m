% validate_omega.m — o amortecimento deve acompanhar a rotação dos rotores?
% =========================================================================
% Três variantes do MESMO vetor identificado (oficial_2026), sem parâmetro novo:
%   constante    L_p, M_q, N_r fixos (modelo atual)
%   ∝ Ω          escalam com ΣΩ_i(t)/ΣΩ_ref     (forma por rotor da NASA, n = 1)
%   ∝ T (∝ Ω²)   escalam com (ΣΩ_i/ΣΩ_ref)²     (k_v ∝ T, n = 2)
% Validação em modo full nos 12 trechos. Ω_ref = RPM de pairado da bancada.
% Uso:  >> validate_omega
% =========================================================================
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
WINS = [605 625; 610 630; 550 600; 400 420; 420 440; 440 460; 460 480; ...
        525 545; 545 565; 565 585; 585 605; 630 650];
P = load(fullfile(paths.outputs,'runs','oficial_2026','P_identified.mat')).P_final(:);
T0 = proj.m*proj.g/4;
rpm_ref = interp1(proj.bench.T_grams*1e-3*proj.g, proj.bench.RPM, T0, 'linear');

L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU) min(L.time_ATT) min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU) max(L.time_ATT) max(L.time_RCOU)]);
tg = (t_lo:0.1:t_hi)';
set_aero_vel(L, tg);
setappdata(0,'damp_form_override','moment');  clear aero_gains
gapA = L.boundaries(end);  gapB = L.log_starts(end);
constants = struct('m', proj.m, 'g', proj.g);
R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

C = struct('nome',{'constante','∝ Ω (NASA, n=1)','∝ T (n=2)'}, 'n',{0,1,2});
res = nan(size(WINS,1), numel(C), 3);  ratio_stats = nan(size(WINS,1),3);
fprintf('\n  Ω_ref (pairado, bancada) = %.0f rpm\n', rpm_ref);
fprintf('\n  %-11s | %-22s | %-22s | %-22s | ΣΩ/ΣΩref (min med max)\n', 'janela', C.nome);
for w = 1:size(WINS,1)
    tw = WINS(w,:);
    if tw(1) < gapB && tw(2) > gapA, continue; end
    idx = tg>=tw(1) & tg<=tw(2);  time = tg(idx);
    ip = @(tt,xx) interp1(tt, xx, time, 'linear');
    pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
    pqr = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
    att = [ip(L.time_ATT,L.roll_deg), ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
    rpm = motor_chain(time, pwm);  rt = sum(rpm,2)/(4*rpm_ref);
    ratio_stats(w,:) = [min(rt) median(rt) max(rt)];
    for k = 1:numel(C)
        if C(k).n == 0, if isappdata(0,'damp_omega'), rmappdata(0,'damp_omega'); end
        else, setappdata(0,'damp_omega', struct('n',C(k).n,'rpm_ref',rpm_ref)); end
        r = sim_window('full', P, time, pwm, pqr, att, constants);
        res(w,k,:) = [R2(pqr(:,1),r.p), R2(pqr(:,2),r.q), R2(pqr(:,3),r.r)];
    end
    fprintf('  %4.0f–%-6.0f', tw);
    for k = 1:numel(C), fprintf(' | %6.3f %6.3f %6.3f', squeeze(res(w,k,:))); end
    fprintf(' | %.2f %.2f %.2f\n', ratio_stats(w,:));
end
if isappdata(0,'damp_omega'), rmappdata(0,'damp_omega'); end
rmappdata(0,'damp_form_override');
ok = all(isfinite(res(:,1,1)),2);
fprintf('  %-11s', 'MÉDIA');
for k = 1:numel(C), fprintf(' | %6.3f %6.3f %6.3f', mean(res(ok,k,:),1)); end, fprintf('\n');
for k = 2:numel(C)
    d = squeeze(res(ok,k,:) - res(ok,1,:));  n = sum(ok);
    fprintf('  Δ vs constante (%s): p %+.4f±%.4f  q %+.4f±%.4f  r %+.4f±%.4f\n', C(k).nome, ...
        mean(d(:,1)), std(d(:,1))/sqrt(n), mean(d(:,2)), std(d(:,2))/sqrt(n), mean(d(:,3)), std(d(:,3))/sqrt(n));
end
save(fullfile(paths.outputs,'validate_omega.mat'), 'res','WINS','ratio_stats');
