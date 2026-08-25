function F = fw_load(which_log, varargin)
%FW_LOAD  Carrega um log de asa fixa de dez/2024 em formato uniforme.
%
%   F = fw_load('long')   voo longitudinal   (Log-DH-longitudinal-dez-24)
%   F = fw_load('lat')    voo látero-direcional
%   F = fw_load(..., 'dt', 0.04)   passo da grade (default 0,04 s = 25 Hz)
%
% Os dois voos foram feitos em MANUAL (ModeNum 0), ou seja, comando do piloto
% direto nas superfícies, sem estabilização. É malha aberta, que é a condição
% boa para identificar. Mapa dos canais, lido do PARM do próprio log:
%   SERVO1 função 4  = aileron     SERVO2 função 19 = profundor
%   SERVO3 função 70 = motor       SERVO4 função 21 = leme
%   todos com MIN 1100, TRIM 1500, MAX 1900, não revertidos
%
% Saída F, tudo na mesma grade de tempo F.t:
%   .V        velocidade verdadeira do pitot (ARSP.Airspeed) [m/s]
%   .alpha    ângulo de ataque estimado pelo EKF (AOA.AOA)   [rad]
%   .beta     ângulo de derrapagem estimado (AOA.SSA)        [rad]
%   .p .q .r  taxas do giroscópio                            [rad/s]
%   .phi .theta .psi   atitude                               [rad]
%   .da .de .dr        deflexões NORMALIZADAS em [-1, 1]     (PWM - 1500)/400
%   .thr      comando de motor normalizado [0, 1]
%   .ax .ay .az        acelerômetro                          [m/s²]
%   .VN .VE .VD        velocidade do EKF                     [m/s]
%   .alt      altitude barométrica [m]
%
% A deflexão FÍSICA em radianos é da·δ_max, com δ_max o curso mecânico da
% superfície, que não está no log. Por isso as deflexões saem normalizadas: nas
% identificações o produto C_ℓδa·δ_max entra como UM parâmetro só.

    opt = struct('dt', 0.04);
    for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
    paths = setup_paths();
    base = fullfile(fileparts(paths.root), 'data', 'DH_modo_AsaFixa_voo_dez_2024');
    switch lower(which_log)
        case {'long','longitudinal'}, fn = 'Log-DH-longitudinal-dez-24.bin-440203.mat';
        case {'lat','laterodirecional','latdir'}, fn = 'Log-DH-laterodirecional-dez-24.bin-449987.mat';
        otherwise, error('fw_load:which', 'use ''long'' ou ''lat''');
    end
    fp = fullfile(base, fn);
    S = load(fp, 'ARSP','AOA','ATT','IMU','RCOU','XKF1','BARO','AETR');

    tt = @(X) X(:,2)/1e6;
    msgs = {S.IMU, S.ATT, S.RCOU, S.ARSP, S.AOA, S.XKF1, S.BARO};
    lo = -Inf; hi = Inf;
    for k = 1:numel(msgs)
        tk = tt(msgs{k});  lo = max(lo, tk(1));  hi = min(hi, tk(end));
    end
    t = (lo:opt.dt:hi)';
    ip = @(X, col) interp_(tt(X), X, col, t);

    F.file = fp;  F.t = t;  F.dt = opt.dt;  F.name = fn;
    F.V     = ip(S.ARSP, 4);                      % Airspeed
    F.alpha = deg2rad(ip(S.AOA, 3));              % AOA
    F.beta  = deg2rad(ip(S.AOA, 4));              % SSA
    F.p = ip(S.IMU, 4);  F.q = ip(S.IMU, 5);  F.r = ip(S.IMU, 6);
    F.ax = ip(S.IMU, 7); F.ay = ip(S.IMU, 8); F.az = ip(S.IMU, 9);
    F.phi = deg2rad(ip(S.ATT, 4));  F.theta = deg2rad(ip(S.ATT, 6));  F.psi = deg2rad(ip(S.ATT, 8));
    F.da = (ip(S.RCOU, 3) - 1500)/400;            % C1 aileron
    F.de = (ip(S.RCOU, 4) - 1500)/400;            % C2 profundor
    F.thr = (ip(S.RCOU, 5) - 1100)/800;           % C3 motor
    F.dr = (ip(S.RCOU, 6) - 1500)/400;            % C4 leme
    F.VN = ip(S.XKF1, 7); F.VE = ip(S.XKF1, 8); F.VD = ip(S.XKF1, 9);
    F.alt = ip(S.BARO, 4);
end

function y = interp_(tx, X, col, t)
    [tu, iu] = unique(tx);
    y = interp1(tu, X(iu, col), t, 'linear');
end
