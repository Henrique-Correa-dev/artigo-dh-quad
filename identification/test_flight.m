%% test_flight.m
%  Junta múltiplos logs de voo num único time-base contínuo e plota:
%  RCOU (motores), pqr (gyro), Acc (IMU), ATT (atitude), RCIN (rádio),
%  MAG (magnetômetro) e GPS (velocidade NED).
%  MAG e GPS são incluídos para alimentar um EKF de atitude (8_atitude_ekf_test/):
%  MAG corrige o yaw; velocidade do GPS corrige o tilt — o que falta a uma
%  estimativa só de IMU.
%
%  Cada log tem tempo absoluto independente (segundos desde boot da Pixhawk).
%  Aqui, cada log é shiftado pra começar imediatamente após o anterior + GAP.

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '1_data'));

% ╔══════════════════════════════════════════════════════════════════╗
% ║  CONFIGURAÇÃO                                                    ║
% ╚══════════════════════════════════════════════════════════════════╝
LOG_DIR = fullfile(fileparts(mfilename('fullpath')), '1_data');

LOGS = { ...
    '4 25-05-2026 09-31-48.log-132954.mat', ...
    '6 31-12-1979 21-00-00.bin-116760.mat', ...
    '7 31-12-1979 21-00-00.bin-113223.mat'};

GAP = 5;   % segundos de "lacuna" entre logs (separação visual)

%% ====================== Junção dos logs ============================
J = struct();   % junto
J.time_IMU  = []; J.gyrX = []; J.gyrY = []; J.gyrZ = [];
J.accX = []; J.accY = []; J.accZ = [];
J.time_ATT  = []; J.roll = []; J.pitch = []; J.yaw = [];
J.time_RCOU = []; J.pwm1 = []; J.pwm2 = []; J.pwm3 = []; J.pwm4 = [];
J.time_RCIN = []; J.rc_roll = []; J.rc_pitch = []; J.rc_thr = []; J.rc_yaw = [];
J.time_MAG = []; J.magX = []; J.magY = []; J.magZ = [];
J.time_GPS = []; J.gps_vn = []; J.gps_ve = []; J.gps_vd = [];
J.gps_spd = []; J.gps_gcrs = []; J.gps_nsats = [];

boundaries = zeros(numel(LOGS)-1, 1);   % tempo (no eixo unificado) onde cada log termina
log_starts = zeros(numel(LOGS), 1);     % onde cada log começa no eixo unificado

t_cursor = 0;
for i = 1:numel(LOGS)
    L = load_log_data(fullfile(LOG_DIR, LOGS{i}));

    % Normaliza: subtrai t_start do próprio log e adiciona o cursor
    t0_log = min([L.time_IMU(1), L.time_ATT(1), L.time_RCOU(1)]);
    t_shift = t_cursor - t0_log;

    log_starts(i) = t_cursor;

    J.time_IMU  = [J.time_IMU;  L.time_IMU  + t_shift];
    J.gyrX = [J.gyrX; L.gyrX_raw]; J.gyrY = [J.gyrY; L.gyrY_raw]; J.gyrZ = [J.gyrZ; L.gyrZ_raw];
    J.accX = [J.accX; L.accX_raw]; J.accY = [J.accY; L.accY_raw]; J.accZ = [J.accZ; L.accZ_raw];

    J.time_ATT  = [J.time_ATT;  L.time_ATT  + t_shift];
    J.roll  = [J.roll;  L.roll_deg];
    J.pitch = [J.pitch; L.pitch_deg];
    J.yaw   = [J.yaw;   L.yaw_deg];

    J.time_RCOU = [J.time_RCOU; L.time_RCOU + t_shift];
    J.pwm1 = [J.pwm1; L.pwm1_raw]; J.pwm2 = [J.pwm2; L.pwm2_raw];
    J.pwm3 = [J.pwm3; L.pwm3_raw]; J.pwm4 = [J.pwm4; L.pwm4_raw];

    if ~isempty(L.time_RCIN)
        J.time_RCIN = [J.time_RCIN; L.time_RCIN + t_shift];
        J.rc_roll  = [J.rc_roll;  L.rcin_roll];
        J.rc_pitch = [J.rc_pitch; L.rcin_pitch];
        J.rc_thr   = [J.rc_thr;   L.rcin_throttle];
        J.rc_yaw   = [J.rc_yaw;   L.rcin_yaw];
    end

    % MAG (magnetômetro) — lê o compass ATIVO direto do bruto.
    % O load_log_data dá MAG_0 (compass 1), MAS neste drone o compass 1 está
    % desligado (COMPASS_USE=0) e mal calibrado (|mag|~968 mGauss); o EKF usa o
    % compass 2 = MAG_1 (externo, |mag|~374 mGauss, casa melhor com o ATT).
    % read_mag_ prefere MAG_1 (fallback MAG_0 / load_log_data).
    Mg = read_mag_(fullfile(LOG_DIR, LOGS{i}));
    if isempty(Mg.time) && isfield(L,'time_MAG') && ~isempty(L.time_MAG)
        Mg.time = L.time_MAG; Mg.x = L.magX; Mg.y = L.magY; Mg.z = L.magZ;
    end
    if ~isempty(Mg.time)
        J.time_MAG = [J.time_MAG; Mg.time + t_shift];
        J.magX = [J.magX; Mg.x]; J.magY = [J.magY; Mg.y]; J.magZ = [J.magZ; Mg.z];
    end

    % GPS (velocidade NED) — lido direto do arquivo bruto; só fixes 3D válidos.
    % A posição absoluta NÃO é concatenada (cada log tem origem distinta), mas a
    % velocidade é relativa ao frame NED local, então é concatenável.
    Gp = read_gps_(fullfile(LOG_DIR, LOGS{i}));
    if ~isempty(Gp.time)
        J.time_GPS  = [J.time_GPS;  Gp.time + t_shift];
        J.gps_vn    = [J.gps_vn;    Gp.vn];
        J.gps_ve    = [J.gps_ve;    Gp.ve];
        J.gps_vd    = [J.gps_vd;    Gp.vd];
        J.gps_spd   = [J.gps_spd;   Gp.spd];
        J.gps_gcrs  = [J.gps_gcrs;  Gp.gcrs];
        J.gps_nsats = [J.gps_nsats; Gp.nsats];
    end

    % Avança cursor pro fim deste log + gap
    t_end_log = max([L.time_IMU(end), L.time_ATT(end), L.time_RCOU(end)]);
    t_cursor = t_end_log + t_shift + GAP;

    if i < numel(LOGS)
        boundaries(i) = t_end_log + t_shift;
    end

    fprintf('Log %d → janela unificada [%.1f, %.1f] s  (dur %.1f s)\n', ...
        i, log_starts(i), t_end_log + t_shift, t_end_log + t_shift - log_starts(i));
end

t_total = t_cursor - GAP;
fprintf('\nTotal: %.1f s de dados concatenados de %d logs.\n', t_total, numel(LOGS));

%% ====================== Salva como .mat reutilizável =================
% Formato compatível com a saída de load_log_data() — pra que qualquer
% script (identify_plant, validate_params, etc.) possa fazer:
%   tmp = load(LOG_FILE);  L = tmp.L;
% no lugar de:
%   L = load_log_data(LOG_FILE);
L = struct();
L.format = 'concat';
L.source = strjoin(LOGS, ' + ');
L.time_IMU = J.time_IMU;
L.gyrX_raw = J.gyrX; L.gyrY_raw = J.gyrY; L.gyrZ_raw = J.gyrZ;
L.accX_raw = J.accX; L.accY_raw = J.accY; L.accZ_raw = J.accZ;
L.time_ATT = J.time_ATT;
L.roll_deg = J.roll; L.pitch_deg = J.pitch; L.yaw_deg = J.yaw;
L.time_RCOU = J.time_RCOU;
L.pwm1_raw = J.pwm1; L.pwm2_raw = J.pwm2;
L.pwm3_raw = J.pwm3; L.pwm4_raw = J.pwm4;
L.time_RCIN = J.time_RCIN;
L.rcin_roll = J.rc_roll; L.rcin_pitch = J.rc_pitch;
L.rcin_throttle = J.rc_thr; L.rcin_yaw = J.rc_yaw;
% Magnetômetro (yaw) e GPS-velocidade (tilt) — para o EKF de atitude
L.time_MAG = J.time_MAG;
L.magX = J.magX; L.magY = J.magY; L.magZ = J.magZ;
L.time_GPS = J.time_GPS;    % só velocidade (posição absoluta não é concatenável)
L.gps_vn = J.gps_vn; L.gps_ve = J.gps_ve; L.gps_vd = J.gps_vd;
L.gps_spd = J.gps_spd; L.gps_gcrs = J.gps_gcrs; L.gps_nsats = J.gps_nsats;
L.boundaries = boundaries; % tempos onde cada log termina (no time-base unificado)
L.log_starts = log_starts; % tempos onde cada log começa
L.log_names  = LOGS;
L.gap        = GAP;

out_file = fullfile(LOG_DIR, 'logs_concat.mat');
save(out_file, 'L', '-v7.3');
fprintf('Salvo: %s\n', out_file);

%% ====================== Helper p/ marcar boundaries =================
mark_logs = @() arrayfun(@(t) xline(t, 'w--', 'Alpha', 0.5, ...
    'HandleVisibility','off'), boundaries);
mark_labels = @() arrayfun(@(i) text(log_starts(i)+5, 0, ...
    sprintf('LOG %d', i+3), 'Color', 'y', 'FontWeight', 'bold', ...
    'VerticalAlignment','top','HandleVisibility','off'), 1:numel(LOGS));

%% =========================================================================
%  Fig 1 — RCOU (PWM dos motores)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 400]);
plot(J.time_RCOU, J.pwm1, 'LineWidth', 0.8); hold on
plot(J.time_RCOU, J.pwm2, 'LineWidth', 0.8);
plot(J.time_RCOU, J.pwm3, 'LineWidth', 0.8);
plot(J.time_RCOU, J.pwm4, 'LineWidth', 0.8);
mark_logs();
xlabel('Tempo unificado [s]'); ylabel('PWM [\mus]');
title('RCOU — PWM enviado aos motores (todos os logs concatenados)');
legend('M1 (FR)','M2 (RL)','M3 (FL)','M4 (RR)','Location','best');
grid on;

%% =========================================================================
%  Fig 2 — RCIN (comando do piloto via rádio)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 600]);
subplot(4,1,1);
plot(J.time_RCIN, J.rc_roll, 'b'); grid on
mark_logs(); ylabel('Roll [\mus]'); ylim([1000 2000]); yline(1500,'k--');
title('RCIN — comando do piloto');
subplot(4,1,2);
plot(J.time_RCIN, J.rc_pitch, 'r'); grid on
mark_logs(); ylabel('Pitch [\mus]'); ylim([1000 2000]); yline(1500,'k--');
subplot(4,1,3);
plot(J.time_RCIN, J.rc_thr, 'Color',[0.1 0.6 0.1]); grid on
mark_logs(); ylabel('Throttle [\mus]'); ylim([1000 2000]); yline(1500,'k--');
subplot(4,1,4);
plot(J.time_RCIN, J.rc_yaw, 'm'); grid on
mark_logs(); ylabel('Yaw [\mus]'); ylim([1000 2000]); yline(1500,'k--');
xlabel('Tempo unificado [s]');

%% =========================================================================
%  Fig 3 — pqr (gyro)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 600]);
subplot(3,1,1)
plot(J.time_IMU, J.gyrX, 'b'); grid on; mark_logs();
title('pqr — velocidades angulares [rad/s]'); legend('p');
subplot(3,1,2)
plot(J.time_IMU, J.gyrY, 'r'); grid on; mark_logs(); legend('q');
subplot(3,1,3)
plot(J.time_IMU, J.gyrZ, 'Color',[0.1 0.6 0.1]); grid on; mark_logs(); legend('r');
xlabel('Tempo unificado [s]');

%% =========================================================================
%  Fig 4 — Acc (acelerações)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 600]);
subplot(3,1,1)
plot(J.time_IMU, J.accX, 'b'); grid on; mark_logs();
title('Aceleração [m/s²]'); legend('AccX');
subplot(3,1,2)
plot(J.time_IMU, J.accY, 'r'); grid on; mark_logs(); legend('AccY');
subplot(3,1,3)
plot(J.time_IMU, J.accZ, 'Color',[0.1 0.6 0.1]); grid on; mark_logs(); legend('AccZ');
xlabel('Tempo unificado [s]');

%% =========================================================================
%  Fig 5 — ATT (atitude em graus)
%  =========================================================================
figure('Color','w','Position',[80 80 1400 600]);
subplot(3,1,1)
plot(J.time_ATT, J.roll, 'b'); grid on; mark_logs();
title('Atitude [deg]'); legend('Roll \phi');
subplot(3,1,2)
plot(J.time_ATT, J.pitch, 'r'); grid on; mark_logs(); legend('Pitch \theta');
subplot(3,1,3)
plot(J.time_ATT, J.yaw, 'Color',[0.1 0.6 0.1]); grid on; mark_logs(); legend('Yaw \psi');
xlabel('Tempo unificado [s]');

%% =========================================================================
%  Fig 6 — MAG (magnetômetro) — alimenta o yaw do EKF
%  =========================================================================
if ~isempty(J.time_MAG)
    figure('Color','w','Position',[80 80 1400 500]);
    subplot(2,1,1); hold on; grid on;
    plot(J.time_MAG, J.magX, 'r', 'DisplayName','MagX');
    plot(J.time_MAG, J.magY, 'g', 'DisplayName','MagY');
    plot(J.time_MAG, J.magZ, 'b', 'DisplayName','MagZ');
    mark_logs(); ylabel('Campo [mGauss]'); legend('Location','best');
    title('MAG — magnetômetro (referência de heading p/ o yaw)');
    subplot(2,1,2); hold on; grid on;
    plot(J.time_MAG, sqrt(J.magX.^2 + J.magY.^2 + J.magZ.^2), 'k');
    mark_logs(); ylabel('|mag| [mGauss]'); xlabel('Tempo unificado [s]');
    title('Norma do campo (deve ser ~constante se bem calibrado)');
end

%% =========================================================================
%  Fig 7 — GPS (velocidade NED) — alimenta o tilt do EKF
%  =========================================================================
if ~isempty(J.time_GPS)
    figure('Color','w','Position',[80 80 1400 500]);
    subplot(2,1,1); hold on; grid on;
    plot(J.time_GPS, J.gps_vn, 'r', 'DisplayName','V_N');
    plot(J.time_GPS, J.gps_ve, 'g', 'DisplayName','V_E');
    plot(J.time_GPS, J.gps_vd, 'b', 'DisplayName','V_D');
    mark_logs(); ylabel('Velocidade [m/s]'); legend('Location','best');
    title('GPS — velocidade NED (só trechos com fix 3D)');
    subplot(2,1,2); hold on; grid on;
    plot(J.time_GPS, J.gps_spd, 'k');
    mark_logs(); ylabel('Speed [m/s]'); xlabel('Tempo unificado [s]');
    title('GPS — velocidade horizontal');
end

%% =========================================================================
%  Resumo
%  =========================================================================
fprintf('\n================== RESUMO LOGS CONCATENADOS ====================\n');
fprintf('Total: %.1f s (%.1f min) | %d logs | gap entre logs = %d s\n', ...
    t_total, t_total/60, numel(LOGS), GAP);
fprintf('IMU:  %d amostras  | dt_med=%.4f s\n', numel(J.time_IMU), median(diff(J.time_IMU(J.time_IMU<boundaries(1)))));
fprintf('ATT:  %d amostras\n', numel(J.time_ATT));
fprintf('RCOU: %d amostras\n', numel(J.time_RCOU));
fprintf('RCIN: %d amostras\n', numel(J.time_RCIN));
fprintf('MAG:  %d amostras\n', numel(J.time_MAG));
fprintf('GPS:  %d amostras (velocidade NED, só fixes 3D)\n', numel(J.time_GPS));
fprintf('===============================================================\n');


%% ====================== Helper: leitura de GPS ======================
function G = read_gps_(filename)
%READ_GPS_  Lê velocidade do GPS de um log (formatos mp_export e legacy).
%  Retorna só fixes 3D (Status>=3). Velocidade em NED [m/s].
    G = struct('time',[],'vn',[],'ve',[],'vd',[],'spd',[],'gcrs',[],'nsats',[]);
    info = whos('-file', filename);
    nm = {info.name};
    if any(strcmp(nm,'GPS_0'))            % mp_export (matriz)
        S = load(filename,'GPS_0');  M = S.GPS_0;
        % cols ArduPilot: TimeUS=2, Status=4, NSats=7, Spd=12, GCrs=13, VZ=14
        t = double(M(:,2))/1e6; status = M(:,4); nsats = M(:,7);
        spd = M(:,12); gcrs = M(:,13); vz = M(:,14);
    elseif any(strcmp(nm,'GPS'))          % legacy (struct)
        S = load(filename,'GPS');  g = S.GPS;
        t = double(g.TimeUS)/1e6; status = double(g.Status);
        if isfield(g,'NSats'), nsats = double(g.NSats); else, nsats = nan(size(t)); end
        spd = double(g.Spd); gcrs = double(g.GCrs); vz = double(g.VZ);
    else
        return;                            % sem GPS neste log
    end
    ok = status >= 3;                      % só 3D fix (descarta logs sem satélite)
    t=t(ok); spd=spd(ok); gcrs=gcrs(ok); vz=vz(ok); nsats=nsats(ok);
    if isempty(t), return; end
    [t, ia] = unique(t, 'stable');         % dedup defensivo de timestamps
    spd=spd(ia); gcrs=gcrs(ia); vz=vz(ia); nsats=nsats(ia);
    G.time  = t;
    G.vn    = spd .* cosd(gcrs);           % VN = Spd*cos(curso)
    G.ve    = spd .* sind(gcrs);           % VE = Spd*sin(curso)
    G.vd    = vz;                          % VZ do ArduPilot (NED, +p/baixo — conferir sinal)
    G.spd   = spd;  G.gcrs = gcrs;  G.nsats = nsats;
end


%% ====================== Helper: leitura de MAG ======================
function Mg = read_mag_(filename)
%READ_MAG_  Lê o magnetômetro ATIVO. Prefere MAG_1 (compass 2 externo, usado
%  pelo EKF neste drone); fallback MAG_0. Formato mp_export (matriz).
    Mg = struct('time',[],'x',[],'y',[],'z',[]);
    info = whos('-file', filename);  nm = {info.name};
    for c = {'MAG_1','MAG_0'}                 % preferência: compass 2 (externo/ativo)
        if any(strcmp(nm, c{1}))
            S = load(filename, c{1});  M = S.(c{1});
            if isstruct(M), continue; end     % legacy struct: deixa o fallback cuidar
            % cols mp_export: TimeUS=2, MagX=4, MagY=5, MagZ=6
            t = double(M(:,2))/1e6;
            [t, ia] = unique(t, 'stable');
            Mg.time = t;
            Mg.x = double(M(ia,4));  Mg.y = double(M(ia,5));  Mg.z = double(M(ia,6));
            return;
        end
    end
end
