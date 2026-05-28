function L = load_log_data(filename)
%LOAD_LOG_DATA  Carrega log do drone em formato unificado, suportando 2 layouts:
%
%   FORMATO LEGADO (log_data.mat original):
%       struct ATT  com campos TimeUS, Roll, Pitch, Yaw
%       struct IMU  com campos TimeUS, GyrX/Y/Z, AccX/Y/Z, I
%       struct RCOU com campos TimeUS, C1, C2, C3, C4
%       struct GPS  com campos TimeUS (opcional)
%
%   FORMATO NOVO (Mission Planner export, ex. "X.bin-NNNNNN.mat"):
%       matriz ATT       Nx9  (col 2=TimeUS, 4=Roll, 6=Pitch, 8=Yaw)
%       matriz IMU_0     Nx16 (col 2=TimeUS, 4=GyrX, 5=GyrY, 6=GyrZ,
%                              7=AccX, 8=AccY, 9=AccZ) — IMU primário
%       matriz IMU_1     Nx16 (mesma, secundário) — opcional
%       matriz RCOU      Nx16 (col 2=TimeUS, 3=C1, 4=C2, 5=C3, 6=C4)
%       matriz GPS_0     Nx16 (col 2=TimeUS) — opcional
%
% Detecta o formato automaticamente pela presença de IMU_0 (matriz) vs IMU (struct).
%
% Saída uniforme (struct L) com tempo já em SEGUNDOS:
%   L.time_IMU   (Nx1)
%   L.gyrX_raw, L.gyrY_raw, L.gyrZ_raw  (Nx1)
%   L.accX_raw, L.accY_raw, L.accZ_raw  (Nx1)
%   L.time_ATT   (Mx1)
%   L.roll_deg, L.pitch_deg, L.yaw_deg  (Mx1)   — em graus
%   L.time_RCOU  (Kx1)
%   L.pwm1_raw, L.pwm2_raw, L.pwm3_raw, L.pwm4_raw  (Kx1)
%   L.time_RCIN  (Jx1, opcional [] se ausente)
%   L.rcin_roll, L.rcin_pitch, L.rcin_throttle, L.rcin_yaw  (Jx1, em PWM us)
%   L.time_GPS   (opcional, [] se não tiver)
%   L.format     string: 'legacy' ou 'mp_export'
%   L.source     filename original

    % Carregar APENAS variáveis válidas (evitar erro de field names como 's---_label')
    info = whos('-file', filename);
    valid_names = {};
    for k = 1:numel(info)
        n = info(k).name;
        if isvarname(n)   % filtra nomes inválidos automaticamente
            valid_names{end+1} = n; %#ok<AGROW>
        end
    end
    S = load(filename, valid_names{:});

    % Detectar formato
    has_concat     = isfield(S, 'L') && isstruct(S.L) && isfield(S.L, 'format') ...
                     && strcmp(S.L.format, 'concat');
    has_struct_IMU = isfield(S, 'IMU') && isstruct(S.IMU);
    has_matrix_IMU = isfield(S, 'IMU_0') && ismatrix(S.IMU_0) && ~isstruct(S.IMU_0);

    if has_concat
        % Arquivo já é pré-processado (gerado por test_flight.m)
        L = S.L;
        L.source = filename;
        % Dedup defensivo (concat antigo pode ter duplicatas dos logs originais)
        L = dedup_timestamps_(L);
        fprintf('load_log_data: %s (concat de %d logs) | %.1fs totais\n', ...
            filename, numel(L.log_names), L.time_IMU(end)-L.time_IMU(1));
        return;
    end

    L = struct();
    L.source = filename;

    if has_struct_IMU
        L.format = 'legacy';
        L = load_legacy_(L, S);
    elseif has_matrix_IMU
        L.format = 'mp_export';
        L = load_mp_export_(L, S);
    else
        error('load_log_data:unknownFormat', ...
            'Não reconheci o formato do arquivo: %s', filename);
    end

    % Deduplica timestamps (logs ArduPilot às vezes têm 2 amostras no mesmo
    % TimeUS — quebra interp1). Mantém primeira ocorrência.
    L = dedup_timestamps_(L);

    fprintf('load_log_data: %s (%s) | IMU %.1fs | ATT %.1fs | RCOU %.1fs\n', ...
        filename, L.format, L.time_IMU(end)-L.time_IMU(1), ...
        L.time_ATT(end)-L.time_ATT(1), L.time_RCOU(end)-L.time_RCOU(1));
end


function L = dedup_timestamps_(L)
% Remove amostras com timestamp duplicado em cada série temporal.
% Mantém a primeira ocorrência. Garante monotonicidade estrita.
    function [t_u, mask] = unique_keep_first(t)
        [t_u, ia, ~] = unique(t, 'stable');
        mask = false(size(t)); mask(ia) = true;
    end

    n0 = numel(L.time_IMU);
    [L.time_IMU, m] = unique_keep_first(L.time_IMU);
    L.gyrX_raw = L.gyrX_raw(m); L.gyrY_raw = L.gyrY_raw(m); L.gyrZ_raw = L.gyrZ_raw(m);
    L.accX_raw = L.accX_raw(m); L.accY_raw = L.accY_raw(m); L.accZ_raw = L.accZ_raw(m);
    if numel(L.time_IMU) < n0
        fprintf('  dedup IMU:  removidas %d duplicatas\n', n0 - numel(L.time_IMU));
    end

    n0 = numel(L.time_ATT);
    [L.time_ATT, m] = unique_keep_first(L.time_ATT);
    L.roll_deg = L.roll_deg(m); L.pitch_deg = L.pitch_deg(m); L.yaw_deg = L.yaw_deg(m);
    if numel(L.time_ATT) < n0
        fprintf('  dedup ATT:  removidas %d duplicatas\n', n0 - numel(L.time_ATT));
    end

    n0 = numel(L.time_RCOU);
    [L.time_RCOU, m] = unique_keep_first(L.time_RCOU);
    L.pwm1_raw = L.pwm1_raw(m); L.pwm2_raw = L.pwm2_raw(m);
    L.pwm3_raw = L.pwm3_raw(m); L.pwm4_raw = L.pwm4_raw(m);
    if numel(L.time_RCOU) < n0
        fprintf('  dedup RCOU: removidas %d duplicatas\n', n0 - numel(L.time_RCOU));
    end

    if isfield(L, 'time_RCIN') && ~isempty(L.time_RCIN)
        n0 = numel(L.time_RCIN);
        [L.time_RCIN, m] = unique_keep_first(L.time_RCIN);
        L.rcin_roll = L.rcin_roll(m); L.rcin_pitch = L.rcin_pitch(m);
        L.rcin_throttle = L.rcin_throttle(m); L.rcin_yaw = L.rcin_yaw(m);
        if numel(L.time_RCIN) < n0
            fprintf('  dedup RCIN: removidas %d duplicatas\n', n0 - numel(L.time_RCIN));
        end
    end
end


function L = load_legacy_(L, S)
% Formato struct (log_data.mat original)
    ATT = S.ATT;  IMU = S.IMU;  RCOU = S.RCOU;

    idx = IMU.I == 0;                 % IMU primário
    L.time_IMU = double(IMU.TimeUS(idx)) / 1e6;
    L.gyrX_raw = IMU.GyrX(idx);  L.gyrY_raw = IMU.GyrY(idx);  L.gyrZ_raw = IMU.GyrZ(idx);
    L.accX_raw = IMU.AccX(idx);  L.accY_raw = IMU.AccY(idx);  L.accZ_raw = IMU.AccZ(idx);

    L.time_ATT  = double(ATT.TimeUS) / 1e6;
    L.roll_deg  = ATT.Roll;
    L.pitch_deg = ATT.Pitch;
    L.yaw_deg   = ATT.Yaw;

    L.time_RCOU = double(RCOU.TimeUS) / 1e6;
    L.pwm1_raw = double(RCOU.C1);  L.pwm2_raw = double(RCOU.C2);
    L.pwm3_raw = double(RCOU.C3);  L.pwm4_raw = double(RCOU.C4);

    if isfield(S, 'RCIN') && isstruct(S.RCIN)
        L.time_RCIN     = double(S.RCIN.TimeUS) / 1e6;
        L.rcin_roll     = double(S.RCIN.C1);
        L.rcin_pitch    = double(S.RCIN.C2);
        L.rcin_throttle = double(S.RCIN.C3);
        L.rcin_yaw      = double(S.RCIN.C4);
    else
        L.time_RCIN = []; L.rcin_roll = []; L.rcin_pitch = [];
        L.rcin_throttle = []; L.rcin_yaw = [];
    end

    if isfield(S, 'GPS') && isstruct(S.GPS)
        L.time_GPS = double(S.GPS.TimeUS) / 1e6;
    else
        L.time_GPS = [];
    end
end


function L = load_mp_export_(L, S)
% Formato matriz (Mission Planner export)
% Colunas (conforme ArduPilot DataFlash):
%   ATT:  [LineNo, TimeUS, DesRoll, Roll, DesPitch, Pitch, DesYaw, Yaw, AEKF]
%   IMU:  [LineNo, TimeUS, I, GyrX, GyrY, GyrZ, AccX, AccY, AccZ, ...]
%   RCOU: [LineNo, TimeUS, C1, C2, C3, C4, C5..C14]
%   GPS:  [LineNo, TimeUS, I, Status, ...]
    IMU  = S.IMU_0;
    ATT  = S.ATT;
    RCOU = S.RCOU;

    L.time_IMU = double(IMU(:,2)) / 1e6;
    L.gyrX_raw = IMU(:,4);  L.gyrY_raw = IMU(:,5);  L.gyrZ_raw = IMU(:,6);
    L.accX_raw = IMU(:,7);  L.accY_raw = IMU(:,8);  L.accZ_raw = IMU(:,9);

    L.time_ATT  = double(ATT(:,2)) / 1e6;
    L.roll_deg  = ATT(:,4);
    L.pitch_deg = ATT(:,6);
    L.yaw_deg   = ATT(:,8);

    L.time_RCOU = double(RCOU(:,2)) / 1e6;
    L.pwm1_raw = double(RCOU(:,3));  L.pwm2_raw = double(RCOU(:,4));
    L.pwm3_raw = double(RCOU(:,5));  L.pwm4_raw = double(RCOU(:,6));

    if isfield(S, 'RCIN')
        L.time_RCIN     = double(S.RCIN(:,2)) / 1e6;
        L.rcin_roll     = double(S.RCIN(:,3));
        L.rcin_pitch    = double(S.RCIN(:,4));
        L.rcin_throttle = double(S.RCIN(:,5));
        L.rcin_yaw      = double(S.RCIN(:,6));
    else
        L.time_RCIN = []; L.rcin_roll = []; L.rcin_pitch = [];
        L.rcin_throttle = []; L.rcin_yaw = [];
    end

    if isfield(S, 'GPS_0')
        L.time_GPS = double(S.GPS_0(:,2)) / 1e6;
    else
        L.time_GPS = [];
    end
end
