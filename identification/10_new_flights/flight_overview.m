function flight_overview(matfile, tag, varargin)
%FLIGHT_OVERVIEW  Visão geral de um log ArduPilot exportado pelo Mission Planner (.mat)
%   flight_overview(matfile, tag)  → outputs/images/new_flights/<tag>_overview.png
%   Painéis: atitude, taxas, saídas de motor/servo (RCOU), velocidades (GPS,
%   EKF, pitot se houver), altitude barométrica e modo de voo. Serve para
%   escolher janelas de identificação/validação ANTES de qualquer processamento.
    opt = struct('tmax', Inf);
    for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
    paths = setup_paths();
    outdir = fullfile(paths.images, 'new_flights');  if ~exist(outdir,'dir'), mkdir(outdir); end
    % alguns exports têm variáveis com nome inválido ('s-s_label'): carrega só o que interessa
    need = {'ATT','IMU','IMU_0','RCOU','BARO','BARO_0','GPS','GPS_0','XKF1','XKF1_0','ARSP','MODE','AOA'};
    need = [need, strcat(need,'_label')];
    w = whos('-file', matfile);  have = intersect({w.name}, need);
    S = load(matfile, have{:});
    G = @(msg) getmsg(S, msg);            % struct com .t e campos por nome (via *_label)

    ATT = G('ATT');  IMU = G('IMU');  RCOU = G('RCOU');  BARO = G('BARO');
    GPS = G('GPS');  XKF1 = G('XKF1'); ARSP = G('ARSP'); MODE = G('MODE'); AOA = G('AOA');
    t0 = ATT.t(1);  tt = @(M) M.t - t0;
    tmax = min(opt.tmax, ATT.t(end)-t0);

    f = figure('Position',[30 30 1300 1150],'Color','w'); try, f.Theme='light'; catch, end
    tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
    % 1 atitude
    nexttile; hold on; grid on;
    plot(tt(ATT), ATT.Roll, 'LineWidth',1); plot(tt(ATT), ATT.Pitch, 'LineWidth',1); plot(tt(ATT), ATT.Yaw/10, 'LineWidth',0.8);
    ylabel('[°]'); legend({'roll','pitch','yaw/10'},'Location','eastoutside'); xlim([0 tmax]);
    ttl(sprintf('%s — %s  (%.0f s)', tag, strrep(basename(matfile),'_','\_'), tmax));
    % 2 taxas
    nexttile; hold on; grid on;
    plot(tt(IMU), IMU.GyrX, 'LineWidth',0.8); plot(tt(IMU), IMU.GyrY, 'LineWidth',0.8); plot(tt(IMU), IMU.GyrZ, 'LineWidth',0.8);
    ylabel('[rad/s]'); legend({'p','q','r'},'Location','eastoutside'); xlim([0 tmax]); ttl('taxas (giroscópio)');
    % 3 RCOU (todos os canais não constantes)
    nexttile; hold on; grid on;
    ch = {}; for c = 1:14
        fn = sprintf('C%d',c); if isfield(RCOU,fn) && range(RCOU.(fn)) > 5, plot(tt(RCOU), RCOU.(fn), 'LineWidth',0.8); ch{end+1} = fn; end %#ok<AGROW>
    end
    ylabel('PWM [µs]'); legend(ch,'Location','eastoutside'); xlim([0 tmax]); ttl('saídas RCOU (motores/servos com atividade)');
    % 4 velocidades
    nexttile; hold on; grid on; lg = {};
    if ~isempty(GPS) && isfield(GPS,'Spd'), plot(tt(GPS), GPS.Spd, 'k-', 'LineWidth',1); lg{end+1} = sprintf('GPS Spd (nsats %d–%d)', min(GPS.NSats), max(GPS.NSats)); end
    if ~isempty(XKF1) && isfield(XKF1,'VN'), plot(tt(XKF1), hypot(XKF1.VN,XKF1.VE), '-', 'Color',[0 0.45 0.7], 'LineWidth',1); lg{end+1} = 'EKF |V_h|'; plot(tt(XKF1), XKF1.VD, '-', 'Color',[0.85 0.37 0.01], 'LineWidth',0.8); lg{end+1} = 'EKF V_D'; end
    if ~isempty(ARSP) && isfield(ARSP,'Airspeed'), plot(tt(ARSP), ARSP.Airspeed, '-', 'Color',[0.3 0.6 0.3], 'LineWidth',1.2); lg{end+1} = 'pitot (ARSP)'; end
    ylabel('[m/s]'); if ~isempty(lg), legend(lg,'Location','eastoutside'); end; xlim([0 tmax]); ttl('velocidades');
    % 5 altitude / AOA
    nexttile; hold on; grid on; lg = {};
    if ~isempty(BARO) && isfield(BARO,'Alt'), plot(tt(BARO), BARO.Alt, 'k-', 'LineWidth',1); lg{end+1} = 'baro Alt [m]'; end
    if ~isempty(AOA) && isfield(AOA,'AOA'), yyaxis right; plot(tt(AOA), AOA.AOA, '-', 'Color',[0.85 0.37 0.01]); ylabel('AOA [°]'); yyaxis left; lg{end+1} = 'AOA (EKF) →'; end
    if ~isempty(lg), legend(lg,'Location','eastoutside'); end; xlim([0 tmax]); ttl('altitude barométrica (e AOA se houver)');
    % 6 modo
    nexttile; hold on; grid on;
    if ~isempty(MODE) && isfield(MODE,'Mode')
        stairs([tt(MODE); tmax], [MODE.Mode; MODE.Mode(end)], 'k-', 'LineWidth',1.5);
        tm = tt(MODE);
        for k = 1:numel(MODE.t), text(tm(k), MODE.Mode(k)+0.3, mode_name(MODE.Mode(k), S), 'FontSize',8); end
    end
    ylabel('modo'); xlim([0 tmax]); xlabel('t [s] (desde o 1º ATT)'); ttl('modo de voo');
    fn = fullfile(outdir, [tag '_overview.png']);
    exportgraphics(f, fn, 'BackgroundColor','white','Resolution',110);
    fprintf('  %s\n', fn);
end

function ttl(s), text(0.005, 1.04, s, 'Units','normalized', 'FontWeight','bold', 'FontSize',9); end
function b = basename(p), [~,b,e] = fileparts(p); b = [b e]; end

function M = getmsg(S, msg)
%GETMSG  Extrai mensagem `msg` (ou `msg_0`) como struct por campo, tempo em s.
    M = [];
    cand = {msg, [msg '_0']};
    for c = 1:2
        if isfield(S, cand{c}) && isnumeric(S.(cand{c}))
            X = S.(cand{c});  lab = [];
            if isfield(S, [msg '_label']), lab = S.([msg '_label']); end
            if isempty(lab), return; end
            M = struct('t', X(:,2)/1e6);
            for j = 1:min(numel(lab), size(X,2))
                nm = matlab.lang.makeValidName(char(lab{j}));
                if any(strcmp(nm, {'LineNo','TimeUS'})), continue; end
                M.(nm) = X(:,j);
            end
            return;
        end
    end
end

function s = mode_name(m, S)
    % Copter/Plane mode numbers diferem; usa a string do MODE se houver ("Mode" numérico + "ModeNum")
    s = sprintf('%d', m);
    if isfield(S,'MODE_label') && any(strcmp(S.MODE_label,'ModeNum')), return; end
end
