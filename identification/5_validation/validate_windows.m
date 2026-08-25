function T = validate_windows(varargin)
%VALIDATE_WINDOWS  R² de p, q, r em todos os trechos de validação combinados.
%
%   T = validate_windows()                        compara momento x taxa
%   T = validate_windows('wins', [610 630; ...])  outras janelas
%
% Simula em modo 'full' (open loop, 9 estados, mesma cadeia de motor) e devolve
% uma tabela com R² por canal e por janela, para cada vetor de parâmetros.
% As duas formas de amortecimento são alternadas pelo appdata damp_form_override,
% de modo que cada vetor é simulado com a forma em que foi identificado.
    opt = struct('wins', [605 625; 610 630; 550 600; 400 420; 420 440; 440 460; ...
                          460 480; 525 545; 545 565; 565 585; 585 605; 630 650]);
    for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
    paths = setup_paths();  proj = parameters();

    % Comparação: sem asa (L_p absorve tudo) × oficial (asa separada), cada um
    % com o seu Θ₀ e o seu identificado. P0 é lido do próprio arquivo de cada rodada.
    f_noaero = fullfile(paths.outputs,'runs','moment_win6_2026','P_identified.mat');
    f_ofic   = fullfile(paths.outputs,'runs','oficial_2026_final','P_identified.mat');
    C = struct('nome', {'SEM ASA: Θ₀','SEM ASA: identificado','OFICIAL: Θ₀','OFICIAL: identificado'}, ...
               'file', {f_noaero, f_noaero, f_ofic, f_ofic}, 'use_P0', {true,false,true,false}, ...
               'form', {'moment','moment','moment','moment'}, 'hu', {false,false,false,false}, ...
               'faero', {'','','measured','measured'});
    if isfield(opt,'only') && ~isempty(opt.only), C = C(opt.only); end

    L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
    t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
    t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
    tg = (t_lo:0.1:t_hi)';
    set_aero_vel(L, tg);      % velocidade exógena dos termos ∝ V
    gapA = L.boundaries(end);  gapB = L.log_starts(end);
    constants = struct('m', proj.m, 'g', proj.g);
    R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

    for k = 1:numel(C)
        Sk = load(C(k).file);
        if isfield(C,'use_P0') && C(k).use_P0, C(k).P = Sk.P0(:); else, C(k).P = Sk.P_final(:); end
    end
    W = opt.wins;  nW = size(W,1);
    res = nan(nW, numel(C), 6);   % p, q, r, a_x, a_y, a_z
    nC = numel(C);
    fprintf('\n  %-11s', 'janela');
    for k = 1:nC, fprintf(' | %-23s', C(k).nome); end, fprintf('\n');
    fprintf('  %-11s', '');
    for k = 1:nC, fprintf(' | %7s %7s %7s', 'R² p','R² q','R² r'); end, fprintf('\n');
    for w = 1:nW
        tw = W(w,:);
        if tw(1) < gapB && tw(2) > gapA
            fprintf('  %4.0f–%-6.0f| (cruza o intervalo entre logs — pulada)\n', tw); continue;
        end
        idx = tg>=tw(1) & tg<=tw(2);  time = tg(idx);
        ip = @(tt,xx) interp1(tt, xx, time, 'linear');
        pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ...
               ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
        pqr = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
        acc = [ip(L.time_IMU,L.accX_raw), ip(L.time_IMU,L.accY_raw), ip(L.time_IMU,L.accZ_raw)];
        att = [ip(L.time_ATT,L.roll_deg), ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
        med = [pqr, acc];
        for k = 1:numel(C)
            setappdata(0,'damp_form_override', C(k).form);  clear aero_gains
            if isfield(C,'hu') && C(k).hu, setappdata(0,'hu_form',true);
            elseif isappdata(0,'hu_form'), rmappdata(0,'hu_form'); end
            if isfield(C,'faero') && ~isempty(C(k).faero), setappdata(0,'faero_on',C(k).faero);
            elseif isappdata(0,'faero_on'), rmappdata(0,'faero_on'); end
            try
                r = sim_window('full', C(k).P, time, pwm, pqr, att, constants);
                sim = [r.p, r.q, r.r, r.accX, r.accY, r.accZ];
                for j = 1:6, res(w,k,j) = R2(med(:,j), sim(:,j)); end
            catch ME
                fprintf('   (janela %d, %s: %s)\n', w, C(k).nome, ME.message);
            end
        end
        fprintf('  %4.0f–%-6.0f', tw);
        for k = 1:nC, fprintf(' | %7.3f %7.3f %7.3f', res(w,k,1), res(w,k,2), res(w,k,3)); end, fprintf('\n');
    end
    rmappdata(0,'damp_form_override');
    if isappdata(0,'hu_form'), rmappdata(0,'hu_form'); end
    if isappdata(0,'faero_on'), rmappdata(0,'faero_on'); end
    ok = all(isfinite(res(:,1,1)),2);
    fprintf('  %-11s', 'MÉDIA');
    for k = 1:nC, fprintf(' | %7.3f %7.3f %7.3f', mean(res(ok,k,1)), mean(res(ok,k,2)), mean(res(ok,k,3))); end
    fprintf('\n');

    % ---- segundo bloco: acelerômetro ----
    fprintf('\n  %-11s', 'janela');
    for k = 1:nC, fprintf(' | %-23s', C(k).nome); end, fprintf('\n');
    fprintf('  %-11s', '');
    for k = 1:nC, fprintf(' | %7s %7s %7s', 'R² a_x','R² a_y','R² a_z'); end, fprintf('\n');
    for w = 1:nW
        if ~isfinite(res(w,1,1)), continue; end
        fprintf('  %4.0f–%-6.0f', W(w,:));
        for k = 1:nC, fprintf(' | %7.3f %7.3f %7.3f', res(w,k,4), res(w,k,5), res(w,k,6)); end, fprintf('\n');
    end
    fprintf('  %-11s', 'MÉDIA');
    for k = 1:nC, fprintf(' | %7.3f %7.3f %7.3f', mean(res(ok,k,4)), mean(res(ok,k,5)), mean(res(ok,k,6))); end
    fprintf('\n');
    T = res;
    save(fullfile(paths.outputs,'validate_windows.mat'), 'res','W');
end
