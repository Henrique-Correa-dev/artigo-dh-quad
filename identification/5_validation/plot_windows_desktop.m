function plot_windows_desktop(varargin)
%PLOT_WINDOWS_DESKTOP  Uma figura por trecho de validação, salva na Área de Trabalho.
%
%   plot_windows_desktop()
%   plot_windows_desktop('wins', [610 630; ...], 'dir', '/caminho')
%
% Cada figura tem seis painéis empilhados (1 coluna, 6 linhas): p, q, r, a_x,
% a_y, a_z, com o medido em preto e as duas formas de amortecimento sobrepostas
% (taxa e momento). Simulação em modo 'full' (malha aberta, 9 estados).
    opt = struct('wins', [605 625; 610 630; 550 600; 400 420; 420 440; 440 460; ...
                          460 480; 525 545; 545 565; 565 585; 585 605; 630 650], ...
                 'dir', fullfile(getenv('HOME'), 'Desktop', 'DH_validacao_trechos'));
    for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
    if ~exist(opt.dir,'dir'), mkdir(opt.dir); end
    paths = setup_paths();  proj = parameters();

    % ANTES: modelo oficial da dissertação — amortecimento como coeficiente de
    %        taxa (c_p, c_q, c_r), multijanelamento de 1, 2 e 3 s.
    % AGORA: amortecimento dentro do vetor de momentos (L_p, M_q, N_r em N·m·s),
    %        chute inicial do balanço físico e multijanelamento até 20 s.
    C = struct('nome', {'antes: taxa, sem aero, janelas 1/2/3 s','agora: momento + asa medida em voo, janelas até 20 s'}, ...
               'file', {fullfile(paths.outputs,'P_identified_rateform_backup.mat'), ...
                        fullfile(paths.outputs,'runs','oficial_2026_final','P_identified.mat')}, ...
               'form', {'rate','moment'}, 'cor', {[0 0.45 0.7],[0.85 0.37 0.01]});
    for k = 1:numel(C), C(k).P = load(C(k).file).P_final(:); end

    L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
    t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
    t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
    tg = (t_lo:0.1:t_hi)';
    set_aero_vel(L, tg);      % velocidade exógena dos termos ∝ V
    setappdata(0,'faero_on','measured');
    gapA = L.boundaries(end);  gapB = L.log_starts(end);
    constants = struct('m', proj.m, 'g', proj.g);
    R2 = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
    % Coeficiente de desigualdade de Theil: 0 = ajuste perfeito, 1 = nenhum.
    % Diferente do R², não é enganado por sinal de variância pequena.
    TIC = @(y,yh) sqrt(mean((y-yh).^2)) / (sqrt(mean(y.^2)) + sqrt(mean(yh.^2)));

    canais = {'p','q','r','accX','accY','accZ'};
    rot    = {'p [rad/s]','q [rad/s]','r [rad/s]','a_x [m/s²]','a_y [m/s²]','a_z [m/s²]'};

    for w = 1:size(opt.wins,1)
        tw = opt.wins(w,:);
        if tw(1) < gapB && tw(2) > gapA
            fprintf('  %4.0f–%-4.0f  pulada (cruza o intervalo entre logs)\n', tw); continue;
        end
        idx = tg>=tw(1) & tg<=tw(2);  time = tg(idx);
        ip = @(tt,xx) interp1(tt, xx, time, 'linear');
        pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ...
               ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
        pqr = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
        acc = [ip(L.time_IMU,L.accX_raw), ip(L.time_IMU,L.accY_raw), ip(L.time_IMU,L.accZ_raw)];
        att = [ip(L.time_ATT,L.roll_deg), ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
        med = [pqr, acc];

        S = cell(numel(C),1);  R = nan(numel(C),6);  Tc = nan(numel(C),6);
        for k = 1:numel(C)
            setappdata(0,'damp_form_override', C(k).form);  clear aero_gains
            S{k} = sim_window('full', C(k).P, time, pwm, pqr, att, constants);
            for j = 1:6
                R(k,j)  = R2( med(:,j), S{k}.(canais{j}));
                Tc(k,j) = TIC(med(:,j), S{k}.(canais{j}));
            end
        end
        if isappdata(0,'damp_form_override'), rmappdata(0,'damp_form_override'); end

        f = figure('Position',[30 30 1150 1250],'Color','w','Visible','off');
        try, f.Theme='light'; catch, end
        tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
        for j = 1:6
            nexttile; hold on; grid on;
            plot(time, med(:,j), 'k-', 'LineWidth', 1.7);
            for k = 1:numel(C)
                plot(time, S{k}.(canais{j}), '-', 'Color', C(k).cor, 'LineWidth', 1.2);
            end
            ylabel(rot{j}); xlim(tw);
            text(0.006, 0.94, sprintf(['R²    antes %6.3f   agora %6.3f   \\Delta %+.3f\n' ...
                                       'TIC   antes %6.3f   agora %6.3f   \\Delta %+.3f'], ...
                R(1,j), R(2,j), R(2,j)-R(1,j), Tc(1,j), Tc(2,j), Tc(2,j)-Tc(1,j)), ...
                'Units','normalized','VerticalAlignment','top','FontSize',8.5, ...
                'FontName','Menlo', 'BackgroundColor',[1 1 1 0.78]);
            if j == 1
                legend([{'medido'}, {C.nome}], 'Location','northoutside', ...
                       'Orientation','horizontal', 'FontSize',9);
            end
            if j == 6, xlabel('tempo [s]'); end
        end
        title(tl, sprintf(['Validação em malha aberta, trecho %.0f a %.0f s   (voo 25/05/2026)\n' ...
              'antes: c_p,c_q,c_r fora dos momentos, sem aerodinâmica   |   agora: L_p,M_q,N_r dentro de M + asa com C_{lp},C_{mq},C_{nr} medidos em voo'], tw), ...
              'FontWeight','bold');
        fn = fullfile(opt.dir, sprintf('antes_x_agora_%04.0f-%04.0f.png', tw));
        exportgraphics(f, fn, 'BackgroundColor','white','Resolution',130);
        close(f);
        fprintf('  %4.0f–%-4.0f | R²  %.3f %.3f %.3f → %.3f %.3f %.3f | TIC %.3f %.3f %.3f → %.3f %.3f %.3f\n', ...
            tw, R(1,1:3), R(2,1:3), Tc(1,1:3), Tc(2,1:3));
    end
    if isappdata(0,'faero_on'), rmappdata(0,'faero_on'); end
    fprintf('\n  %d figuras em %s\n', size(opt.wins,1), opt.dir);
end
