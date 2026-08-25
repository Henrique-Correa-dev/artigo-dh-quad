function plot_oficial_desktop(varargin)
%PLOT_OFICIAL_DESKTOP  Modelo oficial sozinho, uma figura por trecho, na Mesa.
%
%   plot_oficial_desktop()
%
% Seis painéis empilhados (p, q, r, a_x, a_y, a_z): medido em preto, modelo
% identificado (oficial_fw_2026) em laranja, chute inicial Θ₀ em cinza
% tracejado. R² e TIC de cada um no canto. Modo full, 9 estados, malha aberta.
    opt = struct('wins', [605 625; 610 630; 550 600; 400 420; 420 440; 440 460; ...
                          460 480; 525 545; 545 565; 565 585; 585 605; 630 650], ...
                 'dir', fullfile(getenv('HOME'), 'Desktop', 'DH_modelo_oficial'));
    for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
    if ~exist(opt.dir,'dir'), mkdir(opt.dir); end
    paths = setup_paths();  proj = parameters();

    S = load(fullfile(paths.outputs,'runs','oficial_2026_final','P_identified.mat'));
    P  = S.P_final(:);  P0 = S.P0(:);

    L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
    t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
    t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
    tg = (t_lo:0.1:t_hi)';
    set_aero_vel(L, tg);
    setappdata(0,'damp_form_override','moment');  clear aero_gains
    setappdata(0,'faero_on','measured');     % F_aero da asa com C_L0, C_Lα, C_Yβ, C_D0 medidos
    gapA = L.boundaries(end);  gapB = L.log_starts(end);
    constants = struct('m', proj.m, 'g', proj.g);
    R2  = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
    TIC = @(y,yh) sqrt(mean((y-yh).^2)) / (sqrt(mean(y.^2)) + sqrt(mean(yh.^2)));

    canais = {'p','q','r','accX','accY','accZ'};
    rot    = {'p [rad/s]','q [rad/s]','r [rad/s]','a_x [m/s²]','a_y [m/s²]','a_z [m/s²]'};
    fprintf('\n  Modelo: L_p %.4f | M_q %.4f | N_r %.4f N·m·s | C_lp %.4f C_mq %.3f C_nr %.4f\n', ...
        P(13), P(14), P(15), P(17), P(19), P(21));
    for w = 1:size(opt.wins,1)
        tw = opt.wins(w,:);
        if tw(1) < gapB && tw(2) > gapA, continue; end
        idx = tg>=tw(1) & tg<=tw(2);  time = tg(idx);
        ip = @(tt,xx) interp1(tt, xx, time, 'linear');
        pwm = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ...
               ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
        pqr = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
        acc = [ip(L.time_IMU,L.accX_raw), ip(L.time_IMU,L.accY_raw), ip(L.time_IMU,L.accZ_raw)];
        att = [ip(L.time_ATT,L.roll_deg), ip(L.time_ATT,L.pitch_deg), ip(L.time_ATT,L.yaw_deg)];
        med = [pqr, acc];

        r1 = sim_window('full', P,  time, pwm, pqr, att, constants);
        r0 = sim_window('full', P0, time, pwm, pqr, att, constants);
        f = figure('Position',[30 30 1150 1250],'Color','w','Visible','off');
        try, f.Theme='light'; catch, end
        tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
        Rr = zeros(6,2);
        for j = 1:6
            nexttile; hold on; grid on;
            y1 = r1.(canais{j});  y0 = r0.(canais{j});
            plot(time, y0, '--', 'Color',[0.6 0.6 0.6], 'LineWidth', 1.0);
            plot(time, med(:,j), 'k-', 'LineWidth', 1.7);
            plot(time, y1, '-', 'Color',[0.85 0.37 0.01], 'LineWidth', 1.3);
            Rr(j,:) = [R2(med(:,j),y1), TIC(med(:,j),y1)];
            ylabel(rot{j}); xlim(tw);
            text(0.006, 0.94, sprintf('R² %.3f   TIC %.3f   (\\Theta_0: R² %.2f)', Rr(j,1), Rr(j,2), R2(med(:,j),y0)), ...
                'Units','normalized','VerticalAlignment','top','FontSize',9, ...
                'FontName','Menlo','BackgroundColor',[1 1 1 0.78]);
            if j == 1
                legend({'\Theta_0 (chute inicial)','medido','identificado'}, ...
                    'Location','northoutside','Orientation','horizontal','FontSize',9);
            end
            if j == 6, xlabel('tempo [s]'); end
        end
        title(tl, sprintf(['Modelo oficial, validação em malha aberta, trecho %.0f a %.0f s (voo 25/05/2026)\n' ...
            'M = M_{rot} − [L_p p; M_q q; N_r r] + M_{asa}·V_a      F = F_{rot} + F_{asa}·V_a      ' ...
            '(todos os coeficientes da asa medidos em voo de asa fixa)'], tw), ...
            'FontWeight','bold');
        fn = fullfile(opt.dir, sprintf('oficial_%04.0f-%04.0f.png', tw));
        exportgraphics(f, fn, 'BackgroundColor','white','Resolution',130);
        close(f);
        fprintf('  %4.0f–%-4.0f  R² %.3f %.3f %.3f | %.3f %.3f %.3f   TIC %.3f %.3f %.3f | %.3f %.3f %.3f\n', ...
            tw, Rr(:,1), Rr(:,2));
    end
    rmappdata(0,'damp_form_override');  rmappdata(0,'faero_on');
    fprintf('\n  figuras em %s\n', opt.dir);
end
