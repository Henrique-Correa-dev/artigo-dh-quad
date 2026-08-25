function fw_zoom(which_log, wins, tag)
%FW_ZOOM  Amplia janelas do voo de asa fixa para ver a forma do comando.
%
%   fw_zoom('long', [303 312; 325 333; 342 350; 365 374], 'profundor')
%
% Cada coluna é uma janela. Linhas: deflexões, q (ou p e r), θ (ou φ), V.
% Serve para decidir se o comando é um doublet utilizável ou um pico isolado.

    if nargin < 3, tag = ''; end
    paths = setup_paths();
    outdir = fullfile(paths.images,'new_flights');  if ~exist(outdir,'dir'), mkdir(outdir); end
    F = fw_load(which_log);
    is_long = startsWith(lower(which_log),'lon');
    nW = size(wins,1);

    f = figure('Position',[30 30 340*nW+120 820],'Color','w'); try, f.Theme='light'; catch, end
    tl = tiledlayout(4, nW, 'TileSpacing','compact','Padding','compact');
    for c = 1:nW
        ii = F.t >= wins(c,1) & F.t <= wins(c,2);
        tt = F.t(ii);
        nexttile(c); hold on; grid on;
        plot(tt, F.de(ii), '-', 'Color',[0.85 0.37 0.01], 'LineWidth',1.3);
        plot(tt, F.da(ii), '-', 'Color',[0 0.45 0.7], 'LineWidth',0.9);
        plot(tt, F.dr(ii), '-', 'Color',[0.4 0.65 0.2], 'LineWidth',0.9);
        ylabel('deflexão'); xlim(wins(c,:)); ylim([-1.05 1.05]);
        if c == 1, legend({'profundor','aileron','leme'},'Location','southwest','FontSize',7); end
        text(0.02, 1.08, sprintf('%.0f a %.0f s', wins(c,:)), 'Units','normalized','FontWeight','bold');

        nexttile(nW + c); hold on; grid on;
        if is_long
            plot(tt, F.q(ii), 'k-', 'LineWidth',1.3); ylabel('q [rad/s]');
        else
            plot(tt, F.p(ii), 'k-', 'LineWidth',1.3); plot(tt, F.r(ii), '-','Color',[0.6 0.2 0.6],'LineWidth',1.1);
            ylabel('[rad/s]'); if c==1, legend({'p','r'},'FontSize',7,'Location','southwest'); end
        end
        xlim(wins(c,:));

        nexttile(2*nW + c); hold on; grid on;
        plot(tt, rad2deg(F.theta(ii)), '-','Color',[0.85 0.37 0.01],'LineWidth',1.2);
        plot(tt, rad2deg(F.phi(ii)), '-','Color',[0 0.45 0.7],'LineWidth',1.2);
        ylabel('[°]'); xlim(wins(c,:)); if c==1, legend({'\theta','\phi'},'FontSize',7,'Location','southwest'); end

        nexttile(3*nW + c); hold on; grid on;
        plot(tt, F.V(ii), 'k-', 'LineWidth',1.2); ylabel('V [m/s]'); xlabel('t [s]'); xlim(wins(c,:));
        yyaxis right; plot(tt, F.az(ii), '-','Color',[0.5 0.5 0.5]); ylabel('a_z [m/s²]');
    end
    title(tl, sprintf('%s — %s', strrep(F.name,'_','\_'), tag));
    fn = fullfile(outdir, sprintf('fw_%s_zoom.png', lower(which_log)));
    exportgraphics(f, fn, 'BackgroundColor','white','Resolution',120);
    fprintf('  Figura: %s\n', fn);
end
