function W = fw_doublets(which_log, varargin)
%FW_DOUBLETS  Acha e mostra as janelas de doublet nos voos de asa fixa.
%
%   W = fw_doublets('long')   % profundor
%   W = fw_doublets('lat')    % aileron e leme
%
% Critério: a superfície sai do trim por mais de MIN_AMP (normalizado) e volta,
% dentro de um trecho com velocidade de pitot acima de V_MIN (voo estabelecido).
% Cada janela detectada é aberta com 1 s antes e 3 s depois, que é o tempo de a
% resposta livre decair.
%
% Gera outputs/images/new_flights/fw_<log>_doublets.png e devolve a matriz W
% com as janelas [t0 t1] e a superfície de cada uma.

    opt = struct('min_amp', 0.10, 'v_min', 12, 'pre', 1.0, 'post', 3.0, 'min_gap', 1.5);
    for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
    paths = setup_paths();
    outdir = fullfile(paths.images, 'new_flights');
    if ~exist(outdir,'dir'), mkdir(outdir); end

    F = fw_load(which_log);
    is_long = startsWith(lower(which_log), 'lon');
    if is_long, surf = {'de'};  nomes = {'profundor'};
    else,       surf = {'da','dr'};  nomes = {'aileron','leme'};
    end

    ok_v = F.V > opt.v_min;
    W = [];  src = {};
    for s = 1:numel(surf)
        u = F.(surf{s});
        u0 = movmedian(u, round(8/F.dt));          % trim lento (janela de 8 s)
        d  = u - u0;
        hit = abs(d) > opt.min_amp & ok_v;
        idx = find(hit);
        if isempty(idx), continue; end
        brk = [0; find(diff(idx) > round(opt.min_gap/F.dt)); numel(idx)];
        for k = 1:numel(brk)-1
            ii = idx(brk(k)+1 : brk(k+1));
            t0 = F.t(ii(1)) - opt.pre;  t1 = F.t(ii(end)) + opt.post;
            if t1 - t0 < 1.5, continue; end
            W(end+1,:) = [t0 t1]; %#ok<AGROW>
            src{end+1} = surf{s}; %#ok<AGROW>
        end
    end
    fprintf('\n  %s — %d janelas de excitação (V > %g m/s, amplitude > %.2f)\n', ...
        F.name, size(W,1), opt.v_min, opt.min_amp);
    fprintf('  %-4s %-10s %8s %8s %8s %8s %8s %8s\n', '#','superfície','t0','t1','V méd','|Δu|máx','std p','std q');
    for k = 1:size(W,1)
        ii = F.t >= W(k,1) & F.t <= W(k,2);
        u = F.(src{k});
        fprintf('  %-4d %-10s %8.1f %8.1f %8.1f %8.2f %8.3f %8.3f\n', k, src{k}, W(k,1), W(k,2), ...
            mean(F.V(ii)), max(abs(u(ii) - median(u(ii)))), std(F.p(ii)), std(F.q(ii)));
    end

    %% figura: visão geral do voo com as janelas marcadas
    f = figure('Position',[30 30 1350 900],'Color','w'); try, f.Theme='light'; catch, end
    tl = tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
    mark = @() arrayfun(@(k) xregion(W(k,1), W(k,2), 'FaceColor',[0.9 0.5 0.1], 'FaceAlpha',0.15), 1:size(W,1));

    nexttile; hold on; grid on;
    plot(F.t, F.V, 'k-', 'LineWidth',1); yline(opt.v_min,'--','V mínima considerada');
    mark(); ylabel('V pitot [m/s]'); ylim([0 35]);
    text(0.003, 1.06, sprintf('%s — janelas de doublet destacadas', strrep(F.name,'_','\_')), ...
        'Units','normalized','FontWeight','bold','FontSize',10);

    nexttile; hold on; grid on;
    plot(F.t, F.da, 'LineWidth',0.9); plot(F.t, F.de, 'LineWidth',0.9); plot(F.t, F.dr, 'LineWidth',0.9);
    mark(); ylabel('deflexão [-1,1]'); legend({'aileron','profundor','leme'},'Location','eastoutside');

    nexttile; hold on; grid on;
    plot(F.t, F.p, 'LineWidth',0.8); plot(F.t, F.q, 'LineWidth',0.8); plot(F.t, F.r, 'LineWidth',0.8);
    mark(); ylabel('[rad/s]'); legend({'p','q','r'},'Location','eastoutside');

    nexttile; hold on; grid on;
    plot(F.t, rad2deg(F.phi), 'LineWidth',0.9); plot(F.t, rad2deg(F.theta), 'LineWidth',0.9);
    mark(); ylabel('[°]'); legend({'\phi','\theta'},'Location','eastoutside');

    nexttile; hold on; grid on;
    plot(F.t, F.alt, 'k-', 'LineWidth',1); mark();
    ylabel('alt baro [m]'); xlabel('t [s]');
    yyaxis right; plot(F.t, F.thr, '-', 'Color',[0.4 0.4 0.4]); ylabel('motor [0,1]');

    fn = fullfile(outdir, sprintf('fw_%s_doublets.png', lower(which_log)));
    exportgraphics(f, fn, 'BackgroundColor','white','Resolution',120);
    fprintf('  Figura: %s\n', fn);
end
