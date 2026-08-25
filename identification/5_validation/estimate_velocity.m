function VE = estimate_velocity(L, tg, varargin)
%ESTIMATE_VELOCITY  Velocidade estimada com os sensores de bordo (sem GPS).
%
%   VE = estimate_velocity(L, tg)
%   VE = estimate_velocity(L, tg, 'anchor_s', 10, 'kT', 1)
%
%   Método (o mesmo de aero_influence.m, aferido no voo 4 contra o GPS: p95 e máx
%   reproduzidos, corr ≈ 0,3 no tempo → serve para envelope e para escalar
%   termos ∝ V, não como medida instantânea):
%     • vertical  w  : XKF1.VD do EKF (auxiliado pelo barômetro; válido sem GPS),
%                      lido dos logs brutos e levado ao eixo do concat;
%     • horizontal   : a_NED = R_n←b·[0;0;−T/m] + g·e3 com atitude medida e empuxo da
%                      cadeia de motor (motor_chain, k_T = 1), integrada com
%                      ancoragem v = 0 nas bordas de sub-janelas de anchor_s s.
%
%   Saída VE (todos na grade tg, NaN onde não há dado):
%     .VN .VE .VD   NED [m/s]      .Vh = hypot(VN,VE)      .V = hypot(Vh,VD)
%     .u .v .w      corpo (FRD)    .aN .aE  aceleração horizontal quase-estática
%
%   L: struct de load_log_data (concat ou log único). Para o VD do EKF, precisa
%   dos campos L.log_names / L.log_starts / L.boundaries (concat) ou L.source.

    opt = struct('anchor_s', 10, 'kT', 1);
    for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
    paths = setup_paths();  proj = parameters();  g = proj.g;  m = proj.m;
    tg = tg(:);  dt = median(diff(tg));  N = numel(tg);
    ip = @(tt,xx) interp1(tt, xx, tg, 'linear');

    phi = deg2rad(ip(L.time_ATT,L.roll_deg)); th = deg2rad(ip(L.time_ATT,L.pitch_deg)); psi = deg2rad(ip(L.time_ATT,L.yaw_deg));
    W4 = [ip(L.time_RCOU,L.pwm1_raw) ip(L.time_RCOU,L.pwm2_raw) ip(L.time_RCOU,L.pwm3_raw) ip(L.time_RCOU,L.pwm4_raw)];
    [rpm, fT] = motor_chain(tg, W4);
    T = sum(fT(rpm),2) * opt.kT;

    %% VD do EKF (barômetro) → eixo do concat
    VD = nan(N,1);
    if isfield(L,'log_names')
        names = L.log_names;  starts = L.log_starts;  bounds = L.boundaries;
    else
        names = {L.source};  starts = 0;  bounds = [];
    end
    for k = 1:numel(names)
        raw = names{k};  if ~exist(raw,'file'), raw = fullfile(paths.data, raw); end
        try
            if isfield(L,'log_names')
                Lr = load_log_data(raw);
                t0_log = min([Lr.time_IMU(1), Lr.time_ATT(1), Lr.time_RCOU(1)]);
                t_shift = starts(k) - t0_log;
            else
                t_shift = 0;
            end
            X = load(raw, 'XKF1_0').XKF1_0;
            tx = X(:,2)/1e6 + t_shift;  [tx, iu] = unique(tx);
            vd_k = interp1(tx, X(iu,9), tg, 'linear', NaN);
            if isfield(L,'log_names')
                sel = tg >= starts(k) & (k==numel(names) | tg <= bounds(min(k,end)));
            else
                sel = true(N,1);
            end
            VD(sel) = vd_k(sel);
        catch ME
            fprintf('  estimate_velocity: XKF1 do log %d indisponível (%s)\n', k, ME.message);
        end
    end

    %% aceleração horizontal quase-estática
    aN = zeros(N,1); aE = aN;
    for k = 1:N
        cph=cos(phi(k)); sph=sin(phi(k)); cth=cos(th(k)); sth=sin(th(k)); cps=cos(psi(k)); sps=sin(psi(k));
        e3n = [cph*sth*cps + sph*sps;  cph*sth*sps - sph*cps;  cph*cth];
        a = -(T(k)/m)*e3n + [0;0;g];
        aN(k) = a(1); aE(k) = a(2);
    end
    aN(~isfinite(aN)) = 0;  aE(~isfinite(aE)) = 0;

    %% integração ancorada em sub-janelas (em toda a grade)
    n_sub = max(1, round(opt.anchor_s/dt));
    VN = zeros(N,1); VEc = zeros(N,1);
    for s0 = 1:n_sub:N
        s1 = min(s0+n_sub-1, N);  ii = (s0:s1)';
        for c = 1:2
            if c==1, a = aN(ii); else, a = aE(ii); end
            vv = cumtrapz(tg(ii), a);  n = numel(vv);
            if n > 1, vv = vv - vv(1) - (vv(end)-vv(1))*((0:n-1)'/(n-1)); end
            if c==1, VN(ii) = vv; else, VEc(ii) = vv; end
        end
    end
    Vh = hypot(VN, VEc);

    %% corpo
    u = nan(N,1); v = u; w = u;
    for k = 1:N
        cph=cos(phi(k)); sph=sin(phi(k)); cth=cos(th(k)); sth=sin(th(k)); cps=cos(psi(k)); sps=sin(psi(k));
        Rbn = [ cth*cps, cth*sps, -sth; sph*sth*cps-cph*sps, sph*sth*sps+cph*cps, sph*cth; cph*sth*cps+sph*sps, cph*sth*sps-sph*cps, cph*cth ];
        vb = Rbn*[VN(k); VEc(k); VD(k)];  u(k)=vb(1); v(k)=vb(2); w(k)=vb(3);
    end

    VE = struct('t',tg,'VN',VN,'VE',VEc,'VD',VD,'Vh',Vh,'V',hypot(Vh, VD),'u',u,'v',v,'w',w,'aN',aN,'aE',aE, ...
                'anchor_s',opt.anchor_s);
end
