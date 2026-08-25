function check_motor_map(logfile, tw, tag)
%CHECK_MOTOR_MAP  A fiação dos motores é a mesma que o modelo assume?
%
%   check_motor_map(logfile, [t0 t1], tag)
%
% Correlaciona os momentos calculados a partir do PWM com as acelerações
% angulares medidas, para todas as 24 permutações de numeração de motor e para
% os dois sentidos de rotação. A permutação certa dá correlação alta e positiva
% em rolagem e arfagem. Serve para descobrir, num log de outra campanha, se
% M1..M4 estão na mesma posição física do voo de referência.
%
% Modelo de referência (ArduPilot QuadX): M1=FR, M2=RL, M3=FL, M4=RR,
%   Mx = -(Lx_r·T1 - Lx_l·T2 - Lx_l·T3 + Lx_r·T4)
%   My =   Ly_f·T1 - Ly_r·T2 + Ly_f·T3 - Ly_r·T4
%   Mz =   Q1 + Q2 - Q3 - Q4        (M1,M2 CCW | M3,M4 CW)
    if nargin < 3, tag = 'log'; end
    addpath(fileparts(fileparts(mfilename('fullpath'))));
    paths = setup_paths();  pp = parameters();
    Lx_r = pp.arms.Lx_r; Lx_l = pp.arms.Lx_l; Ly_f = pp.arms.Ly_f; Ly_r = pp.arms.Ly_r;

    lp = logfile;  if ~exist(lp,'file'), lp = fullfile(paths.data, logfile); end
    L = load_log_data(lp);
    tg = (max([L.time_IMU(1) L.time_ATT(1) L.time_RCOU(1)]):0.1:min([L.time_IMU(end) L.time_ATT(end) L.time_RCOU(end)]))';
    if nargin >= 2 && ~isempty(tw), tg = tg(tg>=tw(1) & tg<=tw(2)); end
    ip = @(tt,xx) interp1(tt, xx, tg, 'linear');
    W = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ...
         ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
    pqr = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];
    dt = 0.1;
    pd = gradient(movmean(pqr(:,1),5), dt);
    qd = gradient(movmean(pqr(:,2),5), dt);
    rd = gradient(movmean(pqr(:,3),5), dt);

    [rpm, fT, fQ] = motor_chain(tg, W);
    T = fT(rpm);  Q = fQ(rpm);

    cc = @(a,b) local_corr(a,b);
    perms_all = perms(1:4);  perms_all = flipud(perms_all);
    best = struct('c',-Inf);
    fprintf('\n  %s — janela %.0f a %.0f s (%d amostras)\n', tag, tg(1), tg(end), numel(tg));
    fprintf('  ordem      corr(Mx,ṗ)  corr(My,q̇)  corr(Mz,ṙ)   soma rolagem+arfagem\n');
    for i = 1:size(perms_all,1)
        o = perms_all(i,:);
        Tp = T(:,o);  Qp = Q(:,o);
        Mx = -(Lx_r*Tp(:,1) - Lx_l*Tp(:,2) - Lx_l*Tp(:,3) + Lx_r*Tp(:,4));
        My =   Ly_f*Tp(:,1) - Ly_r*Tp(:,2) + Ly_f*Tp(:,3) - Ly_r*Tp(:,4);
        Mz =   Qp(:,1) + Qp(:,2) - Qp(:,3) - Qp(:,4);
        c = [cc(Mx,pd), cc(My,qd), cc(Mz,rd)];
        s = c(1) + c(2);
        if s > best.c, best = struct('c',s, 'o',o, 'cc',c); end
        if i <= 8 || s > 0.8
            fprintf('  [%d %d %d %d]   %+8.3f   %+9.3f   %+9.3f      %+6.3f\n', o, c, s);
        end
    end
    fprintf('  MELHOR: [%d %d %d %d]  →  %+.3f / %+.3f / %+.3f\n', best.o, best.cc);

    % assimetria média (indício de CG deslocado)
    Tm = mean(T,1);
    Mx_m = -(Lx_r*Tm(1) - Lx_l*Tm(2) - Lx_l*Tm(3) + Lx_r*Tm(4));
    My_m =   Ly_f*Tm(1) - Ly_r*Tm(2) + Ly_f*Tm(3) - Ly_r*Tm(4);
    Ttot = sum(Tm);
    fprintf('\n  Empuxos médios por motor [N]: %.2f %.2f %.2f %.2f   (total %.2f N, peso %.2f N)\n', ...
        Tm, Ttot, pp.m*pp.g);
    fprintf('  Momento médio residual: Mx %+0.3f N·m  My %+0.3f N·m\n', Mx_m, My_m);
    % Equilíbrio de momentos no pairado: o CG está onde a soma dá zero.
    %   dx = My/T  (dx > 0 → CG à frente)     dy = -Mx/T  (dy > 0 → CG à direita)
    dx = My_m/max(Ttot,1e-6);  dy = -Mx_m/max(Ttot,1e-6);
    fprintf('  CG estimado pelo equilíbrio: %+0.1f mm longitudinal (%s), %+0.1f mm lateral (%s)\n', ...
        1000*dx, ternary(dx>=0,'à frente','atrás'), 1000*dy, ternary(dy>=0,'à direita','à esquerda'));
end

function c = local_corr(a, b)
    a = a(:) - mean(a(:));  b = b(:) - mean(b(:));
    c = (a.'*b) / max(sqrt((a.'*a)*(b.'*b)), eps);
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
