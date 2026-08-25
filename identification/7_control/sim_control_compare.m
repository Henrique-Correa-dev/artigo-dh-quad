% sim_control_compare.m — Comparação do controle: MODELO LINEAR vs NÃO-LINEAR
%
% Aplica a MESMA lei de controle (control_law) e a MESMA cadeia de atuação
% (alocação forces_to_pwm + SATURAÇÃO de PWM em [1200,2000] µs + curva real do
% motor) a DUAS plantas, isolando o efeito da linearização:
%   • LINEAR     : ẋ = A·x + B·(u_eff - u0),  ḣ = -w        (modelo do Cap.2/5)
%   • NÃO-LINEAR : corpo rígido 6-DOF (vtol_dynamics) + posição inercial NED
% A única diferença entre as curvas é a PLANTA — a atuação e a saturação são
% idênticas nos dois, conforme pedido.
%
% Dois cenários em escada (5 s por degrau) → 2 figuras p/ a §5.4:
%   1) Altitude : 2 → 4 → 6 → 4 → 2 → 0 m
%   2) Guinada  : 0 → 45 → 90 → 45 → 0 °   (em hover a 2 m)
%
% Uso:  >> design_control;  sim_control_compare

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

%% ===================== CARREGAR =====================
lm = load(fullfile(paths.outputs,'linear_model.mat'));
A=lm.A; B=lm.B; u0=lm.u0; P=lm.P;
G = load(fullfile(paths.control,'control_gains.mat')); K=G.K;
proj = parameters();  constants = struct('m',proj.m,'g',proj.g);
[fT,fQ] = motor_models();
dyn = vtol_dynamics('get_handles');
bridge = forces_to_pwm(lm,'nonlinear');

dt = 0.01;  PWM_LIM = [1200 2000];          % idle 1200 µs (rotor armado)
amber=[0.93 0.69 0.13]; graymed=[.5 .5 .5]; blue=[0 0.30 0.85];

%% ===================== CENÁRIOS (escadas, 5 s/degrau) =====================
ds = 5;
% 1) Altitude
t1=(0:dt:numel([2 4 6 4 2 0])*ds)';
sA = struct('name','altitude','h0',0,'t',t1, ...
    'sp_h',stair(t1,[2 4 6 4 2 0],ds),'sp_u',zeros(size(t1)),'sp_psi',zeros(size(t1)));
% 2) Guinada (hover a 2 m)
t2=(0:dt:numel([0 45 90 45 0])*ds)';
sB = struct('name','heading','h0',2,'t',t2, ...
    'sp_h',2*ones(size(t2)),'sp_u',zeros(size(t2)),'sp_psi',deg2rad(stair(t2,[0 45 90 45 0],ds)));
scen = {sA, sB};

%% ===================== RODA + PLOTA =====================
for si = 1:numel(scen)
    s = scen{si};
    RL = run_loop(false, s, A,B,u0,P,K,dyn,bridge,fT,fQ,constants,dt,PWM_LIM);
    RN = run_loop(true,  s, A,B,u0,P,K,dyn,bridge,fT,fQ,constants,dt,PWM_LIM);
    t  = s.t;

    isAlt = strcmp(s.name,'altitude');
    if isAlt
        yL=RL.H; yN=RN.H; spv=s.sp_h; ylab='Altitude (m)'; eun='(m)'; ttl='Altitude';
        scL=RL.U(:,1); scN=RN.U(:,1); sclab='Empuxo T (N)'; scttl='Empuxo total';
    else
        yL=rad2deg(RL.X(:,6)); yN=rad2deg(RN.X(:,6)); spv=rad2deg(s.sp_psi);
        ylab='Guinada \psi (°)'; eun='(°)'; ttl='Guinada (heading)';
        scL=rad2deg(RL.X(:,3)); scN=rad2deg(RN.X(:,3)); sclab='Taxa r (°/s)'; scttl='Taxa de guinada';
    end

    % --- métricas no log ---
    fprintf('==========================================================\n');
    fprintf('  CENÁRIO: escada de %s\n', s.name);
    fprintf('  %-10s | linear: erro máx=%6.3f, final=%7.3f | NL: erro máx=%6.3f, final=%7.3f\n', ...
        ttl, max(abs(spv-yL)), yL(end), max(abs(spv-yN)), yN(end));
    fprintf('  Desvio linear↔NL: máx=%.3f, RMS=%.3f %s | saturação: lin %.1f%%, NL %.1f%%\n', ...
        max(abs(yL-yN)), sqrt(mean((yL-yN).^2)), eun, 100*mean(RL.Usat), 100*mean(RN.Usat));

    % --- figura 2x2 ---
    fig=figure('Name',['sim_compare_' s.name],'Position',[60 60 1250 820]);
    set(fig,'Color','w'); try, fig.Theme='light'; catch, end

    subplot(2,2,1); hold on; grid on;
    plot(t,spv,'--','Color',amber,'LineWidth',1.8,'DisplayName','referência');
    plot(t,yL,'-','Color',graymed,'LineWidth',1.9,'DisplayName','linear');
    plot(t,yN,'-','Color',blue,'LineWidth',1.3,'DisplayName','não-linear');
    ylabel(ylab); legend('Location','best'); title(ttl);

    subplot(2,2,2); hold on; grid on;
    plot(t,spv-yL,'-','Color',graymed,'LineWidth',1.6,'DisplayName','linear');
    plot(t,spv-yN,'-','Color',blue,'LineWidth',1.3,'DisplayName','não-linear');
    yline(0,':','Color',[.7 .7 .7],'HandleVisibility','off');
    ylabel(['Erro de seguimento ' eun]); legend('Location','best'); title('Erro de seguimento');

    subplot(2,2,3); hold on; grid on;
    plot(t,scL,'-','Color',graymed,'LineWidth',1.6,'DisplayName','linear');
    plot(t,scN,'-','Color',blue,'LineWidth',1.3,'DisplayName','não-linear');
    ylabel(sclab); xlabel('t (s)'); legend('Location','best'); title(scttl);

    subplot(2,2,4); hold on; grid on;
    co=lines(4);
    for mi=1:4, plot(t,RN.PWM(:,mi),'Color',co(mi,:),'LineWidth',1.0,'DisplayName',sprintf('M%d',mi)); end
    yline(PWM_LIM(2),'r--','HandleVisibility','off'); yline(PWM_LIM(1),'r--','HandleVisibility','off');
    ylim([PWM_LIM(1)-160 PWM_LIM(2)+160]);
    ylabel('PWM (\mus)'); xlabel('t (s)'); legend('Location','best','Orientation','horizontal');
    title('PWM dos motores (não-linear, saturado)');

    sgtitle(sprintf('Controle em cascata — modelo LINEAR vs NÃO-LINEAR — escada de %s', s.name));
    out=fullfile(paths.images,['sim_compare_' s.name '.png']);
    exportgraphics(fig,out,'BackgroundColor','white');
    fprintf('  Figura: %s\n', out);
end
fprintf('==========================================================\n');


%% ===================== LOOP DE MALHA FECHADA (linear OU NL) =====================
function R = run_loop(isNL, s, A,B,u0,P,K,dyn,bridge,fT,fQ,constants,dt,PWM_LIM)
    N=numel(s.t);
    x=zeros(9,1); h=s.h0; pos=[0;0;-s.h0]; cs=struct('int_h',0,'int_u',0);
    X=zeros(N,9); H=zeros(N,1); U=zeros(N,4); PWM=zeros(N,4); Usat=false(N,1);
    for k=1:N
        if isNL, h=-pos(3); end
        sp=struct('h',s.sp_h(k),'u',s.sp_u(k),'psi',s.sp_psi(k));
        [u_cmd,cs]=control_law(x,h,sp,K,cs,dt);

        % --- alocação + SATURAÇÃO de PWM (idêntica nos dois modelos) ---
        pwm_b = bridge(u_cmd);                                  % u=[T,M] → PWM
        Usat(k)= any(pwm_b>=PWM_LIM(2)-1e-6 | pwm_b<=PWM_LIM(1)+1e-6);
        pwm   = min(max(pwm_b,PWM_LIM(1)),PWM_LIM(2));          % satura [1200,2000]
        [Tt,Mx,My,Mz]=dyn.forces(pwm',P,fT,fQ);                 % PWM saturado → forças
        u_eff=[Tt;Mx;My;Mz];

        X(k,:)=x'; H(k)=h; U(k,:)=u_eff'; PWM(k,:)=pwm';

        % --- integra a planta (RK4) ---
        if isNL
            fx=@(xx) dyn.rigid_body(xx,Tt,Mx,My,Mz,P,constants);
            fp=@(xx) body2ned_vel(xx);
            k1=fx(x);          p1=fp(x);
            k2=fx(x+dt/2*k1);  p2=fp(x+dt/2*k1);
            k3=fx(x+dt/2*k2);  p3=fp(x+dt/2*k2);
            k4=fx(x+dt*k3);    p4=fp(x+dt*k3);
            x   = x   + dt/6*(k1+2*k2+2*k3+k4);
            pos = pos + dt/6*(p1+2*p2+2*p3+p4);
        else
            f=@(xx,hh) deal(A*xx+B*(u_eff-u0), -xx(9));
            [k1x,k1h]=f(x,h);
            [k2x,k2h]=f(x+dt/2*k1x,h+dt/2*k1h);
            [k3x,k3h]=f(x+dt/2*k2x,h+dt/2*k2h);
            [k4x,k4h]=f(x+dt*k3x,  h+dt*k3h);
            x = x + dt/6*(k1x+2*k2x+2*k3x+k4x);
            h = h + dt/6*(k1h+2*k2h+2*k3h+k4h);
        end
    end
    R=struct('X',X,'H',H,'U',U,'PWM',PWM,'Usat',Usat);
end

%% ===================== HELPERS =====================
function sig = stair(t, levels, ds)
    sig = levels(min(floor(t(:)/ds)+1, numel(levels)))';
    sig = sig(:);
end

function vned = body2ned_vel(x)
    phi=x(4); th=x(5); ps=x(6); u=x(7); v=x(8); w=x(9);
    cphi=cos(phi); sphi=sin(phi); cth=cos(th); sth=sin(th); cps=cos(ps); sps=sin(ps);
    R = [ cth*cps,  sphi*sth*cps-cphi*sps,  cphi*sth*cps+sphi*sps; ...
          cth*sps,  sphi*sth*sps+cphi*cps,  cphi*sth*sps-sphi*cps; ...
         -sth,      sphi*cth,               cphi*cth ];
    vned = R*[u;v;w];
end
