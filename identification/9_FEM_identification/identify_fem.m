% identify_fem.m — Identificação por FILTER-ERROR METHOD (taxa + aceleração)
% =========================================================================
% Jategaonkar (2015), Cap. 5 (Filter-Error Method) e §9.8. Do zero, sem toolbox.
%
% ESTADO  x = [p,q,r]                                  (3, dinâmico)
% OBS     y = [p,q,r, accX,accY,accZ]                  (6)
%   - rates  = x (medido direto, C=I nas 3 primeiras)
%   - accX/Y/Z = accelerometer_model(x, ẋ, T/m, r_imu) — ALGÉBRICA (não usa u,v,w)
%     → accZ ≈ −T/m ANCORA k_T → quebra a degenerescência k_T↔inércia
%       → permite liberar Jx,Jy (Jz,Jxz travados no CAD; full-mode do identify_plant)
% ENTRADA u = PWM(4) → [Mx,My,Mz] (taxas) e T/m (acc)
% θ = [Jx,Jy(livres), Jz,Jxz(travados), k_T(4), k_Q(4), Dp,Dq,Dr]
%
% ESTIMADOR (steady-state): ỹ=g(x̃,u) ; ν=z−ỹ ; x̂=x̃+K·ν_rate ; x̃(k+1)=x̂+∫f dt
%   (acc não corrige o estado — é algébrica; entra no CUSTO p/ ancorar k_T)
%   K via Riccati ESCALAR por eixo (taxa). σ_w=0→K=0→OEM.
%
% REGULARIZAÇÃO (= identify_plant): pares k_T/k_Q (yaw fraco) + triângulo inércia.
%
% VALIDAÇÃO: brancura das inovações (taxas) + previsão k-passos + R²/TIC em acc.
%
% Uso:  >> identify_fem
% =========================================================================

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));  paths = setup_paths();
img_dir = fullfile(here,'outputs','images'); if ~exist(img_dir,'dir'), mkdir(img_dir); end
out_dir = fullfile(here,'outputs');           if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% ===================== 0. CONFIG =====================
LOG_FILE = 'logs_concat.mat';
t_trains = {[4,24]; [25,41]; [42,62]; [63,99]; [100,125]};
t_val    = [605, 625];
dt       = 0.1;  SG_ORDER = 2;  SG_FRAME = 7;

SIGMA_V = [0.02; 0.02; 0.02];     % rad/s ruído de medida (gyro) — gain das taxas
SIGMA_W = [0.02; 0.05; 0.02];     % rad/s·√s ruído de processo POR EIXO [p,q,r]
                                  %   (pitch maior: é o canal mais colorido, ~0.65 Hz)
SIGMA_W_SWEEP = [];               % [] = run único; vetor p/ varrer (escala σ_w)
KSTEP_VAL = 10;                   % horizonte de previsão (1 s)

% k_Q: yaw só observa UMA escala (Mz=k_Q·(Q1+Q2−Q3−Q4)). Não dá p/ resolver 4.
SINGLE_KQ = true;                 % true → 1 k_Q compartilhado (identifiável); false → 4

% Jz, Jxz: por padrão travados no CAD (yaw fraco). Libere p/ DOCUMENTAR a tentativa:
% espera-se CRB alto (não-identificáveis por falta de excitação de yaw). Cap. 9.
FREE_JZ_JXZ = false;              % true → tenta identificar Jz,Jxz (provável CRB ruim)

% Remove o amortecimento (Dp=Dq=Dr=0): testa se o process noise do FEM dispensa o
% damping lumped (no rígido ideal não existe esse termo). Espera-se K maior.
REMOVE_D = true;                  % true → zera e trava Dp,Dq,Dr

REG.kt_pair = 100;  REG.kq_pair = 300;  REG.tri = 1e4;
N_RELAX = 3;

%% ===================== 1. CARREGAR + FILTRAR =====================
proj = parameters();  m = proj.m;  r_imu = proj.imu_offset;
[fT, fQ] = motor_models();
arms = [proj.arms.Lx_r, proj.arms.Lx_l, proj.arms.Ly_f, proj.arms.Ly_r];
dyn = vtol_dynamics('get_handles');  moments_fn = dyn.moments;

L = load_log_data(fullfile(paths.data, LOG_FILE));
t0 = max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]);
t1 = min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]);
tg = (t0:dt:t1)';
ip = @(t,x) interp1(t,x,tg,'linear');  sg = @(x) sgolayfilt(x,SG_ORDER,SG_FRAME);
f.p=sg(ip(L.time_IMU,L.gyrX_raw)); f.q=sg(ip(L.time_IMU,L.gyrY_raw)); f.r=sg(ip(L.time_IMU,L.gyrZ_raw));
f.ax=sg(ip(L.time_IMU,L.accX_raw)); f.ay=sg(ip(L.time_IMU,L.accY_raw)); f.az=sg(ip(L.time_IMU,L.accZ_raw));
f.w1=sg(ip(L.time_RCOU,L.pwm1_raw)); f.w2=sg(ip(L.time_RCOU,L.pwm2_raw));
f.w3=sg(ip(L.time_RCOU,L.pwm3_raw)); f.w4=sg(ip(L.time_RCOU,L.pwm4_raw));
segs = build_fem_segs(t_trains, tg, f, fT, fQ);
vseg = build_fem_segs({t_val}, tg, f, fT, fQ);  vseg = vseg{1};

% P0 + bounds: INÉRCIA via parameters.m (Jx,Jy livres ±25%; Jz,Jxz travados no CAD)
P0 = proj.P0_J;  lb = proj.bounds.lb;  ub = proj.bounds.ub;   % <-- NÃO trava 1:4
Pfile = fullfile(paths.outputs,'P_identified.mat');
if exist(Pfile,'file'), d=load(Pfile); P0=d.P_final; fprintf('warm-start: P_final (OEM)\n'); end

% SINGLE_KQ: usa só P(9) como k_Q de todos; trava P(10:12) (inertes) e some o kq_pair.
if SINGLE_KQ
    P0(10:12) = P0(9);  lb(10:12) = P0(9);  ub(10:12) = P0(9);  REG.kq_pair = 0;
end

% FREE_JZ_JXZ: libera Jz (±50% do CAD) e Jxz (livre, pequeno). Triângulo (e_tri) segura.
if FREE_JZ_JXZ
    lb(3) = 0.50*proj.J.Jz;  ub(3) = 1.50*proj.J.Jz;     % Jz ±50% do CAD
    lb(4) = -0.005;          ub(4) =  0.005;             % Jxz livre (±, pequeno)
    fprintf('  >>> Jz, Jxz LIBERADOS (tentativa — espera-se CRB alto p/ yaw fraco)\n');
end

% REMOVE_D: zera e trava o amortecimento (rígido ideal não tem esse termo).
if REMOVE_D
    P0(13:15) = 0;  lb(13:15) = 0;  ub(13:15) = 0;
    fprintf('  >>> Dp,Dq,Dr REMOVIDOS (=0) — process noise do FEM deve compensar\n');
end

pack = struct('segs',{segs},'arms',arms,'moments_fn',moments_fn,'dt',dt, ...
              'sv',SIGMA_V,'reg',REG,'lb',lb,'ub',ub,'nrelax',N_RELAX,'m',m, ...
              'r_imu',r_imu,'single_kq',SINGLE_KQ);

%% ===================== 2. (OPCIONAL) SWEEP de σ_w =====================
if ~isempty(SIGMA_W_SWEEP)
    fprintf('\n=== SWEEP de σ_w (brancura das taxas na validação) ===\n');
    ns=numel(SIGMA_W_SWEEP); Wht=zeros(ns,3); Dmat=zeros(ns,3);
    fprintf('  %8s | p q r (%%branco) | Dp Dq Dr\n','sigma_w');
    for is=1:ns
        sw=SIGMA_W_SWEEP(is);
        [Ps,Ks]=identify_once(P0,pack,sw,40);
        Wht(is,:)=whiteness_val(Ps,vseg,Ks,pack); Dmat(is,:)=Ps(13:15)';
        fprintf('  %8.3f | %3.0f %3.0f %3.0f | %.2f %.2f %.2f\n', sw, Wht(is,:), Ps(13),Ps(14),Ps(15));
    end
    figS=figure('Name','fem sigma sweep','Position',[60 60 1100 720]);
    subplot(2,1,1);hold on;grid on;plot(SIGMA_W_SWEEP,Wht,'-o');yline(90,'k--');ylabel('%branco');legend('p','q','r');title('Brancura vs σ_w');
    subplot(2,1,2);hold on;grid on;plot(SIGMA_W_SWEEP,Dmat,'-o');ylabel('D');xlabel('σ_w');legend('Dp','Dq','Dr');
    saveas(figS,fullfile(img_dir,'fem_sigma_sweep.png'));
end

%% ===================== 3. RUN ÚNICO =====================
fprintf('\n=== RUN FEM | σ_w=[%.3f %.3f %.3f] | obs=[pqr,acc] | Jx,Jy livres | SINGLE_KQ=%d ===\n', SIGMA_W(1),SIGMA_W(2),SIGMA_W(3),SINGLE_KQ);
[P, K] = identify_once(P0, pack, SIGMA_W, 80);
crb_fem(P, pack, K, ones(6,1), proj.param_names);   % usa P com 10:12 travado (bounds ok)
Pdisp = P; if SINGLE_KQ, Pdisp(10:12)=P(9); end     % k_Q único → iguala p/ display/save
print_params(Pdisp, proj.param_names);
tri_check(Pdisp);
P = Pdisp;                                          % daqui pra frente usa o coerente

%% ===================== 4. VALIDAÇÃO =====================
fprintf('\n==========================================================\n');
fprintf('  VALIDAÇÃO [%g-%gs] | K=[%.3f %.3f %.3f]\n', t_val(1),t_val(2),K(1),K(2),K(3));
fprintf('==========================================================\n');
[nu1,~]  = fem_filter(P, vseg, K, pack, 1);
[~,  yk] = fem_filter(P, vseg, K, pack, KSTEP_VAL);   % yk: 6 colunas [pqr,acc]
lbl={'p','q','r','accX','accY','accZ'};
fprintf('  %-5s | %%branco(ν) | prev %d-passos: R²     TIC\n','canal',KSTEP_VAL);
figv=figure('Name','fem validação','Position',[50 50 1300 820]);
for c=1:6
    z=vseg.z(:,c); fm=fit_metrics(z,yk(:,c));
    if c<=3, a=acf_(nu1(:,c),30); b=1.96/sqrt(numel(nu1(:,c))); pin=mean(abs(a(2:end))<=b)*100; ws=sprintf('%4.0f%%',pin);
    else, ws='  -- '; end
    fprintf('  %-5s |    %s   |             %6.3f %6.3f\n', lbl{c}, ws, fm.R2, fm.TIC);
    subplot(3,2,c); hold on; grid on;
    plot(vseg.t,z,'b','LineWidth',1.1,'DisplayName','medido');
    plot(vseg.t,yk(:,c),'r--','LineWidth',1.1,'DisplayName',sprintf('prev %d-passos',KSTEP_VAL));
    ylabel(lbl{c}); if c==1, legend('Location','best'); title(sprintf('Validação FEM [%g-%gs]',t_val(1),t_val(2))); end
    if c>=5, xlabel('t (s)'); end
end
saveas(figv,fullfile(img_dir,'fem_validation.png'));
fprintf('  Figura: %s\n', fullfile(img_dir,'fem_validation.png'));
save(fullfile(out_dir,'P_fem.mat'),'P','K','SIGMA_V','SIGMA_W','t_trains','t_val');
fprintf('\n  Salvo: %s\n', fullfile(out_dir,'P_fem.mat'));

%% ===================== FUNÇÕES LOCAIS =====================
function [P,K] = identify_once(P0, pk, sigma_w, maxit)
    sw3=sigma_w(:); if isscalar(sigma_w), sw3=sigma_w*ones(3,1); end
    W=ones(6,1); P=P0;
    opt=optimoptions('lsqnonlin','Display','off','MaxIterations',maxit, ...
        'MaxFunctionEvaluations',30000,'FunctionTolerance',1e-8,'StepTolerance',1e-10);
    for it=1:pk.nrelax
        K=steady_gain(P,pk.dt,pk.sv,sw3);
        cost=@(Pv) fem_cost(Pv,pk,K,W);
        P=lsqnonlin(cost,P,pk.lb,pk.ub,opt);
        nu=[]; for s=1:numel(pk.segs), nu=[nu; fem_filter(P,pk.segs{s},K,pk,1)]; end %#ok<AGROW>
        W=1./max(var(nu),1e-9)';
    end
end

function wv = whiteness_val(P, vseg, K, pk)
    nu=fem_filter(P,vseg,K,pk,1); wv=zeros(1,3);
    for c=1:3, a=acf_(nu(:,c),30); b=1.96/sqrt(numel(nu(:,c))); wv(c)=mean(abs(a(2:end))<=b)*100; end
end

function segs = build_fem_segs(tw, tg, f, fT, fQ)
    segs=cell(numel(tw),1);
    for s=1:numel(tw)
        mk=(tg>=tw{s}(1))&(tg<=tw{s}(2));
        d.t=tg(mk); d.z=[f.p(mk),f.q(mk),f.r(mk),f.ax(mk),f.ay(mk),f.az(mk)];
        pwm=[f.w1(mk),f.w2(mk),f.w3(mk),f.w4(mk)];
        d.T_ref=[fT(pwm(:,1)),fT(pwm(:,2)),fT(pwm(:,3)),fT(pwm(:,4))];
        d.Q_ref=[fQ(pwm(:,1)),fQ(pwm(:,2)),fQ(pwm(:,3)),fQ(pwm(:,4))];
        d.N=nnz(mk); segs{s}=d;
    end
end

function [G,D] = G_from_P(P)
    Jx=P(1);Jy=P(2);Jz=P(3);Jxz=P(4); g0=Jx*Jz-Jxz^2;
    G.G1=Jxz*(Jx-Jy+Jz)/g0; G.G2=(Jz*(Jz-Jy)+Jxz^2)/g0; G.G3=Jz/g0; G.G4=Jxz/g0;
    G.G5=(Jz-Jx)/Jy; G.G6=Jxz/Jy; G.G7=(Jx*(Jx-Jy)+Jxz^2)/g0; G.G8=Jx/g0; G.invJy=1/Jy;
    D=[P(13);P(14);P(15)];
end

function xd = rate_dot(x,M,G,D)
    p=x(1);q=x(2);r=x(3);
    xd=[ G.G1*p*q-G.G2*q*r+G.G3*M(1)+G.G4*M(3)-D(1)*p;
         G.G5*p*r-G.G6*(p^2-r^2)+G.invJy*M(2)-D(2)*q;
         G.G7*p*q-G.G1*q*r+G.G4*M(1)+G.G8*M(3)-D(3)*r ];
end

function [M,Tm] = forces_at(P,seg,k,pk)
    kq = P(9:12)';  if pk.single_kq, kq = P(9)*[1 1 1 1]; end   % 1 k_Q compartilhado
    Tmr=seg.T_ref(k,:).*P(5:8)'; Qmr=seg.Q_ref(k,:).*kq;
    [Mx,My,Mz]=moments_fn_call(pk.moments_fn,Tmr,Qmr,pk.arms); M=[Mx;My;Mz];
    Tm=sum(Tmr)/pk.m;
end
function [Mx,My,Mz]=moments_fn_call(mf,Tmr,Qmr,arms)
    [Mx,My,Mz]=mf(Tmr,Qmr,arms(1),arms(2),arms(3),arms(4));
end

function yacc = acc_out(x,xd,Tm,r_imu)
    [fx,fy,fz]=accelerometer_model(x(1),x(2),x(3),0,0,0,Tm,xd(1),xd(2),xd(3),r_imu);
    yacc=[fx;fy;fz];
end

function K = steady_gain(P,dt,sv,sw)
    [~,D]=G_from_P(P); K=zeros(3,1);
    for i=1:3
        phi=exp(-D(i)*dt); qd=sw(i)^2*dt; rr=sv(i)^2; Pp=qd; Mc=qd;
        for n=1:300, Mc=phi^2*Pp+qd; Kn=Mc/(Mc+rr); Pn=(1-Kn)*Mc; if abs(Pn-Pp)<1e-13,Pp=Pn;break;end; Pp=Pn; end
        K(i)=Mc/(Mc+rr);
    end
end

function [nu,yhist] = fem_filter(P,seg,K,pk,kstep)
% Estimador FEM. Observação 6-dim [pqr,acc]. Correção só nas TAXAS (acc é algébrica).
    if nargin<5, kstep=1; end
    [G,D]=G_from_P(P); N=seg.N; nu=zeros(N,6); yhist=zeros(N,6); xt=seg.z(1,1:3)';
    for k=1:N
        [M,Tm]=forces_at(P,seg,k,pk);
        xd=rate_dot(xt,M,G,D);
        yt=[xt; acc_out(xt,xd,Tm,pk.r_imu)];          % ỹ = [x̃; acc(x̃)]
        nuk=seg.z(k,:)'-yt; nu(k,:)=nuk'; yhist(k,:)=yt';
        if mod(k-1,kstep)==0, xh=xt+K.*nuk(1:3); else, xh=xt; end   % corrige só taxas
        if k<N
            xt=rk4(@(xx) rate_dot(xx,M,G,D), xh, pk.dt);
            xt=max(min(xt,1e3),-1e3);
        end
    end
end

function xn=rk4(f,x,h), k1=f(x);k2=f(x+h/2*k1);k3=f(x+h/2*k2);k4=f(x+h*k3); xn=x+h/6*(k1+2*k2+2*k3+k4); end

function e = fem_cost(P,pk,K,W)
    sw=sqrt(W(:)); e=[];
    for s=1:numel(pk.segs)
        nu=fem_filter(P,pk.segs{s},K,pk,1);
        nuw=nu.*sw'; e=[e; nuw(:)]; %#ok<AGROW>
    end
    R=pk.reg; kt=P(5:8); kq=P(9:12); Jx=P(1);Jy=P(2);Jz=P(3);
    e=[e; sqrt(R.kt_pair)*(kt(1)-kt(2)); sqrt(R.kt_pair)*(kt(3)-kt(4)); ...
          sqrt(R.kq_pair)*(kq(1)-kq(2)); sqrt(R.kq_pair)*(kq(3)-kq(4)); ...
          sqrt(R.tri)*max(0,Jz-(Jx+Jy)); sqrt(R.tri)*max(0,Jy-(Jx+Jz)); sqrt(R.tri)*max(0,Jx-(Jy+Jz))];
end

function print_params(P,names)
    fprintf('\n--- Parâmetros FEM ---\n');
    for i=1:15, fprintf('  %-6s : %10.5f\n', names{i}, P(i)); end
end

function tri_check(P)
    Jx=P(1);Jy=P(2);Jz=P(3);
    fprintf('\n--- Viabilidade física (triângulo) ---\n');
    chk(Jz,Jx+Jy,'Jz<=Jx+Jy'); chk(Jy,Jx+Jz,'Jy<=Jx+Jz'); chk(Jx,Jy+Jz,'Jx<=Jy+Jz');
end
function chk(a,b,nm)
    if a<=b, fprintf('   %-11s %.4f <= %.4f  OK\n',nm,a,b);
    else, fprintf('   %-11s %.4f  > %.4f  VIOLA <<<\n',nm,a,b); end
end

function crb_fem(P,pk,K,W,names)
    try
        cost=@(Pv) fem_cost(Pv,pk,K,W);
        o=optimoptions('lsqnonlin','Display','off','MaxIterations',0);
        [~,rn,res,~,~,~,J]=lsqnonlin(cost,P,pk.lb,pk.ub,o);
        free=find(pk.ub>pk.lb); Jf=full(J(:,free)); dof=max(numel(res)-numel(free),1);
        C=(rn/dof)*pinv(Jf'*Jf); se=sqrt(max(diag(C),0));
        fprintf('\n--- Qualidade (erro-padrão ≈ Cramér-Rao) ---\n');
        for i=1:numel(free)
            k=free(i); rel=100*se(i)/max(abs(P(k)),1e-9);
            tag='[ok]'; if rel>=20, tag=''; end; if rel>50, tag='<< fraco'; end
            fprintf('  %-6s %10.5f  %8.1f%%  %s\n', names{k}, P(k), rel, tag);
        end
    catch ME, fprintf('  (CRB: %s)\n', ME.message); end
end

function a=acf_(e,maxlag)
    e=e(:)-mean(e); N=numel(e); c0=sum(e.^2)+1e-12; a=zeros(maxlag+1,1);
    for k=0:maxlag, a(k+1)=sum(e(1:N-k).*e(1+k:N))/c0; end
end
