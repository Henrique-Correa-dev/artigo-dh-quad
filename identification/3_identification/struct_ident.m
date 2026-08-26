% struct_ident.m — Identificabilidade estrutural LOCAL ao longo da trajetoria voada
% ================================================================================
% Teste: posto (numerico) da matriz de sensibilidade da saida do OEM aos 13
% parametros estimados, com DADOS SEM RUIDO (a medicao z some na derivada).
% Configuracao IDENTICA a da rodada oficial_2026_final: mesmas janelas, mesma
% grade, mesma cadeia de motor, multi-shooting de 2 s, mesmos pesos.
% Dois cenarios: observacao completa (taxas + acelerometro) e so taxas.
clear; clc;
IDROOT = '/Users/graest/ita-master/artigo/artigo-dh-quad/identification';
addpath(IDROOT); paths = setup_paths();
setappdata(0,'damp_form_override','moment'); clear aero_gains

SS = load(fullfile(paths.outputs,'runs','oficial_2026_teo','summary.mat'));
P  = SS.summary.P_final(:);
fprintf('P_final carregado: %d elementos | L_p=%.4f M_q=%.4f N_r=%.4f\n', numel(P), P(13), P(14), P(15));

t_trains = {[4,24]; [25,41]; [42,62]; [63,99]; [100,125]};
proj = parameters();  m = proj.m;  g = proj.g;

%% dados (identico a identify_plant secao 2)
L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t0 = max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]);
t1 = min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]);
t_common = t0:0.1:t1;  dt = 0.1;
ip = @(tt,xx) interp1(tt,xx,t_common,'linear');
gy = [ip(L.time_IMU,L.gyrX_raw)', ip(L.time_IMU,L.gyrY_raw)', ip(L.time_IMU,L.gyrZ_raw)'];
ac = [ip(L.time_IMU,L.accX_raw)', ip(L.time_IMU,L.accY_raw)', ip(L.time_IMU,L.accZ_raw)'];
at = [ip(L.time_ATT,L.roll_deg)', ip(L.time_ATT,L.pitch_deg)', ip(L.time_ATT,L.yaw_deg)'];
pw = [ip(L.time_RCOU,L.pwm1_raw)', ip(L.time_RCOU,L.pwm2_raw)', ...
      ip(L.time_RCOU,L.pwm3_raw)', ip(L.time_RCOU,L.pwm4_raw)'];
[rpm_lag, fT, fQ, mc] = motor_chain(t_common(:), pw);
VE = estimate_velocity(L, t_common(:));
Vg=VE.V; vg=VE.v; wg=VE.w;
Vg(~isfinite(Vg))=0; vg(~isfinite(vg))=0; wg(~isfinite(wg))=0;
setappdata(0,'aero_vel',struct('t',t_common(:),'V',Vg,'v',vg,'w',wg));

segs = cell(numel(t_trains),1);
for s=1:numel(t_trains)
    ii = t_common>=t_trains{s}(1) & t_common<=t_trains{s}(2);
    sg.time = t_common(ii)'; sg.N=sum(ii);
    sg.pwm = rpm_lag(ii,:); sg.pqr = gy(ii,:); sg.acc = ac(ii,:);
    sg.att_rad = deg2rad(at(ii,:));
    sg.V=Vg(ii); sg.v=vg(ii); sg.w=wg(ii);
    sg.T_ref=zeros(sg.N,4); sg.Q_ref=zeros(sg.N,4);
    for j=1:4, sg.T_ref(:,j)=fT(sg.pwm(:,j)); sg.Q_ref(:,j)=fQ(sg.pwm(:,j)); end
    segs{s}=sg;
end
pqr_cat = cell2mat(cellfun(@(s) s.pqr, segs,'UniformOutput',false));
acc_cat = cell2mat(cellfun(@(s) s.acc, segs,'UniformOutput',false));
wpqr = 1./max(var(pqr_cat),1e-12)';  wacc = 1./max(var(acc_cat),1e-12)';

%% sensibilidades
free  = 1:15;
names = {'Jx','Jy','Jz','Jxz','cT1','cT2','cT3','cT4','cM1','cM2','cM3','cM4','L_p','M_q','N_r'};
win_sec = 2;
modes = {'full','rotational'};
OUT = struct();
for mm=1:2
    mode = modes{mm};
    cf = @(PP) cost_all(PP, segs, win_sec, m, g, dt, wpqr, wacc, mode);
    e0 = cf(P);  Ne = numel(e0);
    Sm = zeros(Ne, numel(free));
    for i=1:numel(free)
        h = 1e-4*max(abs(P(free(i))),1e-3);
        Pp=P; Pp(free(i))=Pp(free(i))+h;
        Pn=P; Pn(free(i))=Pn(free(i))-h;
        Sm(:,i) = (cf(Pp)-cf(Pn))/(2*h);
    end
    % sensibilidade RELATIVA (coluna x |theta|): posto e correlacoes invariantes a escala
    Sr = Sm .* abs(P(free))';
    sv = svd(Sr);
    tolr = max(size(Sr))*eps(sv(1));
    nrank = sum(sv > tolr);
    F = Sr'*Sr;  C = pinv(F);
    sd = sqrt(diag(C));  Rho = C./(sd*sd');
    relsd = sd;                        % ruido unitario: so comparacoes relativas
    [~,imin] = min(sv);  [U2,S2,V2] = svd(Sr); vmin = V2(:,end);
    OUT.(mode) = struct('sv',sv,'cond',sv(1)/sv(end),'nrank',nrank, ...
                        'Rho',Rho,'relsd',relsd,'vmin',vmin,'Ne',Ne);
    fprintf('\n===== MODO %s (Ne=%d) =====\n', upper(mode), Ne);
    fprintf('valores singulares (norm. pelo 1o):\n');
    fprintf('  %.3e', sv/sv(1)); fprintf('\n');
    fprintf('posto numerico: %d de %d | condicionamento: %.3e\n', nrank, numel(free), sv(1)/sv(end));
    fprintf('%-5s %12s %12s\n','par','sd_rel(u.n.)','|vmin|');
    for i=1:numel(free)
        fprintf('%-5s %12.4g %12.3f\n', names{i}, relsd(i), abs(vmin(i)));
    end
    fprintf('pares |rho|>0.9:\n');
    for i=1:numel(free), for j=i+1:numel(free)
        if abs(Rho(i,j))>0.9, fprintf('  %s x %s: %+0.3f\n', names{i}, names{j}, Rho(i,j)); end
    end, end
end
save(fullfile('/Users/graest/ita-master/artigo/artigo-dh-quad/identification/outputs','struct_ident.mat'),'OUT','names','free','P');
fprintf('\nOK struct_ident\n');

function e = cost_all(P, segs, win_sec, m, g, dt, wpqr, wacc, mode)
    e = [];
    for s=1:numel(segs)
        sg = segs{s}; N = sg.N;
        wl = round(win_sec/dt);
        ws = 1:wl:N;
        if numel(ws)>1 && ws(end)==N, ws(end)=[]; end
        we = [ws(2:end)-1, N];
        aero_s = struct('V',sg.V(:),'v',sg.v(:),'w',sg.w(:));
        es = oem_ms_cost_func(P, sg.pqr, sg.acc, sg.att_rad, sg.T_ref, sg.Q_ref, ...
             m, g, dt, N, ws, we, wpqr, wacc, mode, aero_s);
        e = [e; es]; %#ok<AGROW>
    end
end
function e = oem_ms_cost_func(P, pqr, acc, att_rad, T_ref, Q_ref, m, g, dt, N, ...
    win_starts, win_ends, weights_pqr, weights_acc, cost_mode, aero)

    if nargin < 15, cost_mode = 'full'; end
    if nargin < 16, aero = []; end

    % Obter handles centralizados do vtol_dynamics (edite apenas lá!)
    dyn_h = vtol_dynamics('get_handles');
    trans_dot_fn = dyn_h.trans_dot;
    % Acelerômetro: arquivo separado (modelo de sensor, não de planta).
    % IMU offset r_imu lido de parameters().

    % Inércias → constantes G (corpo rígido, consistência garantida)
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
    gamma0 = Jx*Jz - Jxz^2;
    G1 = Jxz*(Jx - Jy + Jz) / gamma0;
    G2 = (Jz*(Jz - Jy) + Jxz^2) / gamma0;
    G3 = Jz / gamma0;
    G4 = Jxz / gamma0;
    G5 = (Jz - Jx) / Jy;
    G6 = Jxz / Jy;
    G7 = (Jx*(Jx - Jy) + Jxz^2) / gamma0;
    G8 = Jx / gamma0;
    invJy = 1 / Jy;

    k_T = P(5:8); k_Q = P(9:12);
    Dp = P(13); Dq = P(14); Dr = P(15);

    % Braços + IMU offset (de parameters() — fonte única)
    proj_p_oem = parameters();
    % FORMA DO AMORTECIMENTO (parameters().damp_form). Em 'moment' os P(13:15)
    % sao L_p, M_q, N_r [N.m.s] e entram DENTRO do vetor de momentos, junto com
    % a aerodinamica; em 'rate' sao c_p, c_q, c_r [1/s] subtraidos da derivada.
    Lpm = 0; Mqm = 0; Nrm = 0;  Dpr = Dp; Dqr = Dq; Drr = Dr;
    if isfield(proj_p_oem,'damp_form') && strcmp(proj_p_oem.damp_form,'moment')
        Lpm = Dp; Mqm = Dq; Nrm = Dr;  Dpr = 0; Dqr = 0; Drr = 0;
    end
    Lx_r = proj_p_oem.arms.Lx_r;
    Lx_l = proj_p_oem.arms.Lx_l;
    Ly_f = proj_p_oem.arms.Ly_f;
    Ly_r = proj_p_oem.arms.Ly_r;
    r_imu = proj_p_oem.imu_offset;

    % --- AERODINÂMICA DA ESTRUTURA (P(17:22), ∝ V) — coeficientes por amostra ---
    % aLp(k) = ¼ρSb²·Cl_p·V(k) multiplica o p SIMULADO (é amortecimento, entra
    % dentro do RK4); aLb(k) = ½ρSb·Cl_β·V(k)·v(k) é exógeno (parcela de
    % derrapagem). Idem para arfagem e guinada. V, v, w vêm de estimate_velocity.
    use_aero = ~isempty(aero) && numel(P) >= 22;
    if use_aero
        kA = aero_gains(proj_p_oem);
        Va = aero.V(:);  vb = aero.v(:);  wb = aero.w(:);
        Va(~isfinite(Va)) = 0;  vb(~isfinite(vb)) = 0;  wb(~isfinite(wb)) = 0;
        aLp = kA.Lp*P(17)*Va;   aLb = kA.Lb*P(18)*Va.*vb;
        aMq = kA.Mq*P(19)*Va;   aMa = kA.Ma*P(20)*Va.*wb;
        aNr = kA.Nr*P(21)*Va;   aNb = kA.Nb*P(22)*Va.*vb;
        if numel(P) >= 25
            if isappdata(0,'hu_form')   % forma do Hu: −½ρV²S·(b,c̄,b)·(c_l,c_m,c_n)
                qS = 0.5*kA.rho*Va.^2*kA.S;
                aLb = aLb - qS*kA.b*P(23);  aMa = aMa - qS*kA.c*P(24);  aNb = aNb - qS*kA.b*P(25);
            else                        % rotor: −L_v·v, −M_w·w, −N_v·v
                aLb = aLb - P(23)*vb;  aMa = aMa - P(24)*wb;  aNb = aNb - P(25)*vb;
            end
        end
    else
        aLp = zeros(N,1); aLb = aLp; aMq = aLp; aMa = aLp; aNr = aLp; aNb = aLp;
    end

    n_sub = 5;
    dt_sub = dt / n_sub;
    h2 = dt_sub / 2;
    MAX_VAL = 50;

    p_sim = zeros(N, 1);
    q_sim = zeros(N, 1);
    r_sim = zeros(N, 1);

    for w = 1:length(win_starts)
        i_s = win_starts(w);
        i_e = win_ends(w);
        ps = pqr(i_s,1); qs = pqr(i_s,2); rs = pqr(i_s,3);
        p_sim(i_s) = ps; q_sim(i_s) = qs; r_sim(i_s) = rs;

        for k = i_s:i_e-1
            % Momentos no início (k) e no fim (k+1) do passo — ArduPilot QuadX
            % (M1,M3=FRONT; M2,M4=REAR). O empuxo vem de motor_chain (RPM contínua
            % interpolada linearmente), então M(t) é interpolado LINEARMENTE dentro
            % do passo, igual ao que ode45/vtol_dynamics fazem na validação.
            % (Antes: M constante por 0,1 s → ≈50 ms de atraso efetivo a mais.)
            Tmr0 = k_T(:)' .* T_ref(k,:);    Qmr0 = k_Q(:)' .* Q_ref(k,:);
            Tmr1 = k_T(:)' .* T_ref(k+1,:);  Qmr1 = k_Q(:)' .* Q_ref(k+1,:);
            Mx0 = -(Lx_r*Tmr0(1) - Lx_l*Tmr0(2) - Lx_l*Tmr0(3) + Lx_r*Tmr0(4));
            My0 =  Ly_f*Tmr0(1) - Ly_r*Tmr0(2) + Ly_f*Tmr0(3) - Ly_r*Tmr0(4);
            Mz0 = Qmr0(1) + Qmr0(2) - Qmr0(3) - Qmr0(4);
            Mx1 = -(Lx_r*Tmr1(1) - Lx_l*Tmr1(2) - Lx_l*Tmr1(3) + Lx_r*Tmr1(4));
            My1 =  Ly_f*Tmr1(1) - Ly_r*Tmr1(2) + Ly_f*Tmr1(3) - Ly_r*Tmr1(4);
            Mz1 = Qmr1(1) + Qmr1(2) - Qmr1(3) - Qmr1(4);
            % parcelas aerodinâmicas EXÓGENAS (derrapagem/incidência) somam-se aos
            % momentos; as parcelas de AMORTECIMENTO (∝ V·p, V·q, V·r) dependem do
            % estado e entram estágio a estágio, logo abaixo.
            Mx0 = Mx0 + aLb(k);  Mx1 = Mx1 + aLb(k+1);
            My0 = My0 + aMa(k);  My1 = My1 + aMa(k+1);
            Mz0 = Mz0 + aNb(k);  Mz1 = Mz1 + aNb(k+1);
            dMx = Mx1-Mx0; dMy = My1-My0; dMz = Mz1-Mz0;
            aLp0 = aLp(k); daLp = aLp(k+1)-aLp(k);
            aMq0 = aMq(k); daMq = aMq(k+1)-aMq(k);
            aNr0 = aNr(k); daNr = aNr(k+1)-aNr(k);

            for si = 1:n_sub
                % RK4 sub-stepping com M(t) linear no passo: frações dos estágios
                fa = (si-1)/n_sub;  fb = fa + 0.5/n_sub;  fc = fa + 1/n_sub;
                Mxa = Mx0+fa*dMx; Mya = My0+fa*dMy; Mza = Mz0+fa*dMz;
                Mxb = Mx0+fb*dMx; Myb = My0+fb*dMy; Mzb = Mz0+fb*dMz;
                Mxc = Mx0+fc*dMx; Myc = My0+fc*dMy; Mzc = Mz0+fc*dMz;
                % coeficientes que multiplicam o estado DENTRO de M: aerodinamica
                % (aLp etc., ja negativa) menos o amortecimento em forma de momento
                La = aLp0+fa*daLp - Lpm; Lb_ = aLp0+fb*daLp - Lpm; Lc = aLp0+fc*daLp - Lpm;
                Qa = aMq0+fa*daMq - Mqm; Qb = aMq0+fb*daMq - Mqm; Qc = aMq0+fc*daMq - Mqm;
                Ra = aNr0+fa*daNr - Nrm; Rb = aNr0+fb*daNr - Nrm; Rc = aNr0+fc*daNr - Nrm;
                % k1
                MxA = Mxa + La*ps;  MyA = Mya + Qa*qs;  MzA = Mza + Ra*rs;
                pd1 = G1*ps*qs - G2*qs*rs + G3*MxA + G4*MzA - Dpr*ps;
                qd1 = G5*ps*rs - G6*(ps^2 - rs^2) + invJy*MyA - Dqr*qs;
                rd1 = G7*ps*qs - G1*qs*rs + G4*MxA + G8*MzA - Drr*rs;
                % k2
                p2 = ps+h2*pd1; q2 = qs+h2*qd1; r2 = rs+h2*rd1;
                MxB = Mxb + Lb_*p2;  MyB = Myb + Qb*q2;  MzB = Mzb + Rb*r2;
                pd2 = G1*p2*q2 - G2*q2*r2 + G3*MxB + G4*MzB - Dpr*p2;
                qd2 = G5*p2*r2 - G6*(p2^2 - r2^2) + invJy*MyB - Dqr*q2;
                rd2 = G7*p2*q2 - G1*q2*r2 + G4*MxB + G8*MzB - Drr*r2;
                % k3
                p3 = ps+h2*pd2; q3 = qs+h2*qd2; r3 = rs+h2*rd2;
                MxC = Mxb + Lb_*p3;  MyC = Myb + Qb*q3;  MzC = Mzb + Rb*r3;
                pd3 = G1*p3*q3 - G2*q3*r3 + G3*MxC + G4*MzC - Dpr*p3;
                qd3 = G5*p3*r3 - G6*(p3^2 - r3^2) + invJy*MyC - Dqr*q3;
                rd3 = G7*p3*q3 - G1*q3*r3 + G4*MxC + G8*MzC - Drr*r3;
                % k4
                p4 = ps+dt_sub*pd3; q4 = qs+dt_sub*qd3; r4 = rs+dt_sub*rd3;
                MxD = Mxc + Lc*p4;  MyD = Myc + Qc*q4;  MzD = Mzc + Rc*r4;
                pd4 = G1*p4*q4 - G2*q4*r4 + G3*MxD + G4*MzD - Dpr*p4;
                qd4 = G5*p4*r4 - G6*(p4^2 - r4^2) + invJy*MyD - Dqr*q4;
                rd4 = G7*p4*q4 - G1*q4*r4 + G4*MxD + G8*MzD - Drr*r4;
                % update
                ps = ps + dt_sub/6*(pd1 + 2*pd2 + 2*pd3 + pd4);
                qs = qs + dt_sub/6*(qd1 + 2*qd2 + 2*qd3 + qd4);
                rs = rs + dt_sub/6*(rd1 + 2*rd2 + 2*rd3 + rd4);
            end

            if ~isfinite(ps), ps = MAX_VAL; end
            if ~isfinite(qs), qs = MAX_VAL; end
            if ~isfinite(rs), rs = MAX_VAL; end
            ps = max(min(ps, MAX_VAL), -MAX_VAL);
            qs = max(min(qs, MAX_VAL), -MAX_VAL);
            rs = max(min(rs, MAX_VAL), -MAX_VAL);

            p_sim(k+1) = ps; q_sim(k+1) = qs; r_sim(k+1) = rs;
        end
    end

    sw_r = sqrt(weights_pqr(:));
    e_rot = [sw_r(1)*(pqr(:,1) - p_sim); ...
             sw_r(2)*(pqr(:,2) - q_sim); ...
             sw_r(3)*(pqr(:,3) - r_sim)];

    T_total = k_T(1)*T_ref(:,1) + k_T(2)*T_ref(:,2) + ...
              k_T(3)*T_ref(:,3) + k_T(4)*T_ref(:,4);
    % influxo no empuxo total: T ← T − 4·k_v·w (k_v medido, w da velocidade estimada)
    kvT = proj_p_oem.k_v_thrust;
    if kvT > 0 && use_aero, T_total = T_total - 4*kvT*wb; end

    phi_r = att_rad(:,1); theta_r = att_rad(:,2);
    gx = -g * sin(theta_r);
    gy = g * cos(theta_r) .* sin(phi_r);
    gz = g * cos(theta_r) .* cos(phi_r);

    % Integração translacional RK4 com sub-stepping (consistente c/ rotacional)
    n_sub_t = n_sub;        % mesmo número de sub-passos
    dt_sub_t = dt / n_sub_t;
    h2_t = dt_sub_t / 2;

    % Arrasto translacional por amostra, Nx3 [1/s]:
    %   x, y : arrasto induzido de rotor (g·C_d, constante) + estrutura ∝ V
    %   z    : só estrutura (o C_d de Beard age em x,y)
    kdm_P = [];
    if numel(P) >= 16, kdm_P = repmat(g*P(16)*[1 1 0], N, 1); end
    if numel(P) >= 25 && any(P(23:25) ~= 0) && use_aero
        kF = aero_gains(proj_p_oem);
        kdm_P = kdm_P + (kF.F/m) * (Va(:) * P(23:25)');
    end

    u_int = zeros(N,1); v_int = zeros(N,1); w_int = zeros(N,1);
    for k = 1:N-1
        pk = pqr(k,1); qk = pqr(k,2); rk = pqr(k,3);
        gxk = gx(k); gyk = gy(k); gzk = gz(k);
        Tk_m = T_total(k)/m;

        us = u_int(k); vs = v_int(k); ws = w_int(k);
        kdm_k = [];  if ~isempty(kdm_P), kdm_k = kdm_P(k,:); end
        for si = 1:n_sub_t
            % k1
            [ud1,vd1,wd1] = trans_dot_fn(pk,qk,rk, us,vs,ws, gxk,gyk,gzk, Tk_m, kdm_k);
            % k2
            u2=us+h2_t*ud1; v2=vs+h2_t*vd1; w2=ws+h2_t*wd1;
            [ud2,vd2,wd2] = trans_dot_fn(pk,qk,rk, u2,v2,w2, gxk,gyk,gzk, Tk_m, kdm_k);
            % k3
            u3=us+h2_t*ud2; v3=vs+h2_t*vd2; w3=ws+h2_t*wd2;
            [ud3,vd3,wd3] = trans_dot_fn(pk,qk,rk, u3,v3,w3, gxk,gyk,gzk, Tk_m, kdm_k);
            % k4
            u4=us+dt_sub_t*ud3; v4=vs+dt_sub_t*vd3; w4=ws+dt_sub_t*wd3;
            [ud4,vd4,wd4] = trans_dot_fn(pk,qk,rk, u4,v4,w4, gxk,gyk,gzk, Tk_m, kdm_k);
            % update
            us = us + dt_sub_t/6*(ud1 + 2*ud2 + 2*ud3 + ud4);
            vs = vs + dt_sub_t/6*(vd1 + 2*vd2 + 2*vd3 + vd4);
            ws = ws + dt_sub_t/6*(wd1 + 2*wd2 + 2*wd3 + wd4);
        end

        if ~isfinite(us), us = MAX_VAL; end
        if ~isfinite(vs), vs = MAX_VAL; end
        if ~isfinite(ws), ws = MAX_VAL; end
        us = max(min(us, MAX_VAL), -MAX_VAL);
        vs = max(min(vs, MAX_VAL), -MAX_VAL);
        ws = max(min(ws, MAX_VAL), -MAX_VAL);

        u_int(k+1) = us; v_int(k+1) = vs; w_int(k+1) = ws;
    end

    % Saída do modelo (acelerômetro como sensor — arquivo separado).
    % α via derivada numérica de p, q, r (no caso da cost: pqr medido).
    p_dot_sig = gradient(pqr(:,1), dt);
    q_dot_sig = gradient(pqr(:,2), dt);
    r_dot_sig = gradient(pqr(:,3), dt);
    [accX_m, accY_m, accZ_m] = accelerometer_model( ...
        pqr(:,1), pqr(:,2), pqr(:,3), ...
        u_int, v_int, w_int, ...
        T_total/m, ...
        p_dot_sig, q_dot_sig, r_dot_sig, ...
        r_imu, kdm_P);

    sw_a = sqrt(weights_acc(:));
    e_acc = [sw_a(1)*(acc(:,1) - accX_m); ...
             sw_a(2)*(acc(:,2) - accY_m); ...
             sw_a(3)*(acc(:,3) - accZ_m)];

    if strcmp(cost_mode, 'rotational')
        e = e_rot;
    elseif strcmp(cost_mode, 'translational')
        e = e_acc;
    else
        e = [e_rot; e_acc];
    end
end
