function dydt = vtol_dynamics(t, y, P, pwm_time, pwm_signals, func_T_ref, func_Q_ref, constants)
    % VTOL_DYNAMICS  Dinâmica do drone — separada em ALOCAÇÃO + CORPO RÍGIDO.
    %
    % ARQUITETURA (corte limpo):
    %   ALOCAÇÃO (motores)        CORPO RÍGIDO (6-DOF)
    %   PWM → [T, Mx, My, Mz]  →  [T, Mx, My, Mz] → ẋ
    %   (usa fT, fQ, k, braços)   (G's, 1/Jy, 1/m, gravidade)
    %
    % A interface EXTERNA continua sendo PWM (pra sysid/validação que reproduz
    % log). Internamente: forces_from_pwm_local (alocação) → rigid_body_local.
    %
    % DISPATCH ESPECIAL:
    %   dyn_h = vtol_dynamics('get_handles')
    %       Retorna struct com handles pras subfunções (fonte única):
    %         .moments    = @moments_local(Tmr, Qmr, Lx_r, Lx_l, Ly_f, Ly_r)
    %         .trans_dot  = @trans_dot_local(p,q,r, u,v,w, gx,gy,gz, T_m)
    %         .rigid_body = @rigid_body_local(y, T, Mx, My, Mz, P, constants)
    %             → CORPO RÍGIDO: recebe forças generalizadas, devolve ẋ.
    %               Usado por linearize_at (entrada [T,M]) e por quem fechar malha.
    %         .forces     = @forces_from_pwm_local(pwm, P, fT, fQ)
    %             → ALOCAÇÃO direta: PWM → [T, Mx, My, Mz]. Vetorizado (aceita Nx4).
    %
    % MODOS NORMAIS (dispatch por tamanho de estado):
    %   length(y) == 3   -> rotacional puro             [p; q; r]
    %   length(y) == 9   -> completo                    [p..r; phi..psi; u..w]
    %
    % P (vetor de 15 elementos — PURAMENTE ROTACIONAL + MOTOR):
    %   P(1:4)   = [Jx, Jy, Jz, Jxz]
    %   P(5:8)   = k_T1..k_T4
    %   P(9:12)  = k_Q1..k_Q4
    %   P(13:15) = Dp, Dq, Dr
    %
    % constants: struct com .m, .g.

    %% =========================================================================
    %  Dispatch especial: retorna handles pras subfunções
    %  =========================================================================
    if ischar(t) && strcmp(t, 'get_handles')
        dydt = struct();
        dydt.moments    = @moments_local;
        dydt.trans_dot  = @trans_dot_local;
        dydt.rigid_body = @rigid_body_local;       % [T,Mx,My,Mz] → ẋ
        dydt.forces     = @forces_from_pwm_local;  % PWM → [T,Mx,My,Mz] (alocação)
        return;
    end

    %% =========================================================================
    %  ALOCAÇÃO: PWM interpolado em t → [T, Mx, My, Mz]
    %  =========================================================================
    current_pwm = zeros(1,4);
    for i = 1:4
        current_pwm(i) = interp1(pwm_time, pwm_signals(:,i), t, 'linear', 'extrap');
    end
    [T_total, Mx, My, Mz] = forces_from_pwm_local(current_pwm, P, func_T_ref, func_Q_ref);

    %% =========================================================================
    %  TESTE: amortecimento escalando com a rotação dos rotores (Nguyen & Webb)
    %  =========================================================================
    %  No modelo da NASA cada rotor contribui ρπR⁵·Ω_i·C_lp,i, somado sobre os
    %  rotores. No balanço físico o influxo faz k_v ∝ √T ∝ Ω. Então, em vez de
    %  L_p constante, L_p(t) = L_p,id · [ΣΩ_i(t)/ΣΩ_i,ref]^n, com Ω_i a RPM
    %  instantânea da cadeia de motor (current_pwm aqui já é RPM) e Ω_ref a do
    %  pairado. Sem parâmetro novo. n = 1 (∝ Ω, forma da NASA) ou n = 2 (∝ T).
    %  Ativado por appdata 'damp_omega' = struct('n', 1, 'rpm_ref', RPM_hover).
    dO = getappdata(0, 'damp_omega');
    if ~isempty(dO) && isstruct(dO)
        ratio = max(sum(current_pwm), 1) / max(4*dO.rpm_ref, 1);
        P(13:15) = P(13:15) * ratio^dO.n;
    end

    %% =========================================================================
    %  AERODINÂMICA DA ESTRUTURA — P(17:22), forma regular ∝ V (aero_gains.m)
    %  =========================================================================
    %  A velocidade que escala estes termos vem, em ordem de preferência:
    %    [1] appdata 'aero_vel' (struct .t .V .v .w) — velocidade ESTIMADA de
    %        bordo (estimate_velocity), exógena. É a mesma que a identificação
    %        usa no custo, então validação e identificação veem o mesmo V.
    %    [2] os próprios estados u, v, w quando y tem 9 estados (simulação
    %        autocontida; atenção: u,v,w integrados derivam sem GPS).
    %    [3] zero (modo rotacional de 3 estados sem appdata) → sem aerodinâmica.
    if numel(P) >= 22 && (any(P(17:22) ~= 0) || (numel(P) >= 25 && any(P(23:25) ~= 0)) || isappdata(0,'faero_on'))
        av = getappdata(0, 'aero_vel');
        if ~isempty(av) && isstruct(av)
            Va = interp1(av.t, av.V, t, 'linear', 'extrap');
            va = interp1(av.t, av.v, t, 'linear', 'extrap');
            wa = interp1(av.t, av.w, t, 'linear', 'extrap');
        else
            % SEM a velocidade exógena a aerodinâmica fica DESLIGADA, e não
            % caindo nos estados u,v,w. Os estados translacionais divergem em
            % malha aberta (u chega a −142 m/s em 20 s), e usá-los aqui
            % multiplica os coeficientes por uma velocidade cem vezes maior que
            % a real, o que arruína a simulação silenciosamente. Quem simula com
            % P(17:22) ≠ 0 TEM de chamar set_aero_vel primeiro.
            Va = 0; va = 0; wa = 0;
            warn_once_no_aero_vel();
        end
        if ~isfinite(Va), Va = 0; end
        if ~isfinite(va), va = 0; end
        if ~isfinite(wa), wa = 0; end
        kA = aero_gains();
        setappdata(0, 'aero_V_now', Va);     % lido pelo rigid_body_local (forças)
        ua = 0; if isfield(av,'u'), ua = interp1(av.t, av.u, t, 'linear', 'extrap'); end
        if ~isfinite(ua), ua = 0; end
        setappdata(0, 'aero_vel_now', [Va va wa ua]);
        Mx = Mx + kA.Lp*P(17)*Va*y(1) + kA.Lb*P(18)*Va*va;
        My = My + kA.Mq*P(19)*Va*y(2) + kA.Ma*P(20)*Va*wa;
        Mz = Mz + kA.Nr*P(21)*Va*y(3) + kA.Nb*P(22)*Va*va;
        %  Termos fixos de equivalência estrutural com os executáveis de asa
        %  fixa (parameters().wing): trim de arfagem C_m0 e cruzadas C_lr, C_np,
        %  ativos apenas quando o eixo correspondente tem aerodinâmica ligada.
        wg = parameters().wing;
        if P(19) ~= 0 || P(20) ~= 0
            My = My + 0.5*kA.rho*kA.S*kA.c*wg.C_m0*Va^2;
        end
        if P(17) ~= 0 || P(18) ~= 0
            Mx = Mx + kA.Lp*wg.C_lr*Va*y(3);     % ¼ρSb²·C_lr·V·r
        end
        if P(21) ~= 0 || P(22) ~= 0
            Mz = Mz + kA.Nr*wg.C_np*Va*y(1);     % ¼ρSb²·C_np·V·p
        end
        if numel(P) >= 25
            if isappdata(0,'hu_form')    % forma do Hu: −½ρV²S·(b,c̄,b)·(c_l,c_m,c_n), um coef por eixo
                qS = 0.5*kA.rho*Va^2*kA.S;
                Mx = Mx - qS*kA.b*P(23);  My = My - qS*kA.c*P(24);  Mz = Mz - qS*kA.b*P(25);
            else                          % rotor: momento por velocidade (NASA C_lv, C_mw, C_nv)
                Mx = Mx - P(23)*va;  My = My - P(24)*wa;  Mz = Mz - P(25)*va;
            end
        end
    end

    %% =========================================================================
    %  INFLUXO NO EMPUXO TOTAL: T ← T − 4·k_v·w   (k_v medido, w exógeno)
    %  =========================================================================
    kvT = parameters().k_v_thrust;
    if kvT > 0
        avn = getappdata(0,'aero_vel_now');
        if isempty(avn)
            av2 = getappdata(0,'aero_vel');
            if ~isempty(av2) && isstruct(av2), avn = [0 0 interp1(av2.t, av2.w, t, 'linear', 0) 0]; end
        end
        if ~isempty(avn) && isfinite(avn(3)), T_total = T_total - 4*kvT*avn(3); end
    end

    %% =========================================================================
    %  CORPO RÍGIDO: [T, Mx, My, Mz] → ẋ
    %  =========================================================================
    % --- TESTE DE DIAGNÓSTICO (não afeta o uso normal): amortecimento ∝ V ---
    %  Se existir appdata 'diag_damp' (posto por identify_plant_diag com
    %  DIAG.damp_mode = 'aero' | 'hybrid'), o termo −Dp·p vira −(c0 + Dp·V(t))·p,
    %  com V(t) a velocidade estimada (estimate_velocity) e c0 o a priori de rotor
    %  (0 em 'aero'). Sem appdata, P(13:15) é usado como sempre.
    dd = getappdata(0,'diag_damp');
    if ~isempty(dd) && isstruct(dd) && any(strcmp(dd.mode,{'aero','hybrid'}))
        Vt = interp1(dd.t, dd.V, t, 'linear', 'extrap');
        P(13:15) = dd.c0(:) + P(13:15).*Vt;
    end
    dydt = rigid_body_local(y, T_total, Mx, My, Mz, P, constants);
end

%% =========================================================================
%  ALOCAÇÃO: PWM → forças generalizadas [T_total, Mx, My, Mz]
%    Vetorizado: aceita pwm 1x4 (escalar no tempo) OU Nx4 (série) → saída casada.
%    Esta é a peça NÃO-LINEAR (curvas fT, fQ + k por motor). Fora do corpo rígido.
%  =========================================================================
function [T_total, Mx, My, Mz] = forces_from_pwm_local(pwm, P, func_T_ref, func_Q_ref)
    k_T = P(5:8);
    k_Q = P(9:12);

    proj_p = parameters();
    Lx_r = proj_p.arms.Lx_r;
    Lx_l = proj_p.arms.Lx_l;
    Ly_f = proj_p.arms.Ly_f;
    Ly_r = proj_p.arms.Ly_r;

    % func_T_ref/func_Q_ref operam elemento-a-elemento (interp1) → aceitam Nx4.
    Tmr = k_T(:)' .* func_T_ref(pwm);   % 1x4 ou Nx4
    Qmr = k_Q(:)' .* func_Q_ref(pwm);

    % moments_local detecta Nx4 e devolve coluna(s); 1x4 → escalares.
    [Mx, My, Mz] = moments_local(Tmr, Qmr, Lx_r, Lx_l, Ly_f, Ly_r);
    T_total = sum(Tmr, 2);
end

%% =========================================================================
%  CORPO RÍGIDO: forças generalizadas [T, Mx, My, Mz] → ẋ
%    LINEAR em [T, Mx, My, Mz] (entram via G3·Mx+G4·Mz, invJy·My, -T/m).
%    As não-linearidades aqui são só de ESTADO (giroscópico pq, gravidade sinθ).
%    Dispatch por tamanho de y: 3 (rotacional) ou 9 (completo).
%  =========================================================================
function dydt = rigid_body_local(y, T_total, Mx, My, Mz, P, constants)
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
    Dp = P(13); Dq = P(14); Dr = P(15);
    % FORMA DO AMORTECIMENTO (parameters().damp_form):
    %  'moment' → P(13:15) = L_p, M_q, N_r [N·m·s], entram no VETOR DE MOMENTOS,
    %             como manda a equação de corpo rígido. Traz os termos cruzados
    %             via J_xz de graça (ṗ recebe −Γ4·N_r·r, ṙ recebe −Γ4·L_p·p).
    %  'rate'   → P(13:15) = c_p, c_q, c_r [1/s], subtraídos da derivada (legado).
    dform = 'moment';
    pp_ = parameters();  if isfield(pp_,'damp_form'), dform = pp_.damp_form; end
    is_moment = strcmp(dform,'moment');
    if is_moment
        Mx = Mx - Dp*y(1);  My = My - Dq*y(2);  Mz = Mz - Dr*y(3);
        Dp = 0; Dq = 0; Dr = 0;
    end

    % Constantes G do corpo rígido (Beard-McLain)
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

    p = y(1); q = y(2); r = y(3);

    %% Acelerações angulares (rotacional)
    p_dot = G1*p*q - G2*q*r + G3*Mx + G4*Mz - Dp*p;
    q_dot = G5*p*r - G6*(p^2 - r^2) + invJy*My - Dq*q;
    r_dot = G7*p*q - G1*q*r + G4*Mx + G8*Mz - Dr*r;

    if length(y) == 3
        dydt = [p_dot; q_dot; r_dot];
        return;
    end

    %% Modo 9 estados: cinemática de atitude + translacional
    phi   = y(4);
    theta = y(5);
    % psi   = y(6);   % não usado nas equações (só sai como integral)
    u = y(7); v = y(8); w = y(9);

    if nargin >= 7 && ~isempty(constants)
        m_body = constants.m;
        g_acc  = constants.g;
    else
        m_body = 1.6011;
        g_acc  = 9.81;
    end

    cos_theta = cos(theta);
    if abs(cos_theta) < 1e-7
        cos_theta = 1e-7 * sign(cos_theta + 1e-12);
    end
    sin_theta = sin(theta);
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    tan_theta = sin_theta / cos_theta;

    %% Cinemática de Euler (taxa de atitude)
    phi_dot   = p + (q*sin_phi + r*cos_phi) * tan_theta;
    theta_dot = q*cos_phi - r*sin_phi;
    psi_dot   = (q*sin_phi + r*cos_phi) / cos_theta;

    %% Gravidade no body frame
    gx_body = -g_acc * sin_theta;
    gy_body =  g_acc * sin_phi * cos_theta;
    gz_body =  g_acc * cos_phi * cos_theta;

    %% Acelerações translacionais (via subfunção — fonte única)
    T_m = T_total / m_body;
    % Arrasto translacional: rotor (C_d, constante) + estrutura (∝ V, P(23:25)).
    % Va vem por variável global posta pelo vtol_dynamics, que é quem tem o t.
    kdm_P = [];
    if numel(P) >= 16, kdm_P = g_acc * P(16) * [1;1;0]; end
    [u_dot, v_dot, w_dot] = trans_dot_local(p, q, r, u, v, w, ...
                                            gx_body, gy_body, gz_body, T_m, kdm_P);
    % F_aero da estrutura (Lombaerts 2020, Eqs. 5, 7, 9): placa plana com
    % C_p = 2, em α e β da velocidade ESTIMADA de bordo (exógena, appdata
    % 'aero_vel'), nunca dos estados u,v,w integrados, que divergem.
    % Ativada por appdata 'faero_on'. Força específica somada a u̇, v̇, ẇ.
    if isappdata(0,'faero_on')
        if strcmp(getappdata(0,'faero_on'),'measured'), fa = faero_measured_local();
        else,                                           fa = faero_flatplate_local(); end
        u_dot = u_dot + fa(1)/m_body;  v_dot = v_dot + fa(2)/m_body;  w_dot = w_dot + fa(3)/m_body;
    end

    %% Monta vetor de saída (9 estados)
    dydt = [p_dot; q_dot; r_dot; ...
            phi_dot; theta_dot; psi_dot; ...
            u_dot; v_dot; w_dot];
end

%% =========================================================================
%  SUBFUNÇÃO: Momentos no body frame (ArduPilot QuadX padrão)
%    M1=FR, M2=RL, M3=FL, M4=RR   |   M1+M2 CCW, M3+M4 CW
%  =========================================================================
function [Mx, My, Mz] = moments_local(Tmr, Qmr, Lx_r, Lx_l, Ly_f, Ly_r)
    % Aceita Tmr, Qmr como vetor 4x1, 1x4 ou Nx4 (vetorizado).
    if size(Tmr, 2) == 4   % formato Nx4 (vetorizado)
        Mx = -(Lx_r*Tmr(:,1) - Lx_l*Tmr(:,2) - Lx_l*Tmr(:,3) + Lx_r*Tmr(:,4));
        My =   Ly_f*Tmr(:,1) - Ly_r*Tmr(:,2) + Ly_f*Tmr(:,3) - Ly_r*Tmr(:,4);
        Mz =   Qmr(:,1) + Qmr(:,2) - Qmr(:,3) - Qmr(:,4);
    else                    % formato escalar (4x1 ou 1x4)
        Mx = -(Lx_r*Tmr(1) - Lx_l*Tmr(2) - Lx_l*Tmr(3) + Lx_r*Tmr(4));
        My =   Ly_f*Tmr(1) - Ly_r*Tmr(2) + Ly_f*Tmr(3) - Ly_r*Tmr(4);
        Mz =   Qmr(1) + Qmr(2) - Qmr(3) - Qmr(4);
    end
end

%% =========================================================================
%  SUBFUNÇÃO: Acelerações translacionais (sem drag)
%    Drag (Xu, Yv, Zw) removido — não-identificável sem GPS velocity.
%    Aceita escalares OU vetores Nx1 (broadcast por .*).
%  =========================================================================
function [u_dot, v_dot, w_dot] = trans_dot_local(p, q, r, u, v, w, gx, gy, gz, T_m, kdm_in)
    % Arrasto induzido do rotor (Beard & McLain, cap. 14.4.2): f_d ≈ −T·C_d·[u;v;0].
    % Com T ≈ mg (pairado), a força específica vira −g·C_d·[u;v;0].
    % Fonte, em ordem: argumento kdm_in (= g·C_d, vindo de P(16)) → appdata
    % 'diag_drag' (teste antigo do identify_plant_diag) → parameters().k_drag.
    if nargin >= 11 && ~isempty(kdm_in)
        kdm = kdm_in;
    else
        pp = parameters();  kdm = pp.k_drag / pp.m;     % k_drag/m
        cd_ = getappdata(0,'diag_drag');
        if ~isempty(cd_), kdm = pp.g * cd_; end         % C_d adimensional → g·C_d [1/s]
    end
    % kdm pode ser ESCALAR (só arrasto de rotor, igual em x e y, nada em z) ou
    % VETOR de 3 (x, y, z), que é o caso quando as forças aerodinâmicas da
    % estrutura P(23:25) estão ativas e cada eixo tem seu próprio amortecimento.
    if size(kdm,2) >= 3        % Nx3 (um por amostra) ou 1x3
        kx = kdm(:,1); ky = kdm(:,2); kz = kdm(:,3);
    elseif numel(kdm) >= 3     % vetor coluna de 3
        kx = kdm(1);   ky = kdm(2);   kz = kdm(3);
    else                        % escalar: só arrasto de rotor, nada em z
        kx = kdm;      ky = kdm;      kz = 0;
    end
    u_dot = r.*v - q.*w + gx - kx.*u;
    v_dot = p.*w - r.*u + gy - ky.*v;
    w_dot = q.*u - p.*v - T_m + gz - kz.*w;
end

function warn_once_no_aero_vel()
%WARN_ONCE_NO_AERO_VEL  Avisa uma vez por sessão que a aerodinâmica está inerte.
    persistent avisado
    if isempty(avisado)
        warning('vtol_dynamics:semAeroVel', ...
            ['P(17:22) nao nulos mas o appdata ''aero_vel'' nao foi posto: ' ...
             'os momentos aerodinamicos ficam DESLIGADOS. Chame set_aero_vel(L, tg).']);
        avisado = true;
    end
end

function fa = faero_measured_local()
%FAERO_MEASURED_LOCAL  Força aerodinâmica da asa com os coeficientes MEDIDOS em voo.
%   Mesma forma regular dos momentos (q̄ substituído, razões abertas):
%     F_x = −½ρS·C_D0·V²          C_D0 = 0,05 (projeto)
%     F_y = +½ρS·C_Yβ·V·v         C_Yβ = −0,196 (doublets de aileron)
%     F_z = −½ρS·(C_L0·V² + C_Lα·V·w)   C_L0 = 0,34, C_Lα = 2,94 (doublets de profundor)
%   V, v, w da velocidade estimada (appdata 'aero_vel_now' = [V v w u]).
    av = getappdata(0,'aero_vel_now');
    fa = [0;0;0];
    if isempty(av) || av(1) < 0.05, return; end
    V = av(1); v = av(2); w = av(3);
    kA = aero_gains();  qS = 0.5*kA.rho*kA.S;
    CD0 = 0.05; CYb = -0.196; CL0 = 0.34; CLa = 2.94;
    fa = [ -qS*CD0*V^2 ;  qS*CYb*V*v ;  -qS*(CL0*V^2 + CLa*V*w) ];
end

function fa = faero_flatplate_local()
%FAERO_FLATPLATE_LOCAL  Força aerodinâmica da estrutura em pairado (placa plana).
%   Lombaerts et al. 2020, Eqs. 9 e 7:  [C_D;C_Y;C_L] = C_p·[senα·cosβ; senβ; senα·cosα],
%   levados ao corpo por R1(α)·R2(β) aplicados a [−C_D; −C_Y; −C_L], vezes ½ρV²S.
%   α, β, V da velocidade estimada publicada em appdata 'aero_vel_now' = [V v w u].
    av = getappdata(0,'aero_vel_now');
    fa = [0;0;0];
    if isempty(av) || av(1) < 0.05, return; end
    V = av(1); v = av(2); w = av(3); u = av(4);
    al = atan2(w, max(u, 1e-6));                 % α = atan(w/u) (u positivo para frente)
    be = asin(max(min(v/V, 1), -1));
    kA = aero_gains();  Cp = 2;
    CD = Cp*sin(al)*cos(be);  CY = Cp*sin(be);  CL = Cp*sin(al)*cos(al);
    R1 = [cos(al) 0 -sin(al); 0 1 0; sin(al) 0 cos(al)];
    R2 = [cos(be) -sin(be) 0; sin(be) cos(be) 0; 0 0 1];
    C = R1*R2*[-CD; -CY; -CL];
    fa = 0.5*kA.rho*V^2*kA.S * C;
end
