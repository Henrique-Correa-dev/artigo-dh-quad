function p = parameters(varargin)
%PARAMETERS  Constantes físicas e chute inicial P0 do drone DH2 (VTOL híbrido).
%
% USO:
%   p = parameters();
%   p = parameters('campaign', 'fev2024');   % fixa a campanha ATIVA do processo
%   p = parameters('campaign', '');          % volta ao default (mai/2026)
%
% CAMPANHAS: cada voo foi feito com uma configuração diferente da aeronave
% (motores, carga paga, posição do CG). Os valores default deste arquivo são os
% da campanha de MAIO DE 2026, que é a da dissertação. Outras campanhas ficam em
% 2_model/campaigns/<nome>.m, cada uma um arquivo que recebe p e devolve p com
% os campos que mudam. A campanha escolhida vale para TODO o processo MATLAB
% (guardada em variável persistente), de modo que as funções internas que
% chamam parameters() sem argumento também a enxergam.
%
%   p.m             % massa total (kg)
%   p.g             % aceleração da gravidade (m/s²)
%   p.P0_J          % chute inicial P_J (15×1)
%   p.bounds.lb     % lower bounds (15×1)
%   p.bounds.ub     % upper bounds (15×1)
%   p.param_names   % nomes (1×15 cell)
%   p.bench         % struct com tabela de bancada
%
% Fonte oficial dos números: slide "Identificacao do modelo aerodinamico"
% (docs/arquitetura_pa.pdf, tabela 1.0).

    persistent ACTIVE_CAMPAIGN
    if nargin >= 1 && ischar(varargin{1}) && strcmpi(varargin{1}, 'campaign')
        if nargin >= 2, ACTIVE_CAMPAIGN = varargin{2}; else, ACTIVE_CAMPAIGN = ''; end
    end
    camp_ = ACTIVE_CAMPAIGN;  if isempty(camp_), camp_ = ''; end

    %% Constantes físicas
    p.m           = 1.993;     % massa total medida (slide oficial)
    p.g           = 9.81;
    p.rho         = 1.225;     % densidade do ar [kg/m³] — nível do mar, ISA.
    %  Todos os coeficientes aerodinâmicos (AVL, XFLR5, prior_damping) foram
    %  calculados com este valor. O ginásio fica a ~600 m (ρ ≈ 1,15), diferença
    %  de 6% que é absorvida pelos próprios coeficientes identificados.

    %% Inércias (CAD SolidWorks — centro de massa, frame de saída, FRD)
    %  Valores oficiais do CAD em g·m² (→ /1000 = kg·m²). Satisfazem as três
    %  desigualdades triangulares (Jz ≤ Jx+Jy com folga de só ~1.2% → drone
    %  quase planar, esperado). O P0 legacy anterior (Jy=250.5!) estava errado.
    %  Produtos Ixy=-1.0, Iyz=0.02 g·m² desprezíveis; só Ixz=1.571 é relevante.
    p.J.Jx  =  43.244 / 1000;   % Ixx
    p.J.Jy  =  84.404 / 1000;   % Iyy
    p.J.Jz  = 126.192 / 1000;   % Izz
    p.J.Jxz =   1.571 / 1000;   % Ixz (único produto relevante)

    %% Geometria (braços dos rotores até CG — fonte única de verdade)
    %  Lx_r e Lx_l são separados pra permitir CG deslocado lateralmente.
    %  Default: simétrico (Lx_r == Lx_l = 0.232 m). Edite se o CG real estiver
    %  lateralmente offset (motores 1,4 do lado direito vs 2,3 do esquerdo).
    p.arms.Lx_r    = 0.232-0.0045;       % direita (motores 1, 4) — CG até motor lado dir
    p.arms.Lx_l    = 0.232+0.0045;       % esquerda (motores 2, 3) — CG até motor lado esq
    p.arms.Ly_f    = 0.323+0.0082;    % CG-frente
    p.arms.Ly_r    = 0.330-0.0082;    % CG-trás

    %% Posição do IMU em relação ao CG (body frame FRD: x frente, y direita, z baixo)
    %  Vetor do CG até o sensor. Causa acoplamento ω_dot ↔ acc linear:
    %     a_imu = a_cg + α × r_imu + ω × (ω × r_imu)
    %  Convenção do drone atual:
    %     rx > 0 → IMU à frente do CG     (drone tem IMU ~2cm à frente)
    %     rz > 0 → IMU abaixo do CG       (z aponta pra baixo)
    p.imu_offset = [+0.10; 0.00; +0.02];  % [rx; ry; rz] em metros
    rxo_ = getappdata(0,'imu_rx_override');        % teste: r_x livre/alterado
    if ~isempty(rxo_) && isnumeric(rxo_), p.imu_offset(1) = rxo_; end

    %% Queda de empuxo com a velocidade axial (influxo): T ← T − 4·k_v·w
    %  k_v MEDIDO no voo (2_model/measure_kv.m): 0,226 ± 0,013 N/(m/s). É o mesmo
    %  mecanismo que dá o termo de influxo do L_p; aqui entra no EMPUXO TOTAL,
    %  no canal vertical, com w da velocidade estimada de bordo (exógena).
    %  Default 0 (desligado); ligado por appdata 'kv_thrust' = k_v.
    p.k_v_thrust = 0;
    kvo_ = getappdata(0,'kv_thrust');
    if ~isempty(kvo_) && isnumeric(kvo_), p.k_v_thrust = kvo_; end

    %% Desalinhamento angular entre o plano dos rotores (eixo do empuxo) e o sensor
    %  Achado do 5_validation/drag_probe.m no trecho com GPS: resíduo de a_x
    %  ∝ −0,0314·(T/m) (R² 0,46) → eps_y = asin(0,0314) ≈ 1,8°. Em a_y não há
    %  correlação com T (eps_x ≈ 0). No pairado T/m ≈ g → −0,31 m/s², que era o
    %  "bias" de x hardcoded no accelerometer_model; agora o bias de x é ~0.
    p.imu_tilt = [0.0; deg2rad(1.8)];     % [eps_x; eps_y] em rad
    p.imu_bias = [0.0; -0.2; +0.4];       % [bx; by; bz] em m/s² (era [-0.3 -0.2 +0.4])

    %% Arrasto translacional de rotor (modelo do livro, Eq. 6.31): f = -k_drag·v_body
    %  Atua em x,y do corpo (não em z). Entra na EOM (amortece u,v → tira o drift da
    %  integração) E na força específica do acelerômetro (accX/accY).
    %  TESTE: valor FIXO aqui. k_drag=0 → desliga (modelo antigo). Depois vira P.
    p.k_drag = 0;     % N/(m/s) — ajustar no teste; 0 desliga

    %% Asa (do slide, ainda não usado mas reservado pra task #67)
    p.wing.S       = 0.27;     % área (m²)
    p.wing.b       = 1.20;     % envergadura (m)
    p.wing.c       = 0.226;    % corda média (m)
    p.wing.airfoil = 'USA-35B';
    p.wing.V_ref   = 15;       % velocidade de cruzeiro (m/s)
    p.wing.X_CG    = -0.230;   % CG longitudinal vs bordo de ataque (m)
    p.wing.Z_CG    = -0.05;    % CG vertical vs plano dos rotores (m)
    %  Termo de trim de arfagem (arqueamento do perfil), identificado nos
    %  doublets de profundor (11_fixed_wing, fw_longitudinal_3win). Mantido
    %  fixo no modelo para equivalência estrutural com os executáveis de asa
    %  fixa (Sato usa 0,109; Ana 0,083). Momento: +½ρSc̄·C_m0·V².
    p.wing.C_m0    = 0.012;
    %  Derivadas cruzadas medidas nos doublets de aileron (fw_lateral), também
    %  fixas, pela mesma equivalência estrutural (Sato e Ana as incluem):
    %    Mx += ¼ρSb²·C_lr·V·r      Mz += ¼ρSb²·C_np·V·p
    p.wing.C_lr    =  0.1420;
    p.wing.C_np    = -0.0319;

    %% Tabela de bancada — motor REAL (SunnySky Angel A2212-15 800 kV), média de 3 testes.
    %  Throttle %→PWM linear [1000,2000] (0%→1000, 100%→2000).
    %  Substitui o motor de REFERÊNCIA antigo (843 g máx), que superestimava o
    %  empuxo em ~30% no PWM de hover e causava o gap de AccZ.
    %  RPM medido → guardado p/ futuro modelo físico T ∝ RPM² (logs não têm ESC RPM).
    %  Torque (Q) NÃO foi medido neste teste → re-amostrado do motor de referência
    %  (mesma FORMA) e marcado A RE-MEDIR; o k_Q da identificação absorve a escala.
    p.bench.pwm     = [1000; 1167; 1333; 1500; 1667; 1833; 2000];
    p.bench.T_grams = [0;    90;   185;  285;  449;  600;  692];
    p.bench.RPM     = [0;    2432; 3378; 4102; 5087; 5791; 6063];
    p.bench.Q_Nm    = [0;    0.028; 0.058; 0.093; 0.134; 0.172; 0.176];  % ref motor (re-grid) — RE-MEDIR

    %% ============ CAMPANHA (sobrescreve o que muda de voo para voo) =======
    %  Aplicada AQUI: depois da geometria e da bancada, antes das constantes
    %  derivadas do motor e do P0, para que tudo abaixo use os valores da campanha.
    p.campaign = 'mai2026';
    if ~isempty(camp_)
        p.campaign = camp_;
        fh = str2func(camp_);          % 2_model/campaigns/<camp_>.m
        p  = fh(p);
    end

    %% ============ Modelo de motor FÍSICO (PWM→RPM→empuxo) ============
    %  Estrutura em 2 estágios (Quan 2017 / rotorcraft), usando o dado de RPM:
    %    [1] estático linear:  RPM_ss = CR·PWM + Omega_b   (ajuste de bancada, s/ idle)
    %    [2] dinâmica:         d(RPM)/dt = (RPM_ss − RPM)/tau_m   (lag 1ª ordem)
    %    [3] estático quadr.:  T = kT_rpm·RPM² ,  Q = kQ_rpm·RPM²
    %  k_T/k_Q identificados (P) entram como ESCALA POR MOTOR sobre o estágio [3].
    %  enable=false → usa o modelo antigo (PWM→empuxo via spline). Só o caminho de
    %  validação (sim_window 'full') checa este flag; NÃO afeta identify_plant/linearize.
    bb_ = p.bench.pwm(2:end);  rr_ = p.bench.RPM(2:end);     % s/ idle (motor parado)
    c_rpm_ = polyfit(bb_, rr_, 1);                            % RPM_ss = CR·PWM + Omega_b
    TN_b_  = p.bench.T_grams * 9.80665/1000;                  % empuxo bancada [N]
    p.motor.CR      = c_rpm_(1);                              % RPM/µs
    p.motor.Omega_b = c_rpm_(2);                              % RPM (intercepto da reta)
    p.motor.kT_rpm  = sum(TN_b_.*p.bench.RPM.^2)/sum(p.bench.RPM.^4);   % N/RPM²
    p.motor.kQ_rpm  = sum(p.bench.Q_Nm.*p.bench.RPM.^2)/sum(p.bench.RPM.^4); % N·m/RPM²
    p.motor.tau_m   = 0.05;     % [s] lag 1ª ordem da RPM (identificar do voo)
    p.motor.delay_s = 0.10;     % [s] atraso puro de transporte PWM→empuxo (ESC+log; = 1 amostra a 10 Hz)
    %  Cadeia de atuação ÚNICA: 2_model/motor_chain.m (atraso + RPM_ss + lag exato + T=kT·RPM²).
    %  identify_plant, sim_window, results_id_figures e compare_* usam TODOS essa função,
    %  para que identificação e validação tenham exatamente o mesmo motor.
    p.motor.enable  = true;    % (legado) modelo físico; a cadeia agora vem de motor_chain

    %% ============ FORMA DO AMORTECIMENTO ANGULAR ============
    %  'moment' (padrão)  P(13:15) = L_p, M_q, N_r em N·m·s, somados DENTRO do
    %                     vetor de momentos:  Mx ← Mx − L_p·p  etc. É a forma
    %                     coerente com a equação de corpo rígido (Beard 3.16),
    %                     em que toda a física entra por M e é multiplicada pela
    %                     inversa do tensor de inércia. Traz de brinde os dois
    %                     termos CRUZADOS via J_xz que a forma antiga não tinha:
    %                        ṗ recebe −Γ4·N_r·r     ṙ recebe −Γ4·L_p·p
    %  'rate' (legado)    P(13:15) = c_p, c_q, c_r em 1/s, subtraídos direto da
    %                     derivada do estado. Equivalência no próprio eixo:
    %                        c_p = Γ3·L_p ,  c_q = M_q/J_y ,  c_r = Γ8·N_r
    %  Para converter vetores antigos use 2_model/damp_convert.m.
    p.damp_form = 'moment';
    %  Override em tempo de execução, para scripts que comparam as duas formas:
    %     setappdata(0,'damp_form_override','rate')  /  rmappdata(0,'damp_form_override')
    dfo_ = getappdata(0,'damp_form_override');
    if ~isempty(dfo_) && ischar(dfo_), p.damp_form = dfo_; end

    %% Chute inicial P_J (15×1) para identificação
    %  Bp, Bq, Br REMOVIDOS — offset de CG capturado via Lx/Ly assimetria.
    %  Xu, Yv, Zw REMOVIDOS — drag translacional não-identificável sem
    %     ground truth de u, v, w (precisa GPS velocity).
    %  Bz REMOVIDO — vai pra modelo de sensor (bias do acelerômetro Z) depois.
    %  Modelo P agora é PURAMENTE ROTACIONAL + parâmetros de motor.
    %  cp, cq, cr: A PRIORI FÍSICO pela hipótese de influxo do rotor
    %  (2_model/prior_damping.m): cp0 = 4·k_v·ly²/Jx, cq0 = 4·k_v·lx²/Jy,
    %  cr0 = 4·k_h·d²/Jz com k_v = 0,239 e k_h = 0,016 N/(m/s) (hélice 1045,
    %  BEMT/Padfield). Substitui o chute arbitrário 0,5. O OEM parte do EEM,
    %  então este valor só afeta a coluna Θ0 (validação com chute inicial).
    %  MODELO ESTENDIDO (22 parâmetros). Os 15 primeiros são o modelo original
    %  (multirrotor puro); os 7 últimos são a AERODINÂMICA DA ESTRUTURA, na forma
    %  REGULAR (multiplicada por V, sem dividir), para não singularizar no pairado:
    %     P(16) C_d   arrasto induzido do rotor (Beard 14.4.2): f = −T·C_d·[u;v;0]
    %     P(17) Cl_p  L += ¼·ρ·S·b²·Cl_p·V·p        P(18) Cl_β  L += ½·ρ·S·b·Cl_β·V·v
    %     P(19) Cm_q  M += ¼·ρ·S·c̄²·Cm_q·V·q       P(20) Cm_α  M += ½·ρ·S·c̄·Cm_α·V·w
    %     P(21) Cn_r  N += ¼·ρ·S·b²·Cn_r·V·r        P(22) Cn_β  N += ½·ρ·S·b·Cn_β·V·v
    %  Ganhos dimensionais em 2_model/aero_gains.m (fonte única).
    %  Chute inicial dos aerodinâmicos: AVL do DH (2_model/avl/dh_st.txt), que
    %  concorda com o XFLR5 do relatório FINEP. C_d do 5_validation/drag_probe.m.
    %  ATENÇÃO: L_p,M_q,N_r (P13:15) e Cl_p,Cm_q,Cn_r (P17,19,21) são quase colineares
    %  no envelope voado (V ≈ 1 m/s); espera-se CR% alto nos aerodinâmicos.
    %  L_p, M_q, N_r: a priori CALCULADO DIRETO EM N·m·s por 2_model/
    %  prior_damping_moment.m, sem passar por inércia nenhuma (é a vantagem da
    %  forma de momento: a derivação física produz momento, e a divisão por J
    %  que gerava c_p é justamente o passo que injetava a incerteza do CAD).
    %  SÓ O INFLUXO, o único mecanismo com o parâmetro MEDIDO no próprio voo
    %  (k_v = 0,2385 N/(m/s) TEORICO, elemento de pá + quantidade de movimento,
    %  Apêndice A.1 da dissertação; o medido em voo era 0,226, measure_kv.m). A aerodinâmica da
    %  asa é termo separado (P(17:22)·V_a). O momento de cubo e a esteira foram
    %  retirados do a priori: o cubo depende de hipóteses sobre a pá (corda,
    %  inclinação de sustentação) sem medida, e a esteira não tem referência
    %  que a sustente como momento. Ficam no damping_budget.m como estimativa.
    %                            L_p       M_q       N_r
    %    influxo 4·k_v·l²      0,0513    0,1017    0,0000   (k_v TEORICO 0,2385, prior_damping_moment)
    %    força H               0,0000    0,0000    0,0100   4·k_h·(l_x²+l_y²), elemento de pá
    %  É o termo ∂Z/∂x·Δy da Eq. 55 de Nguyen & Webb (2025): variação do empuxo
    %  com a velocidade axial induzida pela rotação do corpo, vezes o braço.
    p.P0_J = [p.J.Jx; p.J.Jy; p.J.Jz; p.J.Jxz; ...
              1; 1; 1; 1;       % k_T1..k_T4
              1; 1; 1; 1;       % k_Q1..k_Q4
              0.0513; 0.1017; 0.0100; ...  % L_p, M_q, N_r [N·m·s] (ver nota acima)
              0.000; ...             % C_d = 0 — a aerodinâmica é SÓ ROTACIONAL
              -0.0769; -0.0222; ...  % Cl_p, Cl_β  MEDIDOS em asa fixa
              -5.540;  -0.435; ...   % Cm_q, Cm_α  MEDIDOS
              -0.2037;  0.0472; ...  % Cn_r, Cn_β  MEDIDOS
              0; 0; 0];              % L_v, M_w, N_v [N·s] — rotor, momento por velocidade de
    %  translação (NASA: C_lv, C_mw, C_nv). Análogo de rotor dos estáticos da asa,
    %  sem o fator V. Default ZERO; só entram em teste com v, w estimados.
    %  TODOS OS SEIS MEDIDOS no voo de asa fixa de 17/12/2024 (fw_longitudinal.m,
    %  fw_lateral.m). Longitudinais: média dos 4 doublets de profundor. Látero-
    %  direcionais: média dos doublets 1, 2 e 4 de aileron; o doublet 3
    %  (322–326 s) é excluído para Cl_β e Cn_β porque não tem informação de β
    %  (sai zero com CR de 298% e 8·10¹⁰%), e contaminava a média. Com ele fora:
    %     Cl_β −0,022 ± 0,008    Cn_β +0,047 ± 0,006   (desvio 13% e 34%)
    %  Os de amortecimento mudam < 6% com ou sem o doublet 3 e ficam na média dos 4.
    %  Para referência, AVL: Cl_β −0,062, Cn_β +0,066, Cm_α −0,740 (todos maiores).

    %% Bounds (lb/ub)
    %
    %  INÉRCIAS: bounds APERTADOS em ±10% do CAD. Yaw é o canal fraco do voo →
    %  Jz é quase inobservável; sem isto o optimizer inflava Jz (chegou a 0.235,
    %  violando Jz ≤ Jx+Jy). Como a folga triangular é só ~1.2%, ±10% independente
    %  ainda PODE estourar a desigualdade → o identify_plant adiciona uma
    %  PENALIDADE TRIANGULAR (REG_LAMBDA.tri) pra garantir consistência física.
    %
    %  k_T/k_Q: ±40% em torno de 1.0 — bound largo porque diagnose_mz mostrou
    %  que bancada subestima Q em voo por ~27% (provavelmente também subestima
    %  T por similar ordem). Cada k é INDIVIDUAL por motor.

    % CAD (kg·m²): Jx=0.04324, Jy=0.08440, Jz=0.12619, Jxz=0.00157  (15 elementos)
    %  Jx, Jy: LIVRES em ±25% do CAD — são OBSERVÁVEIS no voo (roll/pitch bem
    %    excitados, CRB ~4%). O dado quer subir um pouco vs o CAD; deixa subir.
    %  Jz, Jxz: TRAVADOS no CAD (lb=ub). Liberá-los foi testado e o optimizer
    %    levou Jz a 0.141 > Jx+Jy (VIOLA o triângulo) e Jxz ao limite (CRB 80%),
    %    sem ganho de ajuste → travar no CAD é mais honesto e mantém o modelo físico.
    p.bounds.lb = [0.03243; 0.06330; 0.126192; 0.001571; ...  % Jx,Jy ±25%; Jz,Jxz=CAD (TRAVADOS)
                   0.40; 0.40; 0.40; 0.40; ...        % k_T
                   -1.40; -1.40; -1.40; -1.40; ...    % k_Q
                   0; 0; 0; ...                       % L_p, M_q, N_r >= 0
                   0.0; ...                           % C_d >= 0
                   -4.0; -0.6; -90.0; -7.4; -0.7; 0.0; ...  % aero momentos: 10× o AVL
                   -1.0; -1.0; -1.0];                    % L_v, M_w, N_v [N·s]

    p.bounds.ub = [0.05406; 0.10551; 0.126192; 0.001571; ...  % Jz,Jxz travados=CAD
                   1.40; 1.40; 1.40; 1.40; ...
                   1.40; 1.40; 1.40; 1.40; ...
                   1.5; 1.5; 1.5; ...                 % L_p, M_q, N_r [N·m·s]
                   0.30; ...                       % C_d  (limite alto de Beard)
                   0.0; 0.0; 0.0; 0.0; 0.0; 0.5; ...  % aero: sinais físicos (Cn_β > 0)
                   1.0; 1.0; 1.0];                    % L_v, M_w, N_v, sinal livre
    %  LIBERAÇÃO OPCIONAL de Jz e Jxz (appdata 'jz_free'): abre os bounds para
    %  o dado poder movê-los, com a desigualdade triangular do custo como
    %  restrição física efetiva. Default: travados (comportamento oficial).
    jzf_ = getappdata(0,'jz_free');
    if ~isempty(jzf_) && isequal(jzf_,1)
        p.bounds.lb(3) = 0.75*0.126192;   p.bounds.ub(3) = 1.25*0.126192;
        p.bounds.lb(4) = 0.80*0.001571;   p.bounds.ub(4) = 2.00*0.001571;
    end

    %  Aerodinâmicos: bounds one-sided pelo SINAL FÍSICO (amortecimento < 0,
    %  Cl_β < 0 diedro, Cm_α < 0 estabilidade estática, Cn_β > 0 deriva),
    %  com folga de 10× o valor do AVL para o dado poder falar.

    %% ====== POR QUE A AERODINÂMICA É SÓ ROTACIONAL NESTE MODELO ======
    %  As forças aerodinâmicas da estrutura (X_u, Y_v, Z_w) e o arrasto induzido
    %  de rotor (C_d) foram TESTADOS e retirados. O motivo não é o termo, é o
    %  ESTADO que ele multiplica: em malha aberta os estados translacionais
    %  divergem, porque o erro das taxas integra na atitude (roll simulado chega
    %  a 126° contra 28° medidos em 20 s), a gravidade projetada fica errada em
    %  g·sen(37°) e a velocidade dispara. Medido na janela de validação:
    %       u simulado −142 m/s | v +82 | w −165      contra ±2,7 m/s reais
    %  Qualquer força proporcional a u, v, w multiplica esse número. Ligando as
    %  forças, o R² de a_z cai de +0,48 para −56 (outputs/runs/aeroF_fix_2026);
    %  ligando só o C_d = 0,010, cai para −125 (runs/aero_avl_2026).
    %  Os MOMENTOS aerodinâmicos não têm esse problema porque são escalados pela
    %  velocidade ESTIMADA de bordo (exógena) e multiplicam p, q, r, que são
    %  estados bem comportados. Por isso P(17:22) ficam e P(16) vai a zero.
    %  A maquinaria das forças continua no código (kdm_of.m, vetor kdm de 3
    %  colunas), inerte com P de 22 elementos, para quando houver realimentação
    %  ou medida de velocidade que impeça a divergência.
    %
    %  ANTIGO BLOCO DE FORÇAS (removido do vetor de parâmetros):
    %  Até aqui a estrutura só entrava nos MOMENTOS. As forças dela não existiam
    %  no modelo: as equações de u̇, v̇, ẇ tinham só gravidade, Coriolis, empuxo e
    %  o arrasto induzido de rotor (C_d). Estes três fecham essa lacuna, na mesma
    %  forma regular ∝ V que os momentos usam, para não singularizar no pairado:
    %
    %     f_x += −½ρS·C_Xu·V·u        f_y += −½ρS·C_Yv·V·v        f_z += −½ρS·C_Zw·V·w
    %
    %  São os clássicos X_u, Y_v e Z_w, que tinham sido removidos do modelo por
    %  não serem identificáveis sem medida de u, v, w. Voltam agora porque a
    %  velocidade V vem de estimate_velocity e os coeficientes têm a priori do
    %  voo de asa fixa:
    %     C_Zw = C_Lα  = 2,94   medido nos doublets de profundor (11_fixed_wing)
    %     C_Yv = |C_Yβ| = 0,196 medido nos doublets de aileron
    %     C_Xu ≈ 2·C_D0 = 0,05  sem medida própria, valor de projeto
    %  Em V = 1,29 m/s isso dá ½ρS·C_Zw·V/m = 0,31 1/s, contra g·C_d = 0,019 1/s
    %  do arrasto de rotor, ou seja, dezessete vezes maior no canal vertical.
    p.param_names = {'Jx','Jy','Jz','Jxz', ...
        'k_T1','k_T2','k_T3','k_T4','k_Q1','k_Q2','k_Q3','k_Q4', ...
        'L_p','M_q','N_r', ...
        'C_d','Cl_p','Cl_b','Cm_q','Cm_a','Cn_r','Cn_b', ...
        'L_v','M_w','N_v'};
    p.n_rot_only = 15;   % nº de parâmetros do modelo original (sem aerodinâmica)
end
