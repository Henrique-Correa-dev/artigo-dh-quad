function p = parameters()
%PARAMETERS  Constantes físicas e chute inicial P0 do drone DH2 (VTOL híbrido).
%
% USO:
%   p = parameters();
%   p.m             % massa total (kg)
%   p.g             % aceleração da gravidade (m/s²)
%   p.tau_motor     % constante de tempo do lag do motor (s)
%   p.P0_J          % chute inicial P_J (15×1)
%   p.bounds.lb     % lower bounds (15×1)
%   p.bounds.ub     % upper bounds (15×1)
%   p.param_names   % nomes (1×15 cell)
%   p.bench         % struct com tabela de bancada
%
% Fonte oficial dos números: slide "Identificacao do modelo aerodinamico"
% (docs/arquitetura_pa.pdf, tabela 1.0).

    %% Constantes físicas
    p.m           = 1.91;     % massa total medida (slide oficial)
    p.g           = 9.81;
    p.tau_motor   = 0;     % constante de tempo BLDC+ESC (s)

    %% Inércias (do slide XFLR5, convenção FRD)
    % Magnitudes do slide: Ix=0.144, Iy=0.116, Iz=0.257 kg·m²
    % NOTA: código atual usa valores de CAD com massa errada — VERIFICAR.
    %       Tarefa #73 ainda pendente: reconciliar.
    p.J.Jx  = 63.244 / 1000;   % CAD legacy — usar este até confirmar mapeamento
    p.J.Jy  = 250.554 / 1000;
    p.J.Jz  = 116.192 / 1000;
    p.J.Jxz = 1.571 / 1000;

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

    %% Asa (do slide, ainda não usado mas reservado pra task #67)
    p.wing.S       = 0.27;     % área (m²)
    p.wing.b       = 1.20;     % envergadura (m)
    p.wing.c       = 0.226;    % corda média (m)
    p.wing.airfoil = 'USA-35B';
    p.wing.V_ref   = 15;       % velocidade de cruzeiro (m/s)
    p.wing.X_CG    = -0.230;   % CG longitudinal vs bordo de ataque (m)
    p.wing.Z_CG    = -0.05;    % CG vertical vs plano dos rotores (m)

    %% Tabela de bancada (motor de referência)
    p.bench.pwm     = [1000; 1200; 1400; 1600; 1800; 2000];
    p.bench.T_grams = [0;    143;  328;  532;  784;  843];
    p.bench.Q_Nm    = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];

    %% Chute inicial P_J (15×1) para identificação
    %  Bp, Bq, Br REMOVIDOS — offset de CG capturado via Lx/Ly assimetria.
    %  Xu, Yv, Zw REMOVIDOS — drag translacional não-identificável sem
    %     ground truth de u, v, w (precisa GPS velocity).
    %  Bz REMOVIDO — vai pra modelo de sensor (bias do acelerômetro Z) depois.
    %  Modelo P agora é PURAMENTE ROTACIONAL + parâmetros de motor.
    p.P0_J = [p.J.Jx; p.J.Jy; p.J.Jz; p.J.Jxz; ...
              1; 1; 1; 1;       % k_T1..k_T4
              1; 1; 1; 1;       % k_Q1..k_Q4
              1; 1; 1];    % Dp, Dq, Dr

    %% Bounds (lb/ub)
    %
    %  INÉRCIAS: bounds APERTADOS em torno do CAD (±20%) — o optimizer estava
    %  saturando Jxz no UB (0.006) e criando acoplamento G4·Mx espúrio no r_dot,
    %  que era compensado por k_Q assimétrico. Jxz ≈ 0 é o esperado pra quad
    %  simétrico.
    %
    %  k_T/k_Q: ±40% em torno de 1.0 — bound largo porque diagnose_mz mostrou
    %  que bancada subestima Q em voo por ~27% (provavelmente também subestima
    %  T por similar ordem). Cada k é INDIVIDUAL por motor.

    % P0_J ref: Jx=0.0632, Jy=0.2506, Jz=0.1162, Jxz=0.00157  (15 elementos)
    p.bounds.lb = [0.040; 0.050; 0.090; 0.00; ...   % inércias ±20% do CAD, Jxz ~0
                   0.40; 0.40; 0.40; 0.40; ...        % k_T
                   -1.40; -1.40; -1.40; -1.40; ...    % k_Q
                   0; 0; 0];                          % Dp, Dq, Dr

    p.bounds.ub = [0.080; 0.310; 0.300; 0.01; ...   % Jxz limitado a 0.003
                   1.40; 1.40; 1.40; 1.40; ...
                   1.40; 1.40; 1.40; 1.40; ...
                   20; 10; 10];

    p.param_names = {'Jx','Jy','Jz','Jxz', ...
        'k_T1','k_T2','k_T3','k_T4','k_Q1','k_Q2','k_Q3','k_Q4', ...
        'Dp','Dq','Dr'};
end
