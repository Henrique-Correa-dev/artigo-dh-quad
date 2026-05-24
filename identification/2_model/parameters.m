function p = parameters()
%PARAMETERS  Constantes físicas e chute inicial P0 do drone DH2 (VTOL híbrido).
%
% USO:
%   p = parameters();
%   p.m             % massa total (kg)
%   p.g             % aceleração da gravidade (m/s²)
%   p.tau_motor     % constante de tempo do lag do motor (s)
%   p.P0_J          % chute inicial P_J (24×1)
%   p.bounds.lb     % lower bounds (24×1)
%   p.bounds.ub     % upper bounds (24×1)
%   p.param_names   % nomes (1×24 cell)
%   p.bench         % struct com tabela de bancada
%
% Fonte oficial dos números: slide "Identificacao do modelo aerodinamico"
% (docs/arquitetura_pa.pdf, tabela 1.0).

    %% Constantes físicas
    p.m           = 2.20;     % massa total medida (slide oficial)
    p.g           = 9.81;
    p.tau_motor   = 0.05;     % constante de tempo BLDC+ESC (s)

    %% Inércias (do slide XFLR5, convenção FRD)
    % Magnitudes do slide: Ix=0.144, Iy=0.116, Iz=0.257 kg·m²
    % NOTA: código atual usa valores de CAD com massa errada — VERIFICAR.
    %       Tarefa #73 ainda pendente: reconciliar.
    p.J.Jx  = 63.244 / 1000;   % CAD legacy — usar este até confirmar mapeamento
    p.J.Jy  = 250.554 / 1000;
    p.J.Jz  = 116.192 / 1000;
    p.J.Jxz = 1.571 / 1000;

    %% Geometria (braços dos rotores até CG — confirmado correto)
    p.arms.Lx_base = 0.232;       % lateral (mesma esq=dir)
    p.arms.Ly_f    = 0.311185;    % CG-frente
    p.arms.Ly_r    = 0.342865;    % CG-trás

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

    %% Chute inicial P_J (24×1) para identificação
    p.P0_J = [p.J.Jx; p.J.Jy; p.J.Jz; p.J.Jxz; ...
              0.55; 0.45; 1.0; 0.75;       % k_T1..k_T4 (chute legado)
              0.55; 0.45; 1.0; 0.75;       % k_Q1..k_Q4
              10; 5; 0.5;                  % Dp, Dq, Dr
              0.7; 1.4; 0.3;               % Bp, Bq, Br
              0; 0;                        % dx_cg, dy_cg
              -4; -4; -0.1; -0.5];         % Xu, Yv, Zw, Bz

    %% Bounds (lb/ub) — atualmente largos demais, ver task #66
    p.bounds.lb = [0.032; 0.125; 0.058; 0.0001; ...
                   0.05; 0.05; 0.05; 0.05; ...
                   0.10; 0.10; 0.10; 0.10; ...
                   0; 0; 0; ...
                   -10; -10; -10; ...
                   -0.08; -0.05; ...
                   -30; -30; -2; -5];

    p.bounds.ub = [0.095; 0.376; 0.174; 0.006; ...
                   5; 5; 5; 5; ...
                   3.0; 3.0; 3.0; 3.0; ...
                   20; 10; 10; ...
                   10; 10; 10; ...
                   0.08; 0.05; ...
                   0; 0; 0; 5];

    p.param_names = {'Jx','Jy','Jz','Jxz', ...
        'k_T1','k_T2','k_T3','k_T4','k_Q1','k_Q2','k_Q3','k_Q4', ...
        'Dp','Dq','Dr','Bp','Bq','Br', ...
        'dx_cg','dy_cg', ...
        'Xu_m','Yv_m','Zw_m','Bz'};
end
