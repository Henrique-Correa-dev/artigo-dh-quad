function params = setup_params()
% SETUP_PARAMS  Configura todos os parametros do modelo VTOL
%
% Reune em uma unica struct:
%   - Propriedades de massa do CAD (SolidWorks)
%   - Geometria do frame H
%   - Mapeamento dos motores (C1..C4)
%   - Modelos de bancada (polinomios PWM->N, PWM->Nm)
%   - Parametros aerodinamicos (iniciais, para identificacao)
%
% Uso:
%   params = setup_params();
%
% =========================================================================

    %% ====================================================================
    %  1. PROPRIEDADES DE MASSA (CAD - SolidWorks)
    % =====================================================================
    params.mass = 1.07201;     % massa total [kg]
    params.g    = 9.80665;     % aceleracao gravitacional [m/s^2]

    % Tensor de inercia COMPLETO no CG [kg.m^2]
    % Fonte: SolidWorks, notacao de tensor positivo
    % Unidade original: gramas * metros^2 -> dividir por 1000
    %
    %   Ixx = 43.244    Ixy = -1.003    Ixz = 1.571
    %   Iyx = -1.003    Iyy = 84.404    Iyz = 0.021
    %   Izx = 1.571     Izy = 0.021     Izz = 126.192
    %
    Jx  =  43.244e-3;   % [kg.m^2]
    Jy  =  84.404e-3;
    Jz  = 126.192e-3;
    Jxy =  -1.003e-3;   % produto de inercia (com sinal do tensor)
    Jxz =   1.571e-3;
    Jyz =   0.021e-3;

    params.J = [Jx,  Jxy, Jxz;
                Jxy, Jy,  Jyz;
                Jxz, Jyz, Jz ];

    params.J_inv = inv(params.J);

    % Valores individuais (para referencia/debug)
    params.Jx  = Jx;
    params.Jy  = Jy;
    params.Jz  = Jz;
    params.Jxy = Jxy;
    params.Jxz = Jxz;
    params.Jyz = Jyz;

    %% ====================================================================
    %  2. GEOMETRIA DO FRAME H
    % =====================================================================
    %
    %              NARIZ (frente, +x)
    %                  |
    %     C3(CCW) o----+----o C1(CW)       <- lx_f = 311.18 mm
    %              |   CG    |
    %     C2(CW)  o---------o C4(CCW)      <- lx_r = 342.87 mm
    %           ly              ly
    %         232mm            232mm
    %
    lx_f = 0.31118;   % braco frontal [m]
    lx_r = 0.34287;   % braco traseiro [m]
    ly   = 0.232;     % braco lateral [m]

    %  Canal   Posicao              x        y       sentido
    %  C1      Frontal-Direito     +lx_f    +ly      CW (+1)
    %  C2      Traseiro-Esquerdo   -lx_r    -ly      CW (+1)
    %  C3      Frontal-Esquerdo    +lx_f    -ly      CCW(-1)
    %  C4      Traseiro-Direito    -lx_r    +ly      CCW(-1)

    params.x_m = [+lx_f; -lx_r; +lx_f; -lx_r];
    params.y_m = [+ly;   -ly;   -ly;   +ly  ];
    params.d_m = [+1;    +1;    -1;    -1   ];

    params.lx_f = lx_f;
    params.lx_r = lx_r;
    params.ly   = ly;

    %% ====================================================================
    %  3. MODELOS DE BANCADA (polinomios grau 3)
    % =====================================================================
    pwm_bench = [1000; 1200; 1400; 1600; 1800; 2000];
    thrust_g  = [   0;  143;  328;  532;  784;  843];
    torque_Nm = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];

    thrust_N = thrust_g * params.g / 1000;

    params.coeffs_T = polyfit(pwm_bench, thrust_N, 3);
    params.coeffs_Q = polyfit(pwm_bench, torque_Nm, 3);

    % Dead zone
    idx_T0 = find(thrust_N > 1e-9, 1, 'first');
    idx_Q0 = find(torque_Nm > 1e-9, 1, 'first');
    pwm_min_T = pwm_bench(idx_T0);
    pwm_min_Q = pwm_bench(idx_Q0);

    cT = params.coeffs_T;
    cQ = params.coeffs_Q;
    params.T_ref = @(pwm) (pwm >= pwm_min_T) .* max(0, polyval(cT, pwm));
    params.Q_ref = @(pwm) (pwm >= pwm_min_Q) .* max(0, polyval(cQ, pwm));

    %% ====================================================================
    %  4. FATORES DE ESCALA DOS MOTORES (chute inicial = 1.0)
    % =====================================================================
    params.k_T = [1.0; 1.0; 1.0; 1.0];   % escala empuxo por motor
    params.k_Q = [1.0; 1.0; 1.0; 1.0];   % escala torque por motor

    %% ====================================================================
    %  5. AMORTECIMENTO E BIAS (chute inicial = 0)
    % =====================================================================
    % Rotacional
    params.Dp = 0;    % amortecimento roll [1/s]
    params.Dq = 0;    % amortecimento pitch [1/s]
    params.Dr = 0;    % amortecimento yaw [1/s]

    params.Bp = 0;    % bias roll [rad/s^2]
    params.Bq = 0;    % bias pitch [rad/s^2]
    params.Br = 0;    % bias yaw [rad/s^2]

    % Translacional
    params.Xu = 0;    % arrasto em x [1/s]
    params.Yv = 0;    % arrasto em y [1/s]
    params.Zw = 0;    % arrasto em z [1/s]
    params.Bz = 0;    % bias vertical [m/s^2]

end
