function varargout = motor_models(varargin)
%MOTOR_MODELS  Modelos PWM→Empuxo (T) e PWM→Contra-torque (Q) via spline Akima.
%
% Une as funções antes separadas em create_thrust_model.m e create_torque_model.m.
%
% USO:
%   [func_T, func_Q] = motor_models()
%       Retorna handles construídos com a TABELA PADRÃO de bancada
%       (1 motor de referência — substitua quando tiver os 4).
%
%   func_T = motor_models('thrust', pwm_pts, T_pts)
%       Constrói só o modelo de empuxo a partir de pontos custom.
%       T_pts em GRAMAS (compat com chamadas antigas que passavam thrust_grams_exp).
%
%   func_Q = motor_models('torque', pwm_pts, Q_pts)
%       Constrói só o modelo de torque. Q_pts em N·m.
%
% Cada handle retornado:
%   - Satura PWM em [1000, 2000] µs
%   - Garante saída >= 0 (max(0, ·))
%   - Interpolação Akima ('makima') — numericamente idêntico ao Akima do Simulink

    if nargin == 0
        % Modo default: usa tabela de bancada padrão
        [pwm_bp, T_grams, Q_Nm] = bench_table();
        varargout{1} = build_thrust_handle(pwm_bp, T_grams);
        varargout{2} = build_torque_handle(pwm_bp, Q_Nm);
        return;
    end

    mode = varargin{1};
    pwm_bp = varargin{2}(:);

    switch lower(mode)
        case 'thrust'
            T_grams = varargin{3}(:);
            varargout{1} = build_thrust_handle(pwm_bp, T_grams);
        case 'torque'
            Q_Nm = varargin{4-1}(:);   % 3o arg (PWM é 2o)
            varargout{1} = build_torque_handle(pwm_bp, Q_Nm);
        otherwise
            error('motor_models:badMode', ...
                'Modo desconhecido: "%s". Use "thrust" ou "torque".', mode);
    end
end


function func = build_thrust_handle(pwm_bp, thrust_grams)
% Spline cúbica Akima. Saída em Newtons.
    pwm_min = 1000; pwm_max = 2000;
    T_N = thrust_grams * 9.80665 / 1000;   % g → N
    func = @(pwm) max(0, interp1(pwm_bp, T_N, ...
        min(max(pwm, pwm_min), pwm_max), 'makima', 'extrap'));

    fprintf('Modelo de Empuxo: Akima sobre %d pontos | PWM [%d,%d] | T>=0\n', ...
        numel(pwm_bp), pwm_min, pwm_max);
end


function func = build_torque_handle(pwm_bp, torque_Nm)
    pwm_min = 1000; pwm_max = 2000;
    func = @(pwm) max(0, interp1(pwm_bp, torque_Nm, ...
        min(max(pwm, pwm_min), pwm_max), 'makima', 'extrap'));

    fprintf('Modelo de Torque: Akima sobre %d pontos | PWM [%d,%d] | Q>=0\n', ...
        numel(pwm_bp), pwm_min, pwm_max);
end


function [pwm_bp, T_grams, Q_Nm] = bench_table()
% Tabela de bancada do motor de referência (RC Benchmark / Tyto Robotics).
% Substituir por dados por motor quando os 4 forem medidos.
    pwm_bp  = [1000; 1200; 1400; 1600; 1800; 2000];
    T_grams = [0;    143;  328;  532;  784;  843];
    Q_Nm    = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176];
end
