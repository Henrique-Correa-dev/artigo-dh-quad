function T = thrustFromPWM(pwm)
% thrustFromPWM Retorna empuxo [N] do Motor Azul para PWM 1000–2000
% pwm: vetor de PWM (ex: [PWM1 PWM2 PWM3 PWM4])
% T: vetor de empuxos em Newtons

% Dados experimentais do motor azul
pwm_vals = [1000, 1200, 1400, 1600, 1800, 2000];
thrust_g = [0, 143, 328, 532, 784, 843];  % em gramas

% Conversão para Newtons
thrust_N = thrust_g * 9.81 / 1000;

% Interpolação
T = interp1(pwm_vals, thrust_N, pwm, 'linear', 'extrap');
end
