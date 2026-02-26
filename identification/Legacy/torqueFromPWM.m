function Q = torqueFromPWM(pwm)
% torqueFromPWM Retorna torque [Nm] do Motor Azul para PWM 1000–2000
% pwm: vetor de PWM (ex: [PWM1 PWM2 PWM3 PWM4])
% Q: vetor de torques em Newton-metro

% Dados experimentais do motor azul
pwm_vals = [1000, 1200, 1400, 1600, 1800, 2000];
torque_vals = [0.000, 0.034, 0.070, 0.115, 0.171, 0.176];

% Interpolação
Q = interp1(pwm_vals, torque_vals, pwm, 'linear', 'extrap');
end
