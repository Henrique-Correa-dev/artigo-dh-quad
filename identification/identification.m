load("log_data.mat")

t1 = 147;    % tempo inicial em segundos
t2 = 157;    % tempo final em segundos

ATT.TimeS = double(ATT.TimeUS) / 1e6;
IMU.TimeS = double(IMU.TimeUS) / 1e6;
RCOU.TimeS = double(RCOU.TimeUS) / 1e6;
GPS.TimeS = double(GPS.TimeUS) / 1e6;

idx = IMU.I == 0;
gyrX_raw = IMU.GyrX(idx);
gyrY_raw = IMU.GyrY(idx);
gyrZ_raw = IMU.GyrZ(idx);
accX_raw = IMU.AccX(idx);
accY_raw = IMU.AccY(idx);
accZ_raw = IMU.AccZ(idx);
time_IMU = IMU.TimeS(idx);

roll_raw = ATT.Roll;
pitch_raw = ATT.Pitch;
yaw_raw = ATT.Yaw;
time_ATT = ATT.TimeS;

lat_raw = GPS.Lat;
lon_raw = GPS.Lng;
alt_raw = GPS.Alt;
time_GPS = GPS.TimeS;

pwm1_raw = double(RCOU.C1);
pwm2_raw = double(RCOU.C2);
pwm3_raw = double(RCOU.C3);
pwm4_raw = double(RCOU.C4);
time_RCOU = RCOU.TimeS;

% 1. Definir tempo comum baseado na sobreposição dos tempos
t_start = max([min(time_IMU), min(time_ATT), min(time_GPS), min(time_RCOU)]);
t_end   = min([max(time_IMU), max(time_ATT), max(time_GPS), max(time_RCOU)]);
t_common = t_start:0.1:t_end;  % vetor de tempo com passo de 100 ms

% GPS
lat_interp = interp1(time_GPS, lat_raw, t_common, 'linear');
lon_interp = interp1(time_GPS, lon_raw, t_common, 'linear');
alt_interp = interp1(time_GPS, alt_raw, t_common, 'linear');

% IMU
gyrX_interp = interp1(time_IMU, gyrX_raw, t_common, 'linear');
gyrY_interp = interp1(time_IMU, gyrY_raw, t_common, 'linear');
gyrZ_interp = interp1(time_IMU, gyrZ_raw, t_common, 'linear');
accX_interp = interp1(time_IMU, accX_raw, t_common, 'linear');
accY_interp = interp1(time_IMU, accY_raw, t_common, 'linear');
accZ_interp = interp1(time_IMU, accZ_raw, t_common, 'linear');

% ATT
roll_interp  = interp1(time_ATT, roll_raw, t_common, 'linear');
pitch_interp = interp1(time_ATT, pitch_raw, t_common, 'linear');
yaw_interp   = interp1(time_ATT, yaw_raw, t_common, 'linear');

% RCOU
pwm1_interp = interp1(time_RCOU, pwm1_raw, t_common, 'linear');
pwm2_interp = interp1(time_RCOU, pwm2_raw, t_common, 'linear');
pwm3_interp = interp1(time_RCOU, pwm3_raw, t_common, 'linear');
pwm4_interp = interp1(time_RCOU, pwm4_raw, t_common, 'linear');

T_mr1_interp = thrustFromPWM(pwm1_interp);
T_mr2_interp = thrustFromPWM(pwm2_interp);
T_mr3_interp = thrustFromPWM(pwm3_interp);
T_mr4_interp = thrustFromPWM(pwm4_interp);

Q_mr1_interp = torqueFromPWM(pwm1_interp);
Q_mr2_interp = torqueFromPWM(pwm2_interp);
Q_mr3_interp = torqueFromPWM(pwm3_interp);
Q_mr4_interp = torqueFromPWM(pwm4_interp);

% Índices do intervalo desejado
idx_range = (t_common >= t1) & (t_common <= t2);

%Plot espaço
% Combinar os vetores de latitude, longitude e altitude em uma matriz
lla = [lat_interp(:), lon_interp(:), alt_interp(:)];
lla = lla(idx_range,:);

% Definir a origem NED como o primeiro ponto
lla0 = [lla(1,1), lla(1,2), lla(1,3)];

% Converter para coordenadas NED usando o método 'flat'
xyzNED = lla2ned(lla, lla0, 'flat');


figure;
plot(t_common(idx_range), pwm1_interp(idx_range)); hold on
plot(t_common(idx_range), pwm2_interp(idx_range)); hold on
plot(t_common(idx_range), pwm3_interp(idx_range)); hold on
plot(t_common(idx_range), pwm4_interp(idx_range));
xlabel('Tempo [s]');
ylabel('PWM');
title('PWM dos motores (trecho de tempo)');
legend('Motor1','Motor2','Motor3','Motor4', 'Interpreter', 'latex');
grid on;

% Plot do trecho filtrado
figure;
plot(t_common(idx_range), T_mr1_interp(idx_range), 'DisplayName', 'Motor 1'); hold on
plot(t_common(idx_range), T_mr2_interp(idx_range), 'DisplayName', 'Motor 2');
plot(t_common(idx_range), T_mr3_interp(idx_range), 'DisplayName', 'Motor 3');
plot(t_common(idx_range), T_mr4_interp(idx_range), 'DisplayName', 'Motor 4');
xlabel('Tempo [s]');
ylabel('Empuxo [N]');
title('Empuxo por motor (trecho de tempo)');
legend('Location','best');
grid on;

figure;
plot(t_common(idx_range), Q_mr1_interp(idx_range), 'DisplayName', 'Motor 1'); hold on
plot(t_common(idx_range), Q_mr2_interp(idx_range), 'DisplayName', 'Motor 2');
plot(t_common(idx_range), Q_mr3_interp(idx_range), 'DisplayName', 'Motor 3');
plot(t_common(idx_range), Q_mr4_interp(idx_range), 'DisplayName', 'Motor 4');
xlabel('Tempo [s]');
ylabel('Contratorque [Nm]');
title('Contratorque (trecho de tempo)');
legend('Location','best');
grid on;

figure;
plot(t_common(idx_range), gyrX_interp(idx_range), 'DisplayName', 'p'); hold on
plot(t_common(idx_range), gyrY_interp(idx_range), 'DisplayName', 'q');
plot(t_common(idx_range), gyrZ_interp(idx_range), 'DisplayName', 'r');
xlabel('Tempo [s]');
ylabel('Velocidade Angular [rad/s]');
title('Velocidade Angular (trecho de tempo)');
legend('Location','best');
grid on;

figure;
plot(t_common(idx_range), accX_interp(idx_range)); hold on
plot(t_common(idx_range), accY_interp(idx_range));
plot(t_common(idx_range), accZ_interp(idx_range));
xlabel('Tempo [s]');
ylabel('Aceleração [m/s^2]');
title('Aceleração (trecho de tempo)');
legend('$\dot{u}$','$\dot{v}$','$\dot{w}$', 'Interpreter', 'latex');
grid on;

figure;
plot(t_common(idx_range), roll_interp(idx_range)); hold on
plot(t_common(idx_range), pitch_interp(idx_range));
plot(t_common(idx_range), yaw_interp(idx_range));
xlabel('Tempo [s]');
ylabel('Atitude [°]');
title('Atitude (trecho de tempo)');
legend('roll','pitch','yaw', 'Interpreter', 'latex');
grid on;

% Dados de entrada
t = t_common(idx_range);
dt = 0.1;
N = length(t);

m = 1072.011/1000;
g = 9.81;
rx(1) = 0.311185;
ry(1) = 0.232;

rx(2) = -0.342865;
ry(2) = -0.232;

rx(3) = 0.311185;
ry(3) = -0.232;

rx(4) = -0.342865;
ry(4) = 0.232;

num_params = 20;
P0 = ones(num_params, 1) * 0.1;

Jx = 63.244/1000; % Ixx = 0.12577 Kg*m2
Jy = 250.554/1000; % Iyy = 0.08781 Kg*m2
1/Jy;
Jz = 116.192/1000; % Izz = 0.21055 Kg*m2
Jxz = 1.571/1000;  % Ixz = -0.00559 Kg*m2

gamma0 = (Jx*Jz)-Jxz^2;
gamma(1) = (Jxz*(Jx-Jy+Jz))/gamma0;
gamma(2) = (Jz*(Jz-Jy)+Jxz^2)/gamma0;
gamma(3) = Jz/gamma0;
gamma(4) = Jxz/gamma0;
gamma(5) = (Jz-Jx)/Jy;
gamma(6) = Jxz/Jy;
gamma(7) = (Jx*(Jx-Jy)+Jxz^2)/gamma0;
gamma(8) = Jx/gamma0;

% Script principal para estimação de parâmetros do VTOL

%% 1. Carregar/Definir Dados Experimentais (Substitua pelos seus dados)
disp('Carregando/Definindo dados experimentais...');
% Dados de entrada
time_vec_exp = t_common(idx_range);
N_points = length(t);

% Sinais PWM (N_points x 4) - Exemplo: entradas senoidais
% Coloque aqui seus dados de PWM. Formato: [PWM1, PWM2, PWM3, PWM4]
pwm_data = zeros(N_points, 4);
pwm_data(:,1) = pwm1_interp(idx_range)'; % PWM Motor 1 (ex: microssegundos)
pwm_data(:,2) = pwm2_interp(idx_range)'; % PWM Motor 2
pwm_data(:,3) = pwm3_interp(idx_range)'; % PWM Motor 3
pwm_data(:,4) = pwm4_interp(idx_range)'; % PWM Motor 4

% Velocidades angulares experimentais [p, q, r] (N_points x 3) - (rad/s)
% Coloque aqui seus dados de p, q, r medidos.
pqr_exp_data = zeros(N_points, 3);
pqr_exp_data(:,1) = gyrX_interp(idx_range)'; % p_exp
pqr_exp_data(:,2) = gyrY_interp(idx_range)'; % q_exp
pqr_exp_data(:,3) = gyrZ_interp(idx_range)'; % r_exp

% Condições iniciais para a simulação [p(0); q(0); r(0)]
initial_conditions = pqr_exp_data(1,:)'; % Usa o primeiro ponto de dados como condição inicial

disp('Dados experimentais definidos (USAR DADOS REAIS!).');
%% 2. Definir Funções de Referência dos Motores (A PARTIR DOS SEUS DADOS)
disp('Criando modelos de referência para empuxo e torque a partir dos dados...');
% Seus dados experimentais para o motor de referência:
pwm_values_exp = [1000; 1200; 1400; 1600; 1800; 2000]; % Coluna "PWM"
thrust_grams_exp = [0; 143; 328; 532; 784; 843];     % Coluna "Empuxo [gramas]"
torque_Nm_exp = [0.000; 0.034; 0.070; 0.115; 0.171; 0.176]; % Coluna "Contra-Torque [Nm]"

% Grau do polinômio para o ajuste (ex: 2 para quadrático)
poly_degree = 3; 

% Criar as funções de modelo usando os dados fornecidos
% Se você salvou as funções em arquivos .m separados, certifique-se que estão no path.
func_T_ref = create_thrust_model(pwm_values_exp, thrust_grams_exp, poly_degree);
func_Q_ref = create_torque_model(pwm_values_exp, torque_Nm_exp, poly_degree);
disp('Funções de referência dos motores func_T_ref e func_Q_ref definidas a partir dos dados.');

%% 3. Chutes Iniciais para os Parâmetros (P0)
% P = [G1,G2,G3,G4,G5,G6,G7,G8, invJy, gamma_ly, gamma_lx, k_T1..4, k_Q1..4] (19 parâmetros)
% Chutes mais específicos (EXEMPLOS, dependem muito da escala do seu sistema):
P0(1) = gamma(1);      % G1 ~ (Iy-Iz)/Ix
P0(2) = gamma(2);     % G2 ~ Ixz/Ix (termo Ixz, geralmente menor)
P0(3) = gamma(3);      % G3 ~ 1/Ix (coeff de Mx para p_dot)
P0(4) = gamma(4);      % G4 ~ Ixz/Ix (coeff de Mz para p_dot)
P0(5) = gamma(5);      % G5 ~ (Iz-Ix)/Iy
P0(6) = gamma(6);      % G6 ~ 1/Iy
P0(7) = gamma(7);      % G7 ~ (Ix-Iy)/Iz
P0(8) = gamma(8);     % G8 ~ Ixz/Iz
P0(9) = 1/Jy;   % invJy (se Jy ~ 0.05 kg.m^2)
P0(10) = 0.55;  % k_T1
P0(11) = 0.45;  % k_T2
P0(12) = 1.0;  % k_T3
P0(13) = 0.75;  % k_T4
P0(14) = P0(10);  % k_Q1
P0(15) = P0(11);  % k_Q2
P0(16) = P0(12);  % k_Q3
P0(17) = P0(13);  % k_Q4
P0(18) = 10;  % Dp
P0(19) = 5;  % Dq
P0(20) = 0.6;  % Dr
P0(21) = 0.7;  % Bp
P0(22) = 1.4;  % Bq
P0(23) = 0.3;  % Br

p_bias = 0.7;
q_bias = 1.4;
r_bias = 0.3;

disp('Chutes iniciais P0 definidos.');

%% 4. Limites dos Parâmetros (lb, ub) - Lower Bounds e Upper Bounds
lb = -Inf * ones(num_params, 1); % Default -Infinito
ub = Inf * ones(num_params, 1);  % Default +Infinito

lb(9) = 1e-7;
lb(10) = 1e-7;
lb(11) = 1e-7;
lb(12:19) = 1e-4;
disp('Limites lb e ub definidos.');
%% 5. Otimização usando lsqnonlin
disp('Iniciando otimização (pode levar tempo)...');
% Handle da função de custo para lsqnonlin
cost_handle = @(P_opt) cost_function_vtol(P_opt, time_vec_exp, pqr_exp_data, initial_conditions, ...
                                           time_vec_exp, pwm_data, func_T_ref, func_Q_ref, ...
                                           time_vec_exp);

% Opções de otimização
options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', 'MaxIterations', 30, ...
                       'StepTolerance', 1e-8, 'FunctionTolerance', 1e-8, ...
                       'UseParallel', true);

% Executar a otimização
[P_estimated, resnorm, residual, exitflag, output] = lsqnonlin(cost_handle, P0, lb, ub, options);
disp('Otimização finalizada.');


disp('Otimização finalizada.');
disp('Parâmetros estimados:');
disp(['G1: ', num2str(P_estimated(1))]);
disp(['G2: ', num2str(P_estimated(2))]);
disp(['G3: ', num2str(P_estimated(3))]);
disp(['G4: ', num2str(P_estimated(4))]);
disp(['G5: ', num2str(P_estimated(5))]);
disp(['G6: ', num2str(P_estimated(6))]);
disp(['G7: ', num2str(P_estimated(7))]);
disp(['G8: ', num2str(P_estimated(8))]);
disp(['InvJy: ', num2str(P_estimated(9))]);

disp(['k_T1: ', num2str(P_estimated(10))]);
disp(['k_T2: ', num2str(P_estimated(11))]);
disp(['k_T3: ', num2str(P_estimated(12))]);
disp(['k_T4: ', num2str(P_estimated(13))]);

disp(['k_Q1: ', num2str(P_estimated(14))]);
disp(['k_Q2: ', num2str(P_estimated(15))]);
disp(['k_Q3: ', num2str(P_estimated(16))]);
disp(['k_Q4: ', num2str(P_estimated(17))]);

disp(['Dp: ', num2str(P_estimated(18))]);
disp(['Dq: ', num2str(P_estimated(19))]);
disp(['Dr: ', num2str(P_estimated(20))]);

disp(['Bp: ', num2str(P_estimated(21))]);
disp(['Bq: ', num2str(P_estimated(22))]);
disp(['Br: ', num2str(P_estimated(23))]);
disp(['Norma residual final (soma dos quadrados): ', num2str(resnorm)]);
disp(['Exit flag: ', num2str(exitflag)]);
disp(output);

%% 6. Pós-Análise: Simular com Parâmetros Estimados e Plotar
disp('Simulando com parâmetros estimados...');
ode_func_estimated = @(t,y) vtol_dynamics(t, y, P_estimated, time_vec_exp, pwm_data, func_T_ref, func_Q_ref);
[t_sim, y_sim] = ode45(ode_func_estimated, time_vec_exp, initial_conditions);

p_sim = y_sim(:,1);
q_sim = y_sim(:,2);
r_sim = y_sim(:,3);

p_exp = pqr_exp_data(:,1);
q_exp = pqr_exp_data(:,2);
r_exp = pqr_exp_data(:,3);

disp('Plotando resultados...');
figure;
subplot(3,1,1);
plot(time_vec_exp', p_exp, 'b-', t_sim, p_sim, 'r--', 'LineWidth', 1.5);
legend('Experimental p', 'Simulado p');
xlabel('Tempo (s)'); ylabel('p (rad/s)');
title('Comparação da Taxa de Rolagem');
grid on;

subplot(3,1,2);
plot(time_vec_exp', q_exp, 'b-', t_sim, q_sim, 'r--', 'LineWidth', 1.5);
legend('Experimental q', 'Simulado q');
xlabel('Tempo (s)'); ylabel('q (rad/s)');
title('Comparação da Taxa de Arfagem');
grid on;

subplot(3,1,3);
plot(time_vec_exp', r_exp, 'b-', t_sim, r_sim, 'r--', 'LineWidth', 1.5);
legend('Experimental r', 'Simulado r');
xlabel('Tempo (s)'); ylabel('r (rad/s)');
title('Comparação da Taxa de Guinada');
grid on;

%sgtitle('Resultados da Estimação de Parâmetros da Dinâmica Rotacional VTOL');

disp('Script finalizado.');
