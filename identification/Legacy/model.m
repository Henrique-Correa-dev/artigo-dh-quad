%[text] ## PASSO 1: Recriação da Atitude com Dados (ODE45)
% SCRIPT PARA CALCULAR ATITUDE A PARTIR DE DADOS DE VOO USANDO ODE45
% Este script simula a atitude (phi, theta, psi) integrando as equações
% cinemáticas com dados experimentais de p, q, r.

% 1. DATA PREPARATION
t_raw = t_common(idx);
dt = 0.1;

p_data = gyrX_interp(idx_range);
q_data = gyrY_interp(idx_range);
r_data = gyrZ_interp(idx_range);

phi_real = roll_interp(idx_range);
theta_real = pitch_interp(idx_range);
psi_real = yaw_interp(idx_range);

% Cria o vetor de tempo para a simulação, começando em zero.
% Isso é uma boa prática para solvers ODE.
time_vector = t_raw - t_raw(1);

% Para comparar no SIMULINK
p_in = [time_vector', p_data'];
q_in = [time_vector', q_data'];
r_in = [time_vector', r_data'];

phi_out = [time_vector', phi_real'];
theta_out = [time_vector', theta_real'];
psi_out = [time_vector', psi_real'];

phi0 = [time_vector',deg2rad(phi_real(1))*ones(100,1)];
theta0 = [time_vector',deg2rad(theta_real(1))*ones(100,1)];
psi0 = [time_vector',deg2rad(psi_real(1))*ones(100,1)];

%% 2. SIMULATION SETUP AND EXECUTION
% Condições iniciais de atitude [phi; theta; psi] em radianos
initial_conditions = [phi0(1,2); theta0(1,2); psi0(1,2)]; 

% Intervalo de tempo para a integração
time_span = [time_vector(1), time_vector(end)];

% Opções do solver para alta precisão
solver_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

disp('Iniciando a simulação com ode45...');
% Chamada do solver usando uma função anônima para passar os dados
[t_sim, y_sim] = ode45(@(t,y) euler_kinematics_ode(t, y, time_vector, p_data, q_data, r_data), time_span, initial_conditions, solver_options);
disp('Simulação com ode45 concluída.');

%% 3. RESULTS PROCESSING AND VISUALIZATION
% Extrai os ângulos e converte para graus
phi_sim   = rad2deg(y_sim(:,1));
theta_sim = rad2deg(y_sim(:,2));
psi_sim   = rad2deg(y_sim(:,3));

% Plot dos resultados com comparação
figure('Name', 'Comparação: Atitude Simulada vs. Real');
sgtitle('Validação da Cinemática: Simulado (ODE45) vs. Dados Reais');

subplot(3,1,1);
plot(t_sim, phi_sim, 'b-', 'LineWidth', 2);
hold on;
plot(time_vector, phi_real, 'k--', 'LineWidth', 1.5);
hold off;
grid on; legend('Simulado (ODE45)', 'Dados Reais');
ylabel('Rolagem \phi (graus)');

subplot(3,1,2);
plot(t_sim, theta_sim, 'r-', 'LineWidth', 2);
hold on;
plot(time_vector, theta_real, 'k--', 'LineWidth', 1.5);
hold off;
grid on; legend('Simulado (ODE45)', 'Dados Reais');
ylabel('Arfagem \theta (graus)');

subplot(3,1,3);
plot(t_sim, psi_sim, 'g-', 'LineWidth', 2);
hold on;
plot(time_vector, psi_real, 'k--', 'LineWidth', 1.5);
hold off;
grid on; legend('Simulado (ODE45)', 'Dados Reais');
ylabel('Guinada \psi (graus)');
xlabel('Tempo (s)');
%% 4. ORDINARY DIFFERENTIAL EQUATION (ODE) FUNCTION
% Função local chamada pelo ode45.
function dydt = euler_kinematics_ode(t, y, time_data, p_data, q_data, r_data)
    % Desempacotar estados de atitude: y = [phi; theta; psi]
    phi   = y(1);
    theta = y(2);

    % Interpolar os dados de entrada para o tempo 't' exato do solver
    p = interp1(time_data, p_data, t);
    q = interp1(time_data, q_data, t);
    r = interp1(time_data, r_data, t);

    % Proteção contra a singularidade de Gimbal Lock
    cos_theta = cos(theta);
    if abs(cos_theta) < 1e-7
        cos_theta = 1e-7 * sign(cos_theta); 
    end
    
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    tan_theta = sin(theta) / cos_theta;

    % Equações da cinemática
    phi_dot   = p + q * sin_phi * tan_theta + r * cos_phi * tan_theta;
    theta_dot =     q * cos_phi             - r * sin_phi;
    psi_dot   =   ( q * sin_phi + r * cos_phi ) / cos_theta;

    % Retornar o vetor de derivadas [phi_dot; theta_dot; psi_dot]
    dydt = [phi_dot; theta_dot; psi_dot];
end
%%
%[text] ## PASSO 2: Recriação da Posição com Dados (ODE45)
% SCRIPT PARA CALCULAR POSIÇÃO A PARTIR DE DADOS DE VOO USANDO ODE45
% Este script simula a posição (pe, pn, pd) integrando as equações
% cinemáticas com dados experimentais de p, q, r.

% 1. DATA PREPARATION
t_raw = t_common(idx_range);
dt = 0.1;

% --- Constantes Físicas (SUBSTITUA PELOS SEUS VALORES) ---
constants.m = 2385.011/1000;      % Massa do VTOL em kg
constants.g = 9.81;    % Aceleração da gravidade em m/s^2
constants.func_T_ref = func_T_ref;

p_data = gyrX_interp(idx_range);
q_data = gyrY_interp(idx_range);
r_data = gyrZ_interp(idx_range);

pwm1_data = pwm1_interp(idx_range);
pwm2_data = pwm2_interp(idx_range);
pwm3_data = pwm3_interp(idx_range);
pwm4_data = pwm4_interp(idx_range);

% DADOS DE ATITUDE SÃO NECESSÁRIOS PARA CALCULAR A GRAVIDADE NO CORPO
% Assumimos que você os tem do passo anterior (em graus)
phi_data_deg = roll_interp(idx_range);
theta_data_deg = pitch_interp(idx_range);
psi_data_deg = yaw_interp(idx_range);

fprintf('Verificando dados de atitude por valores NaN...\n');
if any(isnan(phi_data_deg)), fprintf('AVISO: NaN encontrado em phi_data_deg (roll_interp)\n'); end
if any(isnan(theta_data_deg)), fprintf('AVISO: NaN encontrado em theta_data_deg (pitch_interp)\n'); end
if any(isnan(psi_data_deg)), fprintf('AVISO: NaN encontrado em psi_data_deg (yaw_interp)\n'); end
fprintf('Verificação concluída.\n\n');

% Dados reais para validação
accel_real.x = accX_interp(idx_range); % Aceleração real no corpo
accel_real.y = accY_interp(idx_range);
accel_real.z = accZ_interp(idx_range);
pos_real.x = Pn; % Posição real no mundo
pos_real.y = E;
pos_real.z = -D;

% Cria o vetor de tempo para a simulação, começando em zero.
% Isso é uma boa prática para solvers ODE.
time_vector = t_raw - t_raw(1);

% Para comparar no SIMULINK
min_pwm_for_active_thrust = 1000;
pwm1_in = [time_vector', pwm1_data'];
pwm2_in = [time_vector', pwm2_data'];
pwm3_in = [time_vector', pwm3_data'];
pwm4_in = [time_vector', pwm4_data'];

u_dot_out = [time_vector', accX_interp(idx_range)'];
v_dot_out = [time_vector', accY_interp(idx_range)'];
w_dot_out = [time_vector', accZ_interp(idx_range)'];

p_out = [time_vector', gyrX_interp(idx_range)'];
q_out = [time_vector', gyrY_interp(idx_range)'];
r_out = [time_vector', gyrZ_interp(idx_range)'];

x_out = [time_vector', Pn];
y_out = [time_vector', E];
z_out = [time_vector', -D];

% Agrupa todos os dados de séries temporais em uma única struct para facilitar
time_series_data.time  = time_vector;
time_series_data.p     = p_data;
time_series_data.q     = q_data;
time_series_data.r     = r_data;
time_series_data.phi   = deg2rad(phi_data_deg); % ODE precisa de radianos
time_series_data.theta = deg2rad(theta_data_deg);
time_series_data.psi   = deg2rad(psi_data_deg);
time_series_data.pwm1  = pwm1_data;
time_series_data.pwm2  = pwm2_data;
time_series_data.pwm3  = pwm3_data;
time_series_data.pwm4  = pwm4_data;

%% 2. SIMULATION SETUP AND EXECUTION
% Condições iniciais para os 6 estados: [u, v, w, x, y, z]
% Idealmente, use os primeiros pontos dos seus dados reais.
% u,v,w podem ser 0 se o voo começar do repouso.
initial_conditions = [0; 0; 0]; 
coeffs_thrust = [-1.2349e-08, 5.4906e-05, -0.06998, 27.488];

% Intervalo e opções do solver
time_span = [time_vector(1), time_vector(end)];
solver_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);

disp('Iniciando a simulação da dinâmica de translação...');
[t_sim, y_sim] = ode45(@(t,y) translational_dynamics_ode(t, y, constants, time_series_data), time_span, initial_conditions, solver_options);
disp('Simulação concluída.');

%% 3. RESULTS PROCESSING AND VISUALIZATION
% Extrair estados simulados
u_sim = y_sim(:,1); v_sim = y_sim(:,2); w_sim = y_sim(:,3);
% x_sim = y_sim(:,4); y_sim = y_sim(:,5); z_sim = y_sim(:,6);

% --- Recalcular acelerações simuladas para comparação ---
% Para plotar a aceleração, precisamos recalculá-la para cada ponto da simulação
accel_sim_x = zeros(length(t_sim), 1);
accel_sim_y = zeros(length(t_sim), 1);
accel_sim_z = zeros(length(t_sim), 1);
for k = 1:length(t_sim)
    dydt_k = translational_dynamics_ode(t_sim(k), y_sim(k,:)', constants, time_series_data);
    accel_sim_x(k) = dydt_k(1);
    accel_sim_y(k) = dydt_k(2);
    accel_sim_z(k) = dydt_k(3);
end

% --- Visualização ---
% Comparação das Acelerações no Referencial do Corpo
figure('Name', 'Validação das Acelerações');
sgtitle('Aceleração Simulada vs. Real (no Corpo)');
subplot(3,1,1); plot(t_sim, accel_sim_x, 'b', time_vector, accel_real.x, 'k--'); legend('Simulado', 'Real'); ylabel('u_{dot} (m/s^2)'); grid on;
subplot(3,1,2); plot(t_sim, accel_sim_y, 'r', time_vector, accel_real.y, 'k--'); legend('Simulado', 'Real'); ylabel('v_{dot} (m/s^2)'); grid on;
subplot(3,1,3); plot(t_sim, accel_sim_z-9.81, 'g', time_vector, accel_real.z, 'k--'); legend('Simulado', 'Real'); ylabel('w_{dot} (m/s^2)'); xlabel('Tempo (s)'); grid on;
%% 4. FUNÇÃO DA EQUAÇÃO DIFERENCIAL ORDINÁRIA (EDO)
function dydt = translational_dynamics_ode(t, y, constants, data)
    % --- Desempacotar estados e constantes ---
    u = y(1); v = y(2); w = y(3);
    %x = y(4); y = y(5); z = y(6); % Posições não são necessárias para o cálculo das derivadas

    m = constants.m;
    g = constants.g;
    func_T_ref = constants.func_T_ref;

    % --- Interpolar todas as entradas de séries temporais ---
    p     = interp1(data.time, data.p, t);
    q     = interp1(data.time, data.q, t);
    r     = interp1(data.time, data.r, t);
    phi   = interp1(data.time, data.phi, t);
    theta = interp1(data.time, data.theta, t);
    psi   = interp1(data.time, data.psi, t);
    pwm1  = interp1(data.time, data.pwm1, t);
    pwm2  = interp1(data.time, data.pwm2, t);
    pwm3  = interp1(data.time, data.pwm3, t);
    pwm4  = interp1(data.time, data.pwm4, t);

    % --- Cálculo das Forças (no corpo) ---
    % 1. Empuxo Total (Thrust)
    T_mr1 = func_T_ref(pwm1);
    T_mr2 = func_T_ref(pwm2);
    T_mr3 = func_T_ref(pwm3);
    T_mr4 = func_T_ref(pwm4);
    Tmr = T_mr1 + T_mr2 + T_mr3 + T_mr4;
    Thrust_force_body = [0; 0; -Tmr]; % Empuxo atua na direção -z do corpo

    % 2. Força da Gravidade
    R_bn = [ ... % Matriz de rotação Corpo -> Mundo (Body to Navigation)
        cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
        cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
        -sin(theta),         sin(phi)*cos(theta),                            cos(phi)*cos(theta)
    ];
    Gravity_force_body = R_bn' * [0; 0; m*g]; % Transposta rotaciona do Mundo -> Corpo

    % 3. Força Total no Corpo
    F_body = Thrust_force_body + Gravity_force_body;
    fx = F_body(1); fy = F_body(2); fz = F_body(3);
    
    % --- Cálculo das Derivadas ---
    % 1. Derivadas das velocidades lineares (acelerações no corpo)
    u_dot = r*v - q*w + fx/m;
    v_dot = p*w - r*u + fy/m;
    w_dot = q*u - p*v + fz/m;

    % 2. Derivadas da posição no mundo
    %world_velocity = R_bn * [u; v; w];
    %x_dot = world_velocity(1);
    %y_dot = world_velocity(2);
    %z_dot = world_velocity(3);

    % --- Montar o vetor de saída com as 6 derivadas ---
    dydt = [u_dot; v_dot; w_dot];
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
