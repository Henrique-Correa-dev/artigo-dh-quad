%BUILD_QUAD_MODEL_LINEAR  Constrói (programaticamente) o Simulink linear de hover.
%
% Gera quad_model_linear.slx: modelo de espaço-de-estados linearizado no
% TRIM DE HOVER (equilíbrio), na forma clássica, entrada em FORÇAS:
%   δẋ = A·δx + B·δv,   δv = [T,Mx,My,Mz] - v0,   x = x0 + δx
% (sem termo afim — no equilíbrio c=0.)
%
%   Estados x = [p q r phi theta psi u v w]  (9)
%   Entrada v = [T Mx My Mz]  (empuxo total + 3 momentos)   (4)
%   Saídas: p_sim,q_sim,r_sim (rad/s) | phi_sim,theta_sim,psi_sim (deg) |
%           u_sim,v_sim,w_sim (m/s)
%
% A entrada [T,M] é o que um controlador (atitude/altitude) comanda direto.
% Pra validar contra voo, alimente [T,M](t) derivado do PWM logado (a alocação,
% via vtol_dynamics forces handle) — é o que setup_quad_linear.m faz.
%
% Parâmetros vêm do workspace (popule com setup_quad_linear.m):
%   A_lin (9x9), B_lin (9x4), x0_lin (9x1), u0_lin (4x1=[mg;0;0;0]), dx0_lin (9x1)
%
% Uso:  >> build_quad_model_linear   % cria/salva o .slx

mdl = 'quad_model_linear';
here = fileparts(mfilename('fullpath'));
slx_path = fullfile(here, [mdl '.slx']);

% Recomeça do zero se já existir
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(slx_path, 'file'), delete(slx_path); end

new_system(mdl);
load_system(mdl);

L = 'simulink/';
add = @(src, name, pos, varargin) add_block([L src], [mdl '/' name], ...
        'Position', pos, varargin{:});

%% --- Entradas: 4 From Workspace [T, Mx, My, Mz] → Mux (4-vetor) ---
add('Sources/From Workspace', 'T_in',  [20  20  90  40], 'VariableName','T_in',  'SampleTime','0');
add('Sources/From Workspace', 'Mx_in', [20  60  90  80], 'VariableName','Mx_in', 'SampleTime','0');
add('Sources/From Workspace', 'My_in', [20 100  90 120], 'VariableName','My_in', 'SampleTime','0');
add('Sources/From Workspace', 'Mz_in', [20 140  90 160], 'VariableName','Mz_in', 'SampleTime','0');
add('Signal Routing/Mux', 'MuxForce', [140 20 145 160], 'Inputs','4');

%% --- δv = [T,M] - v0 ---
add('Sources/Constant', 'v0', [140 210 210 240], 'Value','u0_lin');
add('Math Operations/Sum', 'SubV0', [250 80 270 100], 'Inputs','+-');

%% --- State-Space: δẋ = A·δx + B·δv, y = δx (forma clássica, hover) ---
add('Continuous/State-Space', 'LinSS', [390 80 510 150], ...
    'A','A_lin', 'B','B_lin', 'C','eye(9)', 'D','zeros(9,4)', 'X0','dx0_lin');

%% --- x = δx + x0 ---
add('Sources/Constant', 'x0', [390 210 460 240], 'Value','x0_lin');
add('Math Operations/Sum', 'AddX0', [560 100 580 120], 'Inputs','++');

%% --- Demux 9 estados ---
add('Signal Routing/Demux', 'Demux', [630 30 635 410], 'Outputs','9');

%% --- Conversão rad→deg pros ângulos ---
add('Math Operations/Gain', 'rad2deg_phi',   [700 150 740 180], 'Gain','180/pi');
add('Math Operations/Gain', 'rad2deg_theta', [700 190 740 220], 'Gain','180/pi');
add('Math Operations/Gain', 'rad2deg_psi',   [700 230 740 260], 'Gain','180/pi');

%% --- To Workspace (9 saídas, MaxDataPoints=inf pra não truncar) ---
sig = {'p_sim','q_sim','r_sim','phi_sim','theta_sim','psi_sim','u_sim','v_sim','w_sim'};
ypos = 30:42:30+42*8;
for i = 1:9
    add('Sinks/To Workspace', sig{i}, [820 ypos(i) 900 ypos(i)+20], ...
        'VariableName', sig{i}, 'SaveFormat','Timeseries', 'MaxDataPoints','inf');
end

%% --- Conexões ---
add_line(mdl, 'T_in/1',  'MuxForce/1', 'autorouting','on');
add_line(mdl, 'Mx_in/1', 'MuxForce/2', 'autorouting','on');
add_line(mdl, 'My_in/1', 'MuxForce/3', 'autorouting','on');
add_line(mdl, 'Mz_in/1', 'MuxForce/4', 'autorouting','on');
add_line(mdl, 'MuxForce/1', 'SubV0/1', 'autorouting','on');
add_line(mdl, 'v0/1',       'SubV0/2', 'autorouting','on');
add_line(mdl, 'SubV0/1',  'LinSS/1', 'autorouting','on');
add_line(mdl, 'LinSS/1',  'AddX0/1', 'autorouting','on');
add_line(mdl, 'x0/1',     'AddX0/2', 'autorouting','on');
add_line(mdl, 'AddX0/1',  'Demux/1', 'autorouting','on');

% Demux → saídas (4,5,6 passam por rad2deg)
straight = [1 2 3 7 8 9];
for k = straight
    add_line(mdl, sprintf('Demux/%d', k), sprintf('%s/1', sig{k}), 'autorouting','on');
end
add_line(mdl, 'Demux/4', 'rad2deg_phi/1',   'autorouting','on');
add_line(mdl, 'Demux/5', 'rad2deg_theta/1', 'autorouting','on');
add_line(mdl, 'Demux/6', 'rad2deg_psi/1',   'autorouting','on');
add_line(mdl, 'rad2deg_phi/1',   'phi_sim/1',   'autorouting','on');
add_line(mdl, 'rad2deg_theta/1', 'theta_sim/1', 'autorouting','on');
add_line(mdl, 'rad2deg_psi/1',   'psi_sim/1',   'autorouting','on');

%% --- Solver (igual ao v4: fixed-step ode4, dt=0.01) ---
set_param(mdl, 'SolverType','Fixed-step', 'Solver','ode4', 'FixedStep','0.01');

%% --- Salvar ---
save_system(mdl, slx_path);
fprintf('quad_model_linear.slx criado: %s\n', slx_path);
fprintf('  %d blocos | State-Space 9 estados | forma clássica δẋ=Aδx+Bδu (hover)\n', ...
    numel(find_system(mdl,'SearchDepth',1,'Type','Block')));
close_system(mdl, 0);
