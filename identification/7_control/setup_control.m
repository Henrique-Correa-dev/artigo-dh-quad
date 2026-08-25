%SETUP_CONTROL  Popula o workspace pra simular quad_control_loops.slx
%
% Carrega o modelo linear de hover (linear_model.mat) e os ganhos do controle
% cascata (control_gains.mat), montando o que o .slx LÊ do workspace — assim o
% modelo de controle NÃO tem mais planta nem ganhos hardcoded e auto-sincroniza
% com linearize.m → design_control.m (mesma filosofia do setup_quad_v4).
%
% Variáveis populadas no base workspace:
%   A_aug, B_aug, C_aug, D_aug, X0_aug  — planta linear AUMENTADA (10 estados:
%       os 9 do corpo rígido + altitude h, com ḣ = -w). Bloco State-Space
%       "Planta_Linear".
%   K                                    — struct de ganhos (Kp_th, Kd_th, ...),
%       lido pelos blocos Gain dos 5 controladores (Kp_th='K.Kp_th' etc.).
%
% Roda automaticamente via InitFcn do modelo (antes de cada simulação) e também
% pode ser chamado à mão. NÃO faz clear/clc/close (é callback de inicialização).
%
% Uso:  >> setup_control;  sim('quad_control_loops')
%   ou simplesmente   >> sim('quad_control_loops')   (InitFcn dispara este setup)

here  = fileparts(mfilename('fullpath'));
addpath(fileparts(here));            % raiz pra setup_paths
paths = setup_paths();

%% ===================================================================
%  1) Modelo linear de hover (avisa se stale vs P_identified.mat)
%  ===================================================================
lm = load_linear_model(paths);       % A (9x9), B (9x4), u0=[m·g;0;0;0], ...
A9 = lm.A;  B9 = lm.B;

%% ===================================================================
%  2) Planta AUMENTADA com altitude h (estado 10):  ḣ = -w = -x(9)
%     (A,B do corpo rígido vêm 9×9/9×4; só a linha/coluna de h é montada aqui)
%  ===================================================================
A_aug  = [A9, zeros(9,1); zeros(1,8), -1, 0];
B_aug  = [B9; zeros(1,4)];
C_aug  = eye(10);
D_aug  = zeros(10,4);
X0_aug = zeros(10,1);

%% ===================================================================
%  3) Ganhos do controle cascata (de design_control.m)
%  ===================================================================
G = load(fullfile(paths.control, 'control_gains.mat'));
K = G.K;

% Stale guard leve: os ganhos têm de ter sido projetados pra esta massa.
proj_p = parameters();
if isfield(K,'m') && abs(K.m - proj_p.m) > 1e-6
    warning('setup_control:staleGains', ...
        ['control_gains.mat foi projetado com m=%.4f, mas parameters().m=%.4f.\n' ...
         '         >>> Rode design_control.m pra reprojetar os ganhos. <<<'], ...
        K.m, proj_p.m);
end

fprintf('setup_control: planta aumentada (10×10) + ganhos K carregados.\n');
fprintf('  polos da planta (real): '); fprintf('%+.3f ', sort(real(eig(A_aug)))); fprintf('\n');
fprintf('  Kp/Kd: pitch=%.3f/%.3f roll=%.3f/%.3f yaw=%.3f/%.3f | Kp_h=%.3f\n', ...
        K.Kp_th, K.Kd_th, K.Kp_ph, K.Kd_ph, K.Kp_ps, K.Kd_ps, K.Kp_h);
