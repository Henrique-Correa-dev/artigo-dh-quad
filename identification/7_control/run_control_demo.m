% run_control_demo.m — Roda TODA a etapa 7_control em sequência.
%
%   1) design_control        → projeta os ganhos (control_gains.mat)
%   2) sim_control_linear     → malha fechada no MODELO LINEAR  (figura)
%   3) sim_control_nonlinear  → malha fechada na PLANTA NÃO-LINEAR vs linear
%   4) Simulink + 3D          → quad_control_3d.slx, abre e simula
%                               (Mechanics Explorer anima o drone do .stp)
%
% Uso:  >> run_control_demo
%
% Obs.: os scripts de simulação fazem 'clear' no início (por isso este
% orquestrador NÃO guarda variáveis entre eles — chama tudo por nome, já
% que a pasta 7_control está no path via setup_paths).

clear; clc; close all;
addpath(fileparts(mfilename('fullpath')));      % garante 7_control no path
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();

fprintf('\n############ 1/4  PROJETO DE CONTROLE ############\n');
design_control

fprintf('\n############ 2/4  MALHA FECHADA LINEAR ############\n');
sim_control_linear

fprintf('\n############ 3/4  MALHA FECHADA NÃO-LINEAR ############\n');
sim_control_nonlinear

fprintf('\n############ 4/4  SIMULINK (UM BLOCO POR CONTROLADOR) ############\n');
open_system('quad_control_loops');                 % modelo já pronto (.slx)
sim('quad_control_loops');
fprintf('\n>>> Pronto. O modelo "quad_control_loops" tem um bloco por malha:\n');
fprintf('    Ctrl_Altitude | Ctrl_VelFwd | Ctrl_Pitch | Ctrl_Roll | Ctrl_Heading.\n');
