% probe_coupling.m — varredura p/ achar o limiar de divergência linear↔NL
% por ACOPLAMENTO GIROSCÓPICO. Manobra: spin de guinada contínuo (sustenta r) +
% UM degrau de velocidade de avanço em t=3s (pulso de q) → o termo -Γ2·qr induz
% rolagem φ, que o modelo linear (desacoplado) prevê ≡ 0. Mantém amplitudes
% moderadas p/ ficar no regime ESTÁVEL (sem tombamento).
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths=setup_paths();
lm=load(fullfile(paths.outputs,'linear_model.mat'));
G =load(fullfile(paths.control,'control_gains.mat'));
proj=parameters(); [fT,fQ]=motor_models(); dyn=vtol_dynamics('get_handles');
M=struct('A',lm.A,'B',lm.B,'u0',lm.u0,'P',lm.P,'K',G.K,'dyn',dyn, ...
         'bridge',forces_to_pwm(lm,'nonlinear'),'fT',fT,'fQ',fQ, ...
         'constants',struct('m',proj.m,'g',proj.g),'dt',0.01,'PWM_LIM',[1200 2000]);

yaw_rates=[0 15 30 45 60 75 90];   % °/s (spin contínuo)
u_step=3;                          % degrau de vel. de avanço em t=3s (pulso de q)
T_end=10; t=(0:M.dt:T_end)';

fprintf('yawrate |  Δφ    Δθ    Δψ    Δh   |  φ_NLmáx θ_NLmáx (°)  estável?\n');
fprintf('--------+------------------------+-------------------------------\n');
for i=1:numel(yaw_rates)
    yr=yaw_rates(i);
    s=struct('name','coup','h0',3,'t',t, ...
        'sp_h',3*ones(size(t)),'sp_u',u_step*(t>=3),'sp_psi',deg2rad(yr)*t);
    RL=cl_loop(false,s,M); RN=cl_loop(true,s,M);
    dphi=rad2deg(max(abs(RL.X(:,4)-RN.X(:,4))));
    dth =rad2deg(max(abs(RL.X(:,5)-RN.X(:,5))));
    dpsi=rad2deg(max(abs(RL.X(:,6)-RN.X(:,6))));
    dh  =max(abs(RL.H-RN.H));
    phiM=rad2deg(max(abs(RN.X(:,4)))); thM=rad2deg(max(abs(RN.X(:,5))));
    stab = (phiM<60 && thM<60 && dh<2);
    fprintf('%5d°  | %5.2f %5.2f %5.2f %5.3f | %7.1f %7.1f      %s\n', ...
        yr,dphi,dth,dpsi,dh,phiM,thM,string(stab));
end
fprintf('(Δ = max|linear - não-linear|; φ,θ,ψ em °, h em m)\n');
