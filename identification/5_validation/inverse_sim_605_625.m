%INVERSE_SIM_605_625  Inverse Simulation (Jategaonkar 2015, Sec. 11.4, Fig. 11.1)
%  do canal ROTACIONAL, janela 605-625 s.
%
%  Ideia (≠ Output-Error): alimenta o PWM medido no modelo e fecha uma malha de
%  feedback (PI) que força a resposta do modelo (p,q,r) a bater EXATAMENTE com a
%  medida (gyro). A deficiência de modelagem migra pra dentro do controle:
%
%       ΔM = [ΔMx, ΔMy, ΔMz]  = momento EXTRA que o modelo precisou.
%
%  Interpretação (Jategaonkar): ΔM pequeno e centrado em zero → modelo bom.
%  ΔM estruturado (tendência/oscilação) → erro de modelagem NAQUELE eixo, e
%  sugere a extensão de modelo necessária. Mantém a condição de trim.
%
%  Feedback nos MOMENTOS (não no PWM) → interpretação por eixo é direta.

clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath')))); setup_paths();
paths = setup_paths();

% ---------------------------------------------------------------- config
LOG_FILE = 'logs_concat.mat';
t_window = [605, 625];
DELAY_PWM = 1;                 % = identify_plant/validate (alinha sim×medido)
TAU_C = 0.15;                  % [s] constante de tempo alvo da malha (tracking)
TAU_I = 0.60;                  % [s] tempo integral do PI

% P_J (modelo não linear) — mesmo do validate_605_625
P_J = [ 0.050326; 0.097063; 0.126192; 0.001571; ...   % Jx Jy Jz Jxz
        1.06; 1.06; 0.91; 0.93;                  ...   % k_T1..4
        0.65; 0.63; 0.75; 0.67;                  ...   % k_Q1..4
        5.69; 3.94; 0.74];                             % Dp Dq Dr

% ---------------------------------------------------------------- parâmetros do modelo
Pe = P_J_to_simulink(P_J);
G1=Pe(1);G2=Pe(2);G3=Pe(3);G4=Pe(4);G5=Pe(5);G6=Pe(6);G7=Pe(7);G8=Pe(8);
invJy=Pe(9); Dp=Pe(18); Dq=Pe(19); Dr=Pe(20);
kT=P_J(5:8); kQ=P_J(9:12);
pp=parameters();
Lxr=pp.arms.Lx_r; Lxl=pp.arms.Lx_l; Lyf=pp.arms.Ly_f; Lyr=pp.arms.Ly_r;
fT=@(x)max(0,interp1(pp.bench.pwm, pp.bench.T_grams*9.80665/1000, x,'makima','extrap'));
fQ=@(x)max(0,interp1(pp.bench.pwm, pp.bench.Q_Nm,                 x,'makima','extrap'));

% ---------------------------------------------------------------- log + janela
L = load_log_data(fullfile(paths.data, LOG_FILE));
t_lo=max([min(L.time_IMU) min(L.time_RCOU)]); t_hi=min([max(L.time_IMU) max(L.time_RCOU)]);
tc=(t_lo:0.1:t_hi)'; idx=(tc>=t_window(1))&(tc<=t_window(2)); t=tc(idx); N=numel(t); dt=0.1;
pwm=[interp1(L.time_RCOU,L.pwm1_raw,tc),interp1(L.time_RCOU,L.pwm2_raw,tc),...
     interp1(L.time_RCOU,L.pwm3_raw,tc),interp1(L.time_RCOU,L.pwm4_raw,tc)];
if DELAY_PWM>0, dl=@(x)[repmat(x(1),DELAY_PWM,1);x(1:end-DELAY_PWM)]; for c=1:4,pwm(:,c)=dl(pwm(:,c));end; end
pwm=pwm(idx,:);
pqr=[interp1(L.time_IMU,L.gyrX_raw,tc),interp1(L.time_IMU,L.gyrY_raw,tc),interp1(L.time_IMU,L.gyrZ_raw,tc)];
pqr=pqr(idx,:);

% ---------------------------------------------------------------- momentos do modelo (PWM→M)
T=(kT(:)').*fT(pwm); Q=(kQ(:)').*fQ(pwm);            % Nx4
Mx_m=-(Lxr*T(:,1)-Lxl*T(:,2)-Lxl*T(:,3)+Lxr*T(:,4));
My_m=  Lyf*T(:,1)-Lyr*T(:,2)+Lyf*T(:,3)-Lyr*T(:,4);
Mz_m=  Q(:,1)+Q(:,2)-Q(:,3)-Q(:,4);

% ---------------------------------------------------------------- ganhos PI (auto-tune pela efetividade)
% ∂ṗ/∂Mx=G3 , ∂q̇/∂My=invJy , ∂ṙ/∂Mz=G8  → Kp = 1/(TAU_C·B)
Kp=[1/(TAU_C*G3), 1/(TAU_C*invJy), 1/(TAU_C*G8)];
Ki=Kp/TAU_I;

% ---------------------------------------------------------------- malha inversa (PI força p,q,r ≡ medido)
p=pqr(1,1); q=pqr(1,2); r=pqr(1,3); I=[0 0 0];
dM=zeros(N,3); pqr_mod=zeros(N,3); nsub=10; h=dt/nsub;
for k=1:N
    pqr_mod(k,:)=[p q r];
    e=pqr(k,:)-[p q r];                  % resíduo de SAÍDA (z - y)
    I=I+e*dt;
    dMk=Kp.*e + Ki.*I;                   % PI → correção de momento (ΔM)
    dM(k,:)=dMk;
    Mx=Mx_m(k)+dMk(1); My=My_m(k)+dMk(2); Mz=Mz_m(k)+dMk(3);
    for s=1:nsub                         % integra o modelo dt à frente
        pd=G1*p*q-G2*q*r+G3*Mx+G4*Mz-Dp*p;
        qd=G5*p*r-G6*(p^2-r^2)+invJy*My-Dq*q;
        rd=G7*p*q-G1*q*r+G4*Mx+G8*Mz-Dr*r;
        p=p+h*pd; q=q+h*qd; r=r+h*rd;
    end
end

% ---------------------------------------------------------------- métricas
R2=@(y,yh)1-sum((y-yh).^2)/sum((y-mean(y)).^2);
trk=[R2(pqr(:,1),pqr_mod(:,1)) R2(pqr(:,2),pqr_mod(:,2)) R2(pqr(:,3),pqr_mod(:,3))];
rmsM =@(x)sqrt(mean(x.^2));
fprintf('\n=== INVERSE SIMULATION (Jategaonkar §11.4) — janela %g-%g s ===\n',t_window);
fprintf('  Tracking do modelo (deve ~1): R² p=%.3f q=%.3f r=%.3f\n',trk);
fprintf('  ΔM (deficiência): eixos | RMS(ΔM) | RMS(M_modelo) | razão %%\n');
nm={'Mx (roll)','My (pitch)','Mz (yaw) '}; Mm=[Mx_m My_m Mz_m];
for a=1:3
  fprintf('    %s | %7.4f N·m | %7.4f N·m | %5.1f%%  (média ΔM=%+.4f)\n', ...
    nm{a}, rmsM(dM(:,a)), rmsM(Mm(:,a)), 100*rmsM(dM(:,a))/max(rmsM(Mm(:,a)),1e-9), mean(dM(:,a)));
end

% ---------------------------------------------------------------- plot ΔM (a deficiência)
fig=figure('Position',[60 60 1150 720],'Color','w');
ax_t=t-t(1); ttl={'\DeltaM_x  (roll)','\DeltaM_y  (pitch)','\DeltaM_z  (yaw)'};
for a=1:3
  subplot(3,1,a); hold on; grid on;
  plot(ax_t, Mm(:,a),'Color',[.6 .6 .6],'LineWidth',.8);          % M_modelo (escala)
  plot(ax_t, dM(:,a),'r','LineWidth',1.3);                        % ΔM (deficiência)
  yline(0,'k:');
  ylabel('N·m'); legend({'M_{modelo}','\DeltaM (deficiência)'},'Location','best');
  title(sprintf('%s   |   RMS(\\DeltaM)=%.4f N·m  (%.1f%% de M_{mod})', ...
        ttl{a}, rmsM(dM(:,a)), 100*rmsM(dM(:,a))/max(rmsM(Mm(:,a)),1e-9)));
end
xlabel('t [s]');
sgtitle('Inverse Simulation: \DeltaM = controle extra que o modelo exigiu (deficiência de modelagem)');
out=fullfile(paths.images,'inverse_sim_605_625.png');
if ~exist(paths.images,'dir'),mkdir(paths.images);end
exportgraphics(fig,out,'Resolution',130);
fprintf('\nFigura: %s\n',out);
