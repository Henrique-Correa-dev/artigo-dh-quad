% Varredura de tau_m: isola o efeito do LAG (tau->0 = só curva RPM²) no R² de p,q,r.
addpath(fileparts(fileparts(mfilename('fullpath')))); setup_paths();
p=parameters(); L=load_log_data(fullfile(setup_paths().data,'logs_concat.mat'));
P=[0.050326;0.097063;0.126192;0.001571; 1.06;1.06;0.91;0.93; 0.65;0.63;0.75;0.67; 5.69;3.94;0.74];
Jx=P(1);Jy=P(2);Jz=P(3);Jxz=P(4);Dp=P(13);Dq=P(14);Dr=P(15);
g0=Jx*Jz-Jxz^2; G=[Jxz*(Jx-Jy+Jz)/g0;(Jz*(Jz-Jy)+Jxz^2)/g0;Jz/g0;Jxz/g0;(Jz-Jx)/Jy;Jxz/Jy;(Jx*(Jx-Jy)+Jxz^2)/g0;Jx/g0;1/Jy];
mo=p.motor; ar=p.arms; kt=P(5:8); kq=P(9:12);
tlo=max(min(L.time_IMU),min(L.time_RCOU)); thi=min(max(L.time_IMU),max(L.time_RCOU));
tc=(tlo:0.1:thi)'; id=tc>=610&tc<=620; t=tc(id);
pw=[interp1(L.time_RCOU,L.pwm1_raw,tc),interp1(L.time_RCOU,L.pwm2_raw,tc),interp1(L.time_RCOU,L.pwm3_raw,tc),interp1(L.time_RCOU,L.pwm4_raw,tc)];
dl=@(x)[x(1);x(1:end-1)]; for c=1:4,pw(:,c)=dl(pw(:,c));end; pw=pw(id,:);
pqr=[interp1(L.time_IMU,L.gyrX_raw,tc),interp1(L.time_IMU,L.gyrY_raw,tc),interp1(L.time_IMU,L.gyrZ_raw,tc)]; pqr=pqr(id,:);
R2=@(y,yh)1-sum((y-yh).^2)/sum((y-mean(y)).^2);
fprintf('\ntau_m    R2_p     R2_q     R2_r\n');
for tau=[0.001 0.02 0.05 0.08 0.12 0.20]
  ode=@(tt,y) rotmotor(tt,y,t,pw,mo,kt,kq,ar,G,Dp,Dq,Dr,tau);
  w0=max(0,mo.CR*pw(1,:)+mo.Omega_b); y0=[pqr(1,:)';w0(:)];
  [ts,ys]=ode45(ode,t,y0,odeset('RelTol',1e-6,'AbsTol',1e-9)); yo=interp1(ts,ys,t);
  fprintf('%.3f   %+.4f  %+.4f  %+.4f\n',tau,R2(pqr(:,1),yo(:,1)),R2(pqr(:,2),yo(:,2)),R2(pqr(:,3),yo(:,3)));
end

function dy=rotmotor(tt,y,tg,pw,mo,kt,kq,ar,G,Dp,Dq,Dr,tau)
  pp=y(1);qq=y(2);rr=y(3);w=y(4:7)';
  pwm=[interp1(tg,pw(:,1),tt,'linear','extrap'),interp1(tg,pw(:,2),tt,'linear','extrap'),interp1(tg,pw(:,3),tt,'linear','extrap'),interp1(tg,pw(:,4),tt,'linear','extrap')];
  wss=max(0,mo.CR*pwm+mo.Omega_b); wd=(wss-w)/tau;
  Tmr=kt(:)'.*(mo.kT_rpm*w.^2); Qmr=kq(:)'.*(mo.kQ_rpm*w.^2);
  Mx=-(ar.Lx_r*Tmr(1)-ar.Lx_l*Tmr(2)-ar.Lx_l*Tmr(3)+ar.Lx_r*Tmr(4));
  My= ar.Ly_f*Tmr(1)-ar.Ly_r*Tmr(2)+ar.Ly_f*Tmr(3)-ar.Ly_r*Tmr(4);
  Mz= Qmr(1)+Qmr(2)-Qmr(3)-Qmr(4);
  pd=G(1)*pp*qq-G(2)*qq*rr+G(3)*Mx+G(4)*Mz-Dp*pp;
  qd=G(5)*pp*rr-G(6)*(pp^2-rr^2)+G(9)*My-Dq*qq;
  rd=G(7)*pp*qq-G(1)*qq*rr+G(4)*Mx+G(8)*Mz-Dr*rr;
  dy=[pd;qd;rd;wd(:)];
end
