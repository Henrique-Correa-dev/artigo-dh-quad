% tilt_check.m — R² de a_x, a_y por janela, com e sem o desalinhamento eps_y (só sensor + T da cadeia)
clear; clc; addpath(fileparts(fileparts(mfilename('fullpath')))); paths=setup_paths(); proj=parameters(); g=proj.g;
L=load_log_data(fullfile(paths.data,'logs_concat.mat')); dt=0.1;
t_lo=max([min(L.time_IMU),min(L.time_ATT),min(L.time_RCOU)]); t_hi=min([max(L.time_IMU),max(L.time_ATT),max(L.time_RCOU)]);
tg=(t_lo:dt:t_hi)'; ip=@(tt,xx) interp1(tt,xx,tg,'linear');
p=ip(L.time_IMU,L.gyrX_raw); q=ip(L.time_IMU,L.gyrY_raw); r=ip(L.time_IMU,L.gyrZ_raw);
ax=ip(L.time_IMU,L.accX_raw); ay=ip(L.time_IMU,L.accY_raw); az=ip(L.time_IMU,L.accZ_raw);
pd=gradient(p,dt); qd=gradient(q,dt); rd=gradient(r,dt);
prog=load(fullfile(paths.outputs,'P_identified.mat')); P=prog.P_final(:);
W=[ip(L.time_RCOU,L.pwm1_raw) ip(L.time_RCOU,L.pwm2_raw) ip(L.time_RCOU,L.pwm3_raw) ip(L.time_RCOU,L.pwm4_raw)];
[rpm,fT]=motor_chain(tg,W); T_m=sum(fT(rpm).*P(5:8)',2)/proj.m;
rimu=proj.imu_offset; rx=rimu(1); ry=rimu(2); rz=rimu(3);
eul_x=qd.*rz-rd.*ry; eul_y=rd.*rx-pd.*rz; cen_x=p.*q.*ry+p.*r.*rz-(q.^2+r.^2).*rx; cen_y=p.*q.*rx+q.*r.*rz-(p.^2+r.^2).*ry;
R2=@(y,yh) 1-sum((y-yh).^2)/sum((y-mean(y)).^2);
WINS=[4 125; 130 260; 550 600; 605 625; 610 630];
fprintf('%-10s | %8s %8s | %8s %8s | %8s %8s\n','janela','ax e=0','ax e=1.8','ay e=0','ay e=1.8','ax bias','ax tilt-fit');
for i=1:size(WINS,1)
  m=tg>=WINS(i,1)&tg<=WINS(i,2);
  fx0=eul_x+cen_x-0.3;                              % modelo antigo: bias -0.3
  fx1=eul_x+cen_x-T_m*sin(deg2rad(1.8))+0.0;        % novo: tilt 1,8°, bias 0
  fy0=eul_y+cen_y-0.2; fy1=fy0;                     % y não muda
  % melhor eps por janela (regressão a_x - lever = c0 - T_m*s)
  c=[ones(nnz(m),1) -T_m(m)]\(ax(m)-eul_x(m)-cen_x(m));
  fprintf('%3d–%-3d    | %8.3f %8.3f | %8.3f %8.3f | c0=%+.2f eps_fit=%.2f°\n',WINS(i,:),R2(ax(m),fx0(m)),R2(ax(m),fx1(m)),R2(ay(m),fy0(m)),R2(ay(m),fy1(m)),c(1),rad2deg(asin(c(2))));
end
