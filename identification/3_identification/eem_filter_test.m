% eem_filter_test.m — o valor de c_p do EEM depende do filtro da derivada?
% =========================================================================
% No identify_plant o EEM monta o resíduo com ṗ = d/dt(movmean(p,5)) mas usa o
% p CRU no termo de amortecimento. Como a média móvel de 0,5 s tem ganho 0,3 em
% 1,5 Hz, o lado esquerdo fica atenuado e o direito não, o que puxa TODOS os
% parâmetros de escala (k_T e c juntos) para baixo. Este teste refaz o mesmo
% ajuste com três tratamentos do filtro para medir o tamanho do viés:
%   (a) assimétrico  — como está hoje
%   (b) simétrico    — mesmo filtro nos dois lados (correto)
%   (c) sem filtro   — derivada crua
%
% Uso:  >> eem_filter_test
% =========================================================================
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  proj = parameters();
T_TRAINS = {[4,24];[25,41];[42,62];[63,99];[100,125]};
dt = 0.1;  SW = 5;

L = load_log_data(fullfile(paths.data,'logs_concat.mat'));
t_lo = max([min(L.time_IMU), min(L.time_ATT), min(L.time_RCOU)]);
t_hi = min([max(L.time_IMU), max(L.time_ATT), max(L.time_RCOU)]);
tg = (t_lo:dt:t_hi)';
ip = @(tt,xx) interp1(tt, xx, tg, 'linear');
W4 = [ip(L.time_RCOU,L.pwm1_raw), ip(L.time_RCOU,L.pwm2_raw), ip(L.time_RCOU,L.pwm3_raw), ip(L.time_RCOU,L.pwm4_raw)];
[rpm, fT, fQ] = motor_chain(tg, W4);
G = [ip(L.time_IMU,L.gyrX_raw), ip(L.time_IMU,L.gyrY_raw), ip(L.time_IMU,L.gyrZ_raw)];

MODOS = {'assimétrico (atual)','simétrico','sem filtro'};
P0 = proj.P0_J(:);  lb = proj.bounds.lb(:); ub = proj.bounds.ub(:);
lb(16:22) = P0(16:22); ub(16:22) = P0(16:22);      % aerodinâmica travada (só o bloco rotacional)
RL = struct('kt_pair',100,'kq_pair',100,'tri',1e4);
fprintf('\n  %-22s %8s %8s %8s %8s %8s\n','filtro','k_T méd','c_p','c_q','c_r','Jx');
for md = 1:3
    pa=[]; qa=[]; ra=[]; pd=[]; qd=[]; rd=[]; Tr=[]; Qr=[];
    for s = 1:numel(T_TRAINS)
        ii = tg>=T_TRAINS{s}(1) & tg<=T_TRAINS{s}(2);
        ps = G(ii,1); qs = G(ii,2); rs = G(ii,3);
        switch md
            case 1, pf=ps; qf=qs; rf=rs;  ds=@(x) gradient(movmean(x,SW),dt);
            case 2, pf=movmean(ps,SW); qf=movmean(qs,SW); rf=movmean(rs,SW); ds=@(x) gradient(movmean(x,SW),dt);
            case 3, pf=ps; qf=qs; rf=rs;  ds=@(x) gradient(x,dt);
        end
        pa=[pa;pf]; qa=[qa;qf]; ra=[ra;rf];
        pd=[pd;ds(ps)]; qd=[qd;ds(qs)]; rd=[rd;ds(rs)];
        Tr=[Tr; fT(rpm(ii,:))]; Qr=[Qr; fQ(rpm(ii,:))];
    end
    wts = [1/var(pd); 1/var(qd); 1/var(rd)];
    cost = @(P) eem_cost_function(ones(22,1), P, wts, pa,qa,ra, pd,qd,rd, Tr,Qr, RL);
    o = optimoptions('lsqnonlin','Display','off','MaxIterations',800,'MaxFunctionEvaluations',4e4, ...
        'FunctionTolerance',1e-14,'StepTolerance',1e-14);
    P = lsqnonlin(cost, P0, lb, ub, o);
    fprintf('  %-22s %8.3f %8.3f %8.3f %8.3f %8.4f\n', MODOS{md}, mean(P(5:8)), P(13), P(14), P(15), P(1));
end
fprintf('\n  (OEM oficial: k_T 1.06 | c_p 5.50 | c_q 3.86 | c_r 0.78 | Jx 0.0475)\n');
