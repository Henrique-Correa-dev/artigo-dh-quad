function R = cl_loop(isNL, s, M)
%CL_LOOP  Integra a malha fechada (controlador + atuação saturada) numa planta.
%
%  Mesma lei (control_law) + mesma cadeia de atuação (forces_to_pwm + saturação
%  de PWM em M.PWM_LIM + curva real do motor) aplicada a uma das duas plantas:
%    isNL=false → LINEAR    : ẋ = A·x + B·(u_eff - u0),  ḣ = -w
%    isNL=true  → NÃO-LINEAR: corpo rígido 6-DOF (vtol_dynamics) + posição NED
%
%  s : struct do cenário  (.t, .h0, .sp_h, .sp_u, .sp_psi)
%  M : struct do modelo   (.A,.B,.u0,.P,.K,.dyn,.bridge,.fT,.fQ,.constants,
%                          .dt,.PWM_LIM)
%
%  R : struct com X(N×9) [p q r φ θ ψ u v w], H, U(N×4), PWM(N×4), Usat(N×1)

    N=numel(s.t); dt=M.dt; PWM_LIM=M.PWM_LIM;
    x=zeros(9,1); h=s.h0; pos=[0;0;-s.h0]; cs=struct('int_h',0,'int_u',0);
    X=zeros(N,9); H=zeros(N,1); U=zeros(N,4); PWM=zeros(N,4); Usat=false(N,1);
    for k=1:N
        if isNL, h=-pos(3); end
        sp=struct('h',s.sp_h(k),'u',s.sp_u(k),'psi',s.sp_psi(k));
        [u_cmd,cs]=control_law(x,h,sp,M.K,cs,dt);

        % --- alocação + SATURAÇÃO de PWM (idêntica nos dois modelos) ---
        pwm_b = M.bridge(u_cmd);
        Usat(k)= any(pwm_b>=PWM_LIM(2)-1e-6 | pwm_b<=PWM_LIM(1)+1e-6);
        pwm   = min(max(pwm_b,PWM_LIM(1)),PWM_LIM(2));
        [Tt,Mx,My,Mz]=M.dyn.forces(pwm',M.P,M.fT,M.fQ);
        u_eff=[Tt;Mx;My;Mz];

        X(k,:)=x'; H(k)=h; U(k,:)=u_eff'; PWM(k,:)=pwm';

        % --- integra a planta (RK4) ---
        if isNL
            fx=@(xx) M.dyn.rigid_body(xx,Tt,Mx,My,Mz,M.P,M.constants);
            fp=@(xx) body2ned_vel(xx);
            k1=fx(x);          p1=fp(x);
            k2=fx(x+dt/2*k1);  p2=fp(x+dt/2*k1);
            k3=fx(x+dt/2*k2);  p3=fp(x+dt/2*k2);
            k4=fx(x+dt*k3);    p4=fp(x+dt*k3);
            x   = x   + dt/6*(k1+2*k2+2*k3+k4);
            pos = pos + dt/6*(p1+2*p2+2*p3+p4);
        else
            f=@(xx,hh) deal(M.A*xx+M.B*(u_eff-M.u0), -xx(9));
            [k1x,k1h]=f(x,h);
            [k2x,k2h]=f(x+dt/2*k1x,h+dt/2*k1h);
            [k3x,k3h]=f(x+dt/2*k2x,h+dt/2*k2h);
            [k4x,k4h]=f(x+dt*k3x,  h+dt*k3h);
            x = x + dt/6*(k1x+2*k2x+2*k3x+k4x);
            h = h + dt/6*(k1h+2*k2h+2*k3h+k4h);
        end
    end
    R=struct('X',X,'H',H,'U',U,'PWM',PWM,'Usat',Usat);
end

function vned = body2ned_vel(x)
    phi=x(4); th=x(5); ps=x(6); u=x(7); v=x(8); w=x(9);
    cphi=cos(phi); sphi=sin(phi); cth=cos(th); sth=sin(th); cps=cos(ps); sps=sin(ps);
    R = [ cth*cps,  sphi*sth*cps-cphi*sps,  cphi*sth*cps+sphi*sps; ...
          cth*sps,  sphi*sth*sps+cphi*cps,  cphi*sth*sps-sphi*cps; ...
         -sth,      sphi*cth,               cphi*cth ];
    vned = R*[u;v;w];
end
