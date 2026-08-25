function VE = set_aero_vel(L, tg)
%SET_AERO_VEL  Publica a velocidade estimada de bordo para os termos ∝ V.
%
%   VE = set_aero_vel(L, tg)
%
% Os momentos aerodinâmicos da estrutura, P(17:22), são escalados por V(t), que
% é EXÓGENA: vem de estimate_velocity (EKF vertical com barômetro mais
% integração ancorada da aceleração quase-estática). Ela é publicada por appdata
% porque o vtol_dynamics é chamado pelo ode45 e não tem por onde recebê-la.
%
% Sem esta chamada o vtol_dynamics DESLIGA a aerodinâmica e emite aviso. Nunca
% use os estados u, v, w para isso: eles divergem em malha aberta.
    VE = estimate_velocity(L, tg(:));
    V = VE.V;  v = VE.v;  w = VE.w;  u = VE.u;
    V(~isfinite(V)) = 0;  v(~isfinite(v)) = 0;  w(~isfinite(w)) = 0;  u(~isfinite(u)) = 0;
    setappdata(0, 'aero_vel', struct('t', tg(:), 'V', V, 'v', v, 'w', w, 'u', u));
end
