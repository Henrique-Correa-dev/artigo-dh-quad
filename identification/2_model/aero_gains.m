function k = aero_gains(prm)
%AERO_GAINS  Ganhos dimensionais dos momentos aerodinâmicos da estrutura.
%
%   k = aero_gains()      lê a geometria de parameters()
%   k = aero_gains(prm)   usa um struct de parameters() já carregado
%
%   MODELO (Salahudden et al. 2023, Eq. 13, com η_aero→hover = 1, e as razões
%   adimensionais abertas em β = v/V, α = w/V, p̄ = p·b/(2V), q̄ = q·c̄/(2V),
%   r̄ = r·b/(2V)), com q̄∞ = ½ρV²:
%
%     L_aero = q̄∞·S·b·(Cl_β·β + Cl_p·p̄) = ¼ρSb²·Cl_p·V·p + ½ρSb·Cl_β·V·v
%     M_aero = q̄∞·S·c̄·(Cm_α·α + Cm_q·q̄) = ¼ρSc̄²·Cm_q·V·q + ½ρSc̄·Cm_α·V·w
%     N_aero = q̄∞·S·b·(Cn_β·β + Cn_r·r̄) = ¼ρSb²·Cn_r·V·r + ½ρSb·Cn_β·V·v
%
%   O ponto central: V aparece na PRIMEIRA potência. A forma é REGULAR em
%   V → 0 (todos os termos valem zero no pairado), sem a singularidade 1/V da
%   forma adimensional. É o que permite usar a mesma equação no pairado e em
%   voo à frente, sem função de blend.
%
%   Uso típico (todos escalares ou vetores da mesma altura):
%     k = aero_gains();
%     L_a = k.Lp*Cl_p*V.*p + k.Lb*Cl_b*V.*v;
%     M_a = k.Mq*Cm_q*V.*q + k.Ma*Cm_a*V.*w;
%     N_a = k.Nr*Cn_r*V.*r + k.Nb*Cn_b*V.*v;
%
%   CONVENÇÃO dos adimensionais: a mesma da AVL e do XFLR5 (pb/2V, qc̄/2V,
%   rb/2V), portanto os chutes iniciais de dh_st.txt entram direto.
%   ATENÇÃO (correção): os scripts de diagnóstico anteriores (make_diag.py,
%   prior_damping.m) usavam ½ρSc̄²·Cm_q no arfagem, isto é, normalização por
%   q·c̄/V, 2× maior que o correto. Aqui está ¼, coerente com a AVL.

    %  Cache: a geometria não muda em tempo de execução e esta função é chamada
    %  dentro do ode45. Depois de editar parameters(), rode  clear aero_gains.
    persistent K
    if nargin < 1 && ~isempty(K), k = K; return; end
    if nargin < 1 || isempty(prm), prm = parameters(); end
    rho = prm.rho;  S = prm.wing.S;  b = prm.wing.b;  c = prm.wing.c;

    k = struct( ...
        'Lp', 0.25*rho*S*b^2, ...   % N·m por (Cl_p · V · p)
        'Lb', 0.50*rho*S*b,   ...   % N·m por (Cl_β · V · v)
        'Mq', 0.25*rho*S*c^2, ...   % N·m por (Cm_q · V · q)
        'Ma', 0.50*rho*S*c,   ...   % N·m por (Cm_α · V · w)
        'Nr', 0.25*rho*S*b^2, ...   % N·m por (Cn_r · V · r)
        'Nb', 0.50*rho*S*b,   ...   % N·m por (Cn_β · V · v)
        'F',  0.50*rho*S,     ...   % N por (C_Xu·V·u), (C_Yv·V·v), (C_Zw·V·w)
        'rho', rho, 'S', S, 'b', b, 'c', c);
    K = k;
end
