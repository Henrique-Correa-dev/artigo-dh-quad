function K = kdm_of(P, time, m)
%KDM_OF  Amortecimento translacional por amostra [1/s], colunas x, y, z.
%
%   K = kdm_of(P, time, m)   devolve N×3
%
% Reúne as duas fontes que agem em u, v, w:
%   x, y : arrasto induzido de rotor, g·C_d com C_d = P(16), constante (Beard 14.4.2)
%          MAIS a força aerodinâmica da estrutura ½ρS·C_Xu·V/m e ½ρS·C_Yv·V/m
%   z    : só a estrutura, ½ρS·C_Zw·V/m  (o arrasto de rotor age em x,y)
%
% A velocidade V vem do appdata 'aero_vel' (estimate_velocity), interpolada em
% `time`. Sem appdata, ou com P de menos de 25 elementos, sobra só o C_d.
    N = numel(time);  K = zeros(N,3);
    if numel(P) >= 16
        pp = parameters();
        K = repmat(pp.g*P(16)*[1 1 0], N, 1);
    end
    if numel(P) >= 25 && any(P(23:25) ~= 0)
        av = getappdata(0,'aero_vel');
        if ~isempty(av) && isstruct(av)
            V = interp1(av.t, av.V, time(:), 'linear', 0);
            V(~isfinite(V)) = 0;
            kF = aero_gains();
            K = K + (kF.F/m) * (V * P(23:25)');
        end
    end
end
