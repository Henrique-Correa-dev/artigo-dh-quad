function Q = damp_convert(P, to_form)
%DAMP_CONVERT  Converte o bloco de amortecimento P(13:15) entre as duas formas.
%
%   Q = damp_convert(P, 'moment')   c_p,c_q,c_r [1/s]  →  L_p,M_q,N_r [N·m·s]
%   Q = damp_convert(P, 'rate')     L_p,M_q,N_r        →  c_p,c_q,c_r
%
% A equivalência é EXATA no próprio eixo e usa as inércias do próprio vetor P:
%     c_p = Γ3·L_p        c_q = M_q/J_y        c_r = Γ8·N_r
% com Γ3 = J_z/(J_x·J_z − J_xz²) e Γ8 = J_x/(J_x·J_z − J_xz²).
%
% ATENÇÃO: a conversão NÃO é uma equivalência completa do modelo. A forma de
% momento acrescenta os dois termos cruzados via J_xz (−Γ4·N_r·r em ṗ e
% −Γ4·L_p·p em ṙ) que a forma de taxa não tem. Converter e simular dá
% resultados ligeiramente diferentes, e essa diferença é justamente o efeito
% que se quer medir.
    Jx = P(1); Jy = P(2); Jz = P(3); Jxz = P(4);
    g0 = Jx*Jz - Jxz^2;  G3 = Jz/g0;  G8 = Jx/g0;
    Q = P(:);
    switch lower(to_form)
        case 'moment', Q(13) = P(13)/G3;  Q(14) = P(14)*Jy;  Q(15) = P(15)/G8;
        case 'rate',   Q(13) = P(13)*G3;  Q(14) = P(14)/Jy;  Q(15) = P(15)*G8;
        otherwise, error('damp_convert:form','use ''moment'' ou ''rate''');
    end
end
