function align_yyaxis_zero(ax)
%ALIGN_YYAXIS_ZERO  Alinha o zero dos dois eixos Y de um subplot com yyaxis.
%
% Garante que "0" fica na mesma posição visual em ambos os eixos.
% Estende o lado negativo do eixo com fração menor (nunca corta dado).
%
% USO:
%   yyaxis left;  plot(...);
%   yyaxis right; plot(...);
%   align_yyaxis_zero(gca);

    yyaxis(ax, 'left');  yl = ylim(ax);
    yyaxis(ax, 'right'); yr = ylim(ax);

    % Garante que zero está dentro de cada faixa
    yl(1) = min(yl(1), 0); yl(2) = max(yl(2), 0);
    yr(1) = min(yr(1), 0); yr(2) = max(yr(2), 0);

    % Fração da altura que está abaixo de zero
    f_L = -yl(1) / max(yl(2) - yl(1), eps);
    f_R = -yr(1) / max(yr(2) - yr(1), eps);
    f = max(f_L, f_R);

    if f <= 1e-6 || f >= 1 - 1e-6
        return;  % caso degenerado (zero está numa extremidade)
    end

    % Estende lado negativo do eixo com fração menor (preserva yl(2)/yr(2))
    if f_L < f, yl(1) = -f * yl(2) / (1 - f); end
    if f_R < f, yr(1) = -f * yr(2) / (1 - f); end

    yyaxis(ax, 'left');  ylim(ax, yl);
    yyaxis(ax, 'right'); ylim(ax, yr);
end
