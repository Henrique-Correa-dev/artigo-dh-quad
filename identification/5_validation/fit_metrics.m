function m = fit_metrics(y, yhat)
%FIT_METRICS  Métricas de ajuste medido vs simulado para validação de sysid.
%
%   m = fit_metrics(y, yhat)
%
% Calcula, entre o sinal MEDIDO y e o SIMULADO yhat:
%   .R2     coeficiente de determinação (1=perfeito, 0=média, <0=ruim)
%   .RMSE   raiz do erro quadrático médio (UNIDADES FÍSICAS do sinal)
%   .NRMSE  RMSE normalizado pelo desvio-padrão do medido
%   .TIC    Theil Inequality Coefficient ∈ [0,1] (0=perfeito)
%           Critério Klein & Morelli: TIC < 0.25–0.30 = validação aceitável.
%   .U_bias .U_var .U_cov   decomposição de Theil (somam 1):
%           U_bias = erro de MÉDIA (offset/viés)        → quer ~0
%           U_var  = erro de AMPLITUDE/escala            → quer ~0
%           U_cov  = parte ALEATÓRIA (não-modelável)     → quer ~1
%           Se U_bias ou U_var grandes → há erro SISTEMÁTICO a corrigir.
%
% Sem dependência de toolbox (correlação calculada à mão).

    y = y(:);  yhat = yhat(:);
    e = y - yhat;
    mse = mean(e.^2);

    % R²
    ss_tot = sum((y - mean(y)).^2);
    m.R2 = 1 - sum(e.^2) / max(ss_tot, 1e-12);

    % RMSE / NRMSE
    m.RMSE  = sqrt(mse);
    m.NRMSE = m.RMSE / max(std(y), 1e-12);

    % Theil Inequality Coefficient
    rms_y  = sqrt(mean(y.^2));
    rms_yh = sqrt(mean(yhat.^2));
    m.TIC = sqrt(mse) / max(rms_y + rms_yh, 1e-12);

    % Decomposição de Theil: MSE = (μy-μŷ)² + (σy-σŷ)² + 2(1-ρ)σy σŷ
    my = mean(y);  myh = mean(yhat);
    sy  = sqrt(mean((y - my).^2));         % desvio populacional
    syh = sqrt(mean((yhat - myh).^2));
    cov_yyh = mean((y - my) .* (yhat - myh));
    rho = cov_yyh / max(sy * syh, 1e-12);
    if ~isfinite(rho), rho = 0; end

    m.U_bias = (my - myh)^2          / max(mse, 1e-12);
    m.U_var  = (sy - syh)^2          / max(mse, 1e-12);
    m.U_cov  = 2*(1 - rho)*sy*syh    / max(mse, 1e-12);
end
