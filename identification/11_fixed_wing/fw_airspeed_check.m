% fw_airspeed_check.m — o pitot está calibrado?
% =========================================================================
% Os coeficientes identificados nos doublets saem MENORES que os da AVL por um
% fator quase constante dentro de cada eixo. Isso é a assinatura de um erro de
% ESCALA comum, e o candidato mais simples é a pressão dinâmica q̄ = ½ρV², ou
% seja, o pitot. O parâmetro ARSPD_RATIO estava em 2 (valor default, sem
% calibração em voo) e ARSPD_AUTOCAL desligado.
%
% Teste: em voo nivelado com vento constante, a velocidade em relação ao solo
% do EKF gira com a proa enquanto a velocidade do ar não. Ajustando
%     V_solo_NED(t) = V_ar·[cosψ; senψ; 0] + W
% por mínimos quadrados num trecho com proas variadas, saem de uma vez a
% velocidade do ar verdadeira e o vento. Comparando com o pitot sai o fator.
%
% Uso:  >> fw_airspeed_check
% =========================================================================
clear; clc; close all;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

% RESULTADO DESTE TESTE: os dois voos de dez/2024 foram feitos SEM GPS. A
% mensagem GPS não existe no log e o EKF não estimou nem velocidade em relação
% ao solo (XKF1.VN/VE = 0) nem vento (XKF2.VWN/VWE = 0). Logo o triângulo do
% vento não fecha e a calibração do pitot NÃO PODE ser verificada com este dado.
% O script fica aqui documentando a tentativa e detecta o caso.
for which = {'long','lat'}
    F = fw_load(which{1});
    if max(abs(F.VN)) < 1e-6 && max(abs(F.VE)) < 1e-6
        fprintf(['\n  === %s: o EKF não tem velocidade em relação ao solo ' ...
                 '(voo sem GPS).\n      Pitot médio %.2f m/s. Sem referência ' ...
                 'independente, o fator de escala do pitot fica em aberto.\n'], ...
                 which{1}, mean(F.V(F.V>12)));
        continue;
    end
    ok = F.V > 12 & isfinite(F.VN) & isfinite(F.VE);
    t = F.t(ok);  psi = F.psi(ok);  VN = F.VN(ok);  VE = F.VE(ok);  Vp = F.V(ok);
    Vg = hypot(VN, VE);
    fprintf('\n  === %s: %d amostras em voo (%.0f a %.0f s)\n', which{1}, numel(t), t(1), t(end));
    fprintf('  pitot: média %.2f  |  solo (EKF): média %.2f  |  razão %.3f\n', ...
        mean(Vp), mean(Vg), mean(Vg)/mean(Vp));
    fprintf('  proa varrida: %.0f°\n', rad2deg(max(psi)-min(psi)));

    % ajuste do triângulo do vento: [cosψ senψ 1 0; -senψ cosψ 0 1]·[Va; Wn; We]
    A = [cos(psi), ones(size(psi)), zeros(size(psi)); ...
         sin(psi), zeros(size(psi)), ones(size(psi))];
    bvec = [VN; VE];
    x = A\bvec;                      % [Va, WN, WE]
    res = bvec - A*x;
    fprintf('  ajuste do triângulo do vento: V_ar = %.2f m/s | vento (%.2f, %.2f) m/s = %.2f m/s de %.0f°\n', ...
        x(1), x(2), x(3), hypot(x(2),x(3)), mod(rad2deg(atan2(x(3),x(2)))+180,360));
    fprintf('  resíduo RMS %.2f m/s   |   FATOR pitot: V_ar/V_pitot = %.3f  → q̄ escala por %.3f\n', ...
        sqrt(mean(res.^2)), x(1)/mean(Vp), (x(1)/mean(Vp))^2);

    if strcmp(which{1},'long'), Fl = F; xl = x; end
end
if ~exist('xl','var')
    fprintf('\n  Nenhum dos voos tem velocidade de solo: nada a plotar.\n');
    return;
end

%% figura: velocidade de solo contra proa, com o ajuste
F = Fl;  x = xl;
ok = F.V > 12 & isfinite(F.VN);
psi = F.psi(ok);  Vg = hypot(F.VN(ok), F.VE(ok));  Vp = F.V(ok);
psi_g = linspace(-pi, pi, 200)';
Vfit = hypot(x(1)*cos(psi_g) + x(2), x(1)*sin(psi_g) + x(3));

f = figure('Position',[40 40 1150 450],'Color','w'); try, f.Theme='light'; catch, end
subplot(1,2,1); hold on; grid on;
scatter(rad2deg(psi), Vg, 8, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerFaceColor',[0 0.45 0.7]);
plot(rad2deg(psi_g), Vfit, 'r-', 'LineWidth',2);
yline(mean(Vp), 'k--', 'média do pitot', 'LineWidth',1.3);
yline(x(1), '-', sprintf('V_{ar} ajustada = %.1f m/s', x(1)), 'Color',[0.85 0.37 0.01], 'LineWidth',1.5);
xlabel('proa \psi [°]'); ylabel('velocidade em relação ao solo [m/s]');
title('triângulo do vento (voo longitudinal)');
legend({'EKF','ajuste','pitot'}, 'Location','south');

subplot(1,2,2); hold on; grid on;
plot(F.t(ok), Vp, 'k-', 'LineWidth',1); plot(F.t(ok), Vg, '-', 'Color',[0 0.45 0.7], 'LineWidth',1);
yline(x(1), '-', 'V_{ar} ajustada', 'Color',[0.85 0.37 0.01], 'LineWidth',1.5);
xlabel('t [s]'); ylabel('[m/s]'); legend({'pitot','solo (EKF)'},'Location','best');
title('série temporal');
fn = fullfile(paths.images,'new_flights','fw_airspeed_check.png');
exportgraphics(f, fn, 'BackgroundColor','white','Resolution',130);
fprintf('\n  Figura: %s\n', fn);
