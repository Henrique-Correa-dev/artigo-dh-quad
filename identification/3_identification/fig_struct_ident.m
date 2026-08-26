% fig_ident.m — figura da analise de identificabilidade estrutural local
load('/Users/graest/ita-master/artigo/artigo-dh-quad/identification/outputs/struct_ident.mat');
labels = {'J_x','J_y','J_z','J_{xz}','c_{T1}','c_{T2}','c_{T3}','c_{T4}','c_{M1}','c_{M2}','c_{M3}','c_{M4}','L_p','M_q','N_r'};
svF = OUT.full.sv/OUT.full.sv(1);   svR = OUT.rotational.sv/OUT.rotational.sv(1);
sdF = OUT.full.relsd;               sdR = OUT.rotational.relsd;

fig = figure("Position",[100 100 820 430],"Color","w"); try, theme(fig,"light"); catch, end

b = bar(1:15, [sdF(:) sdR(:)], 'grouped');
set(gca,'YScale','log'); grid on; ylim([1e-2 12]);
set(gca,'XTick',1:15,'XTickLabel',labels,'FontSize',12);
ylabel('incerteza mínima relativa de cada parâmetro');
legend('observação completa (giroscópios e acelerômetro)','observação restrita aos giroscópios','Location','northeast','FontSize',11);

out = '/Users/graest/ita-master/artigo/artigo-dh-quad/dissertacao/Cap3/identificabilidade_local.png';
exportgraphics(fig, out, 'Resolution', 200);
fprintf('salvo %s\n', out);
