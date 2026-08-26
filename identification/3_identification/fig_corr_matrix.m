load('/Users/graest/ita-master/artigo/artigo-dh-quad/identification/outputs/struct_ident.mat');
labels={'J_x','J_y','J_z','J_{xz}','c_{T1}','c_{T2}','c_{T3}','c_{T4}','c_{M1}','c_{M2}','c_{M3}','c_{M4}','L_p','M_q','N_r'};
R=OUT.full.Rho; n=15;
fig=figure('Position',[80 80 900 780],'Color','w'); try, theme(fig,'light'); catch, end
imagesc(R,[-1 1]); axis square;
cmap=[linspace(0.75,1,128)' linspace(0.30,1,128)' linspace(0.20,1,128)'; ...
      linspace(1,0.10,128)' linspace(1,0.35,128)' linspace(1,0.65,128)'];
colormap(cmap); colorbar;
set(gca,'XTick',1:n,'XTickLabel',labels,'YTick',1:n,'YTickLabel',labels,'FontSize',11,'TickLabelInterpreter','tex');
for i=1:n, for j=1:n
    v=R(i,j);
    if i==j, txt='1'; else, txt=sprintf('%.2f',v); end
    col='k'; if abs(v)>0.75, col='w'; end
    text(j,i,txt,'HorizontalAlignment','center','FontSize',8,'Color',col);
end, end
exportgraphics(fig,'/Users/graest/ita-master/artigo/artigo-dh-quad/dissertacao/Cap5/corr_matrix.png','Resolution',180);
fprintf('salvo\n');
