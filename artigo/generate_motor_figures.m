%GENERATE_MOTOR_FIGURES  Figuras do modelo do motor para o artigo (Secao III).
%
% Gera, a partir da tabela de bancada em parameters.m (fonte unica):
%   images/motor_pwm_rpm.png       - PWM normalizado (sigma) x Omega + ajuste linear (C_R, Omega_b)
%   images/motor_thrust_torque.png - Omega x Empuxo e Omega x Contra-torque + ajustes quadraticos
%
% Estilo academico com cor: fonte Times, circulos azuis para dados,
% linha continua vermelho-escuro para os ajustes, anotacoes em preto.
%
% Mesmas contas do modelo fisico de parameters.m: ajuste linear sem o ponto
% de motor parado; kT/kQ por minimos quadrados ponderados por Omega^2.

mdl = fullfile(fileparts(mfilename('fullpath')), '..', 'identification', '2_model');
img = fullfile(fileparts(mfilename('fullpath')), 'images');
addpath(mdl); p = parameters();

pwm = p.bench.pwm; rpm = p.bench.RPM;
TN  = p.bench.T_grams*9.80665/1000; Q = p.bench.Q_Nm;
sig = (pwm-1000)/1000;

% Ajuste linear sigma->RPM (regime permanente), sem o ponto de motor parado
cf = polyfit(sig(2:end), rpm(2:end), 1); CR = cf(1); Ob = cf(2);
kT = p.motor.kT_rpm; kQ = p.motor.kQ_rpm;
fprintf('CR_norm = %.1f RPM | Omega_b = %.1f RPM | kT = %.3e N/RPM^2 | kQ = %.3e Nm/RPM^2\n', CR, Ob, kT, kQ);

fonte = 'Times New Roman';
azul  = [0 0.447 0.741];
verm  = [0.75 0.15 0.10];

estilo = @(ax) set(ax, 'FontName',fonte, 'FontSize',9, 'Color','w', ...
    'XColor','k', 'YColor','k', 'LineWidth',0.6, ...
    'GridColor',[0.7 0.7 0.7], 'GridAlpha',0.5, 'GridLineStyle',':', ...
    'TickDir','in', 'Layer','top');

% ================= Figura A: PWM normalizado x Omega =================
f1 = figure('Visible','off','Units','inches','Position',[1 1 3.5 2.5],'Color','w');
try, f1.Theme = 'light'; catch, end
ax = axes(f1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
ss = linspace(0,1,50);
hF = plot(ax, ss, CR*ss+Ob, '-', 'Color',verm,'LineWidth',1.2);
hB = plot(ax, sig(2:end), rpm(2:end), 'o','MarkerSize',5,'MarkerFaceColor',azul,'MarkerEdgeColor',azul,'LineStyle','none');
xlabel(ax,'\sigma (PWM normalizado)','FontName',fonte); ylabel(ax,'\Omega [RPM]','FontName',fonte);
estilo(ax); xlim(ax,[0 1]); ylim(ax,[0 7000]);
text(ax,0.05,6400,sprintf('\\Omega_{ss} = %.0f\\sigma + %.0f [RPM]',CR,Ob), ...
    'FontSize',9,'Color','k','FontName',fonte);
lg = legend(ax,[hB hF],{'Medições de bancada','Ajuste linear'}, ...
    'Location','southeast','FontSize',8,'TextColor','k','Color','w','EdgeColor','k','FontName',fonte);
lg.ItemTokenSize = [16 18];
exportgraphics(f1, fullfile(img,'motor_pwm_rpm.png'), 'Resolution',300, 'BackgroundColor','white');
close(f1);

% ============ Figura B: Omega x Empuxo e Omega x Contra-torque ============
f2 = figure('Visible','off','Units','inches','Position',[1 1 3.5 4.2],'Color','w');
try, f2.Theme = 'light'; catch, end
tl = tiledlayout(f2,2,1,'TileSpacing','compact','Padding','compact');
ww = linspace(0,6500,200);

ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
hF1 = plot(ax1, ww, kT*ww.^2, '-','Color',verm,'LineWidth',1.2);
hB1 = plot(ax1, rpm, TN, 'o','MarkerSize',5,'MarkerFaceColor',azul,'MarkerEdgeColor',azul,'LineStyle','none');
ylabel(ax1,'T [N]','FontName',fonte);
estilo(ax1); xlim(ax1,[0 6500]); ylim(ax1,[0 8]);
text(ax1,300,7.40,'(a)','FontSize',9,'Color','k','FontName',fonte);
text(ax1,300,6.50,'T = c_T\Omega^2','FontSize',9,'Color','k','FontName',fonte);
text(ax1,300,5.60,'c_T = 1,76\times10^{-7} N/RPM^2','FontSize',8,'Color','k','FontName',fonte);
lg1 = legend(ax1,[hB1 hF1],{'Bancada','Ajuste quadrático'}, ...
    'Location','southeast','FontSize',8,'TextColor','k','Color','w','EdgeColor','k','FontName',fonte);
lg1.ItemTokenSize = [16 18];

ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
hF2 = plot(ax2, ww, kQ*ww.^2, '-','Color',verm,'LineWidth',1.2);
hB2 = plot(ax2, rpm, Q, 'o','MarkerSize',5,'MarkerFaceColor',azul,'MarkerEdgeColor',azul,'LineStyle','none');
ylabel(ax2,'Q [N\cdotm]','FontName',fonte); xlabel(ax2,'\Omega [RPM]','FontName',fonte);
estilo(ax2); xlim(ax2,[0 6500]); ylim(ax2,[0 0.21]);
text(ax2,300,0.194,'(b)','FontSize',9,'Color','k','FontName',fonte);
text(ax2,300,0.170,'Q = c_M\Omega^2','FontSize',9,'Color','k','FontName',fonte);
text(ax2,300,0.146,'c_M = 5,04\times10^{-9} N\cdotm/RPM^2','FontSize',8,'Color','k','FontName',fonte);
lg2 = legend(ax2,[hB2 hF2],{'Bancada','Ajuste quadrático'}, ...
    'Location','southeast','FontSize',8,'TextColor','k','Color','w','EdgeColor','k','FontName',fonte);
lg2.ItemTokenSize = [16 18];

exportgraphics(f2, fullfile(img,'motor_thrust_torque.png'), 'Resolution',300, 'BackgroundColor','white');
close(f2);
disp('Figuras geradas em images/.');
