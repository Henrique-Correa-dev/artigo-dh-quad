% inertia_lumped.m — Inércia a priori: CAD (estrutura) + massas concentradas
% =========================================================================
% MOTIVO
%   O tensor de inércia usado hoje vem do SolidWorks ("modelo completo_V2"):
%       Ixx = 43,244   Iyy = 84,404   Izz = 126,192   Ixz = 1,571   [g·m²]
%   Só que aquele modelo pesa 1,072 kg, e a aeronave voa com 1,993 kg. Falta
%   0,921 kg de componentes (bateria, ESCs, eletrônica, cabos, fixações) que
%   não estão no CAD, e cuja contribuição à inércia depende de ONDE estão.
%
%   Este script recompõe o tensor pelo teorema dos eixos paralelos:
%       J = J_CAD + Σ m_i ( |r_i|² I − r_i r_iᵀ )
%   e compara com o valor identificado, para dar um a priori defensável e
%   quantificar quanto da diferença CAD × identificado a massa faltante explica.
%
% CONVENÇÃO
%   Corpo FRD: x para frente, y para a direita, z para baixo. Origem no CG do
%   CAD. Posições em metros. Massas em kg.
%
% FONTES
%   CAD          tela de propriedades de massa do SolidWorks (modelo completo_V2)
%   componentes  Tabela II de DeLucena2025 (pesos) + geometria da Figura 2
%   braços       parameters.m (Lx_* = lateral, Ly_* = longitudinal)
%
% ATENÇÃO
%   As POSIÇÕES abaixo são estimativas a partir do desenho. Corrija-as com as
%   medidas reais: a inércia é sensível ao quadrado da distância. A coluna
%   "no_cad" marca o que JÁ está no CAD (não recontabilizar).
%
% Uso:  >> inertia_lumped
% =========================================================================
clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();  p = parameters();

%% ---------------- CAD (SolidWorks) ----------------
m_cad = 1.07201;                       % kg — massa do modelo CAD
J_cad = [ 43.244  -(-1.003)  -(1.571) ; ...
          -(-1.003)  84.404  -(0.021) ; ...
          -(1.571)  -(0.021) 126.192 ] / 1000;   % kg·m², tensor com produtos NEGATIVOS
% (a tela do SolidWorks usa "notação de tensor positivo": Lxy, Lxz, Lyz são os
%  PRODUTOS de inércia, e o tensor de inércia tem esses termos com sinal trocado)
m_real = p.m;                          % kg — massa medida da aeronave (1,993)
m_falta = m_real - m_cad;

%% ---------------- COMPONENTES (nome, massa [kg], x, y, z [m], já no CAD?) ----------------
ly = mean([p.arms.Lx_r, p.arms.Lx_l]);   % 0,232 m — braço lateral
lx = mean([p.arms.Ly_f, p.arms.Ly_r]);   % 0,327 m — braço longitudinal
z_rot = -0.05;                            % rotores acima do CG (z para baixo)

%            nome              massa    x      y      z     no_cad
C = { ...
    'motor 1 (A2212+ESC+hél)', 0.092,  +lx,   +ly,  z_rot,  true ; ...
    'motor 2 (A2212+ESC+hél)', 0.092,  +lx,   -ly,  z_rot,  true ; ...
    'motor 3 (A2212+ESC+hél)', 0.092,  -lx,   -ly,  z_rot,  true ; ...
    'motor 4 (A2212+ESC+hél)', 0.092,  -lx,   +ly,  z_rot,  true ; ...
    'motor pusher D3536+hél',  0.131,  -0.50,  0.00, 0.00,  true ; ...
    'bateria LiPo 3S 2200',    0.186,  +0.05,  0.00, +0.03, false; ...
    'Pixhawk 4 + suporte',     0.030,   0.00,  0.00, +0.01, false; ...
    'GPS NEO-M8N + mastro',    0.021,  -0.10,  0.00, -0.08, false; ...
    'receptor + telemetria',   0.030,  -0.05,  0.00, +0.02, false; ...
    'ESC do pusher',           0.026,  -0.40,  0.00,  0.00, false; ...
    'placa de distribuição',   0.005,   0.00,  0.00, +0.02, false; ...
    'sensor de velocidade',    0.025,  +0.30,  0.00, -0.02, false; ...
    '3 servos + hastes',       0.015,  -0.45,  0.00,  0.00, false; ...
    'cabeamento (estimado)',   0.060,   0.00,  0.00,  0.00, false; ...
};

nome = C(:,1);  mi = cell2mat(C(:,2));  ri = cell2mat(C(:,3:5));  no_cad = cell2mat(C(:,6));
novos = ~no_cad;

fprintf('\n=============== BALANÇO DE MASSA ===============\n');
fprintf('  CAD (modelo completo_V2) : %6.3f kg\n', m_cad);
fprintf('  aeronave (medida)        : %6.3f kg\n', m_real);
fprintf('  faltando no CAD          : %6.3f kg  (%.0f%% da massa real)\n', m_falta, 100*m_falta/m_real);
fprintf('  componentes listados como fora do CAD: %6.3f kg  (%s)\n', sum(mi(novos)), ...
    ternario(abs(sum(mi(novos))-m_falta) < 0.05, 'fecha com o balanço', 'NÃO fecha — revisar a lista'));

%% ---------------- TENSOR RECOMPOSTO ----------------
Jadd = zeros(3);
for k = find(novos)'
    r = ri(k,:)';
    Jadd = Jadd + mi(k)*( (r'*r)*eye(3) - r*r' );
end
J_tot = J_cad + Jadd;

% deslocamento do CG causado pelos componentes novos (deveria ser pequeno)
r_cg = sum(mi(novos).*ri(novos,:), 1)' / m_real;

fprintf('\n=============== TENSOR ===============\n');
fprintf('  %-14s %8s %8s %8s %9s\n', '', 'Jx', 'Jy', 'Jz', 'Jxz');
fprintf('  %-14s %8.4f %8.4f %8.4f %9.5f\n', 'CAD', J_cad(1,1), J_cad(2,2), J_cad(3,3), -J_cad(1,3));
fprintf('  %-14s %8.4f %8.4f %8.4f %9.5f\n', '+ componentes', J_tot(1,1), J_tot(2,2), J_tot(3,3), -J_tot(1,3));
lmv = load(fullfile(paths.outputs,'linear_model.mat'));  Pid = lmv.P(:);
fprintf('  %-14s %8.4f %8.4f %8.4f %9.5f   (Jz, Jxz travados no CAD)\n', 'identificado', Pid(1), Pid(2), Pid(3), Pid(4));
fprintf('  %-14s %8.4f %8.4f %8.4f %9.5f   (Tabela 5.3)\n', 'dissertação', 0.0469, 0.0984, 0.1453, 0.00261);
fprintf('\n  desvio CG pelos componentes novos: [%+.3f %+.3f %+.3f] m\n', r_cg);
fprintf('  folga triangular Jx+Jy-Jz: CAD %+.5f | + componentes %+.5f\n', ...
    J_cad(1,1)+J_cad(2,2)-J_cad(3,3), J_tot(1,1)+J_tot(2,2)-J_tot(3,3));

%% ---------------- QUE RAIO EXPLICA O IDENTIFICADO? ----------------
% Se a massa faltante estivesse concentrada a um raio efetivo k de cada eixo,
% J_id = J_CAD + m_falta·k² → k = sqrt((J_id − J_CAD)/m_falta).
fprintf('\n=============== RAIO EFETIVO IMPLICADO ===============\n');
eixo = {'x (rolagem)','y (arfagem)'};
Jc = [J_cad(1,1), J_cad(2,2)];  Ji = [Pid(1), Pid(2)];
for i = 1:2
    dJ = Ji(i) - Jc(i);
    k  = sqrt(max(dJ,0)/m_falta);
    fprintf('  %-12s  CAD %6.4f → id %6.4f  (%+.1f%%)  ⇒ raio efetivo %.3f m\n', ...
        eixo{i}, Jc(i), Ji(i), 100*dJ/Jc(i), k);
end
fprintf('  Leitura: raio de 5 a 8 cm é compatível com bateria e eletrônica dentro\n');
fprintf('  da fuselagem, perto do CG. Raio grande (>20 cm) indicaria massa nas asas.\n');

%% ---------------- CONTRIBUIÇÃO DE CADA COMPONENTE ----------------
fprintf('\n=============== CONTRIBUIÇÃO POR COMPONENTE (só os fora do CAD) ===============\n');
fprintf('  %-26s %7s %8s %8s %8s\n', 'componente', 'massa', 'dJx', 'dJy', 'dJz');
for k = find(novos)'
    r = ri(k,:)';  dJ = mi(k)*( (r'*r)*eye(3) - r*r' );
    fprintf('  %-26s %7.3f %8.5f %8.5f %8.5f\n', nome{k}, mi(k), dJ(1,1), dJ(2,2), dJ(3,3));
end

fprintf('\n  Sugestão de a priori (P0_J), se o balanço fechar:\n');
fprintf('    p.J.Jx = %.3f / 1000;  p.J.Jy = %.3f / 1000;  p.J.Jz = %.3f / 1000;\n', ...
    1000*J_tot(1,1), 1000*J_tot(2,2), 1000*J_tot(3,3));

function s = ternario(c, a, b), if c, s = a; else, s = b; end, end
