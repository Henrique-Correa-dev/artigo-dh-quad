function p = fev2024(p)
%FEV2024  Campanha de 17/02/2024 (voo "quarto teste", log VooVert_001).
%
%   Aplicada por parameters('campaign','fev2024'). Recebe o struct com os
%   valores default (campanha de maio de 2026, a da dissertação) e devolve só
%   o que muda. NÃO reidentifique nada aqui: este arquivo guarda o que se sabe
%   da aeronave NAQUELE dia.
%
% O QUE MUDA, e como foi obtido
% -----------------------------------------------------------------------
% [1] POSIÇÃO DO CG.  É a diferença dominante. Em voo pairado a soma dos
%     momentos de empuxo em torno do CG tem de ser zero, logo a assimetria
%     média de empuxo mede diretamente onde o CG está. Com os empuxos médios
%     da janela 141 a 211 s (10_new_flights/check_motor_map.m):
%        empuxos médios  [4,30  6,06  3,11  4,04] N   (total 17,5 N)
%        momento residual no CG do CAD:  Mx +0,273 N·m   My −0,796 N·m
%        → CG 45,4 mm ATRÁS e 15,6 mm à ESQUERDA do ponto do CAD
%     Em maio de 2026 esses mesmos números dão 1,4 mm e 1,2 mm, ou seja, a
%     aeronave estava praticamente balanceada. Sem esta correção o modelo
%     carrega um momento de arfagem constante de 0,8 N·m, que dividido por Jy
%     vale 9,4 rad/s²: em 3 s de simulação em malha aberta q dispara, e é por
%     isso que a identificação nesta campanha dava R² próximo de zero.
% [2] MOTORES.  São outros (os de 2026 são os SunnySky A2212-15 de bancada).
%     Não há ensaio de bancada dos antigos, então a curva de empuxo continua
%     sendo a de 2026 e o k_T da identificação absorve a diferença de escala.
%     Espera-se k_T em torno de 1,1 a 1,3 (no PWM voado a curva de 2026 dá
%     17,5 N contra 19,6 N de peso). Se aparecer bancada dos motores antigos,
%     é aqui que a tabela entra.
% [3] MASSA.  Não foi pesada nesse dia. Mantida em 1,993 kg (valor medido em
%     2026). Massa e escala de empuxo são degeneradas no canal de aceleração
%     vertical, então o que sobrar de erro vai para o k_T.
% [4] SATURAÇÃO.  MOT_SPIN_MAX = 0,95 com MOT_PWM_MIN/MAX 1000/2000 limita a
%     saída em 1950 µs. O motor 2 voou com média 1888 µs e bateu 1949 µs em
%     boa parte do voo, ou seja, esteve saturado. Nesses trechos o empuxo real
%     não acompanha mais o comando e o dado não serve para identificar.
% -----------------------------------------------------------------------

    % [1] CG estimado pelo equilíbrio (dx > 0 é CG à frente, dy > 0 é CG à direita)
    dx = -0.0454;      % m — CG 45,4 mm atrás do ponto de referência do CAD
    dy = -0.0156;      % m — CG 15,6 mm à esquerda
    Lx_nom = 0.232;  Ly_f_nom = 0.323;  Ly_r_nom = 0.330;   % geometria nominal
    p.arms.Lx_r = Lx_nom   - dy;     % CG até os motores da direita (1 e 4)
    p.arms.Lx_l = Lx_nom   + dy;     % CG até os motores da esquerda (2 e 3)
    p.arms.Ly_f = Ly_f_nom - dx;     % CG até os motores da frente (1 e 3)
    p.arms.Ly_r = Ly_r_nom + dx;     % CG até os motores de trás (2 e 4)
    p.cg_offset = [dx; dy; 0];

    % [4] limite de saída do ArduPilot nesta campanha (MOT_SPIN_MAX = 0,95)
    p.motor.pwm_sat = 1950;

    % O tensor de inércia continua o do CAD: mesma estrutura, e não há medida
    % independente para esse dia. O deslocamento de CG de 45 mm muda o tensor
    % pelo teorema dos eixos paralelos em ~m·dx² = 1,993·0,0454² = 0,004 kg·m²
    % em Jy e Jz (5% e 3%), dentro da faixa que a identificação já explora.
end
