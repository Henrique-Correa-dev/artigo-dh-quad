function [A, B, c] = linearize_at(x0, v0, P, constants)
%LINEARIZE_AT  Lineariza o CORPO RÍGIDO em torno de (x0, v0), entrada v=[T,Mx,My,Mz].
%
%   [A, B, c] = linearize_at(x0, v0, P, constants)
%
%   Aproximação de 1ª ordem:  ẋ ≈ f(x0,v0) + A·(x-x0) + B·(v-v0)
%       A = ∂f/∂x |_(x0,v0)   (9×9)   — Jacobiano numérico (diferenças centrais)
%       B = ∂f/∂v |_(x0,v0)   (9×4)   — CONSTANTE (corpo rígido é linear em [T,M])
%       c = f(x0,v0)          (9×1)   — termo afim (taxa de estado no ponto)
%
%   ENTRADA v = [T, Mx, My, Mz]  (empuxo total + 3 momentos), NÃO PWM.
%   A não-linearidade dos motores (curvas fT/fQ) fica FORA daqui — é a alocação.
%   Logo este B é geometria/inércia pura, INDEPENDENTE do ponto de operação:
%       ∂ṗ/∂Mx=G3, ∂ṗ/∂Mz=G4, ∂q̇/∂My=1/Jy, ∂ṙ/∂Mx=G4, ∂ṙ/∂Mz=G8, ∂ẇ/∂T=-1/m.
%   (Como f é exatamente linear em v, a diferença finita dá B exato p/ qualquer h.)
%
%   Usa o handle .rigid_body de vtol_dynamics (fonte única — não reimplanta G's).
%
%   Se (x0,v0) for equilíbrio (ex.: v0=[m·g;0;0;0] no hover), c≈0 e recai na
%   forma clássica δẋ = A·δx + B·δv.
%
%   Estados x = [p q r phi theta psi u v w] (9), entrada v = [T Mx My Mz] (4).

    x0 = x0(:);  v0 = v0(:);
    nx = numel(x0);  nv = numel(v0);

    % Perturbações de estado (escala física); v é exato em qualquer h (linear).
    h_x = [1e-4; 1e-4; 1e-4;   ... % p, q, r        (rad/s)
           1e-4; 1e-4; 1e-4;   ... % phi, theta, psi (rad)
           1e-3; 1e-3; 1e-3];      % u, v, w         (m/s)
    h_v = [1e-2; 1e-3; 1e-3; 1e-3];  % T (N), Mx, My, Mz (N·m)

    dyn = vtol_dynamics('get_handles');
    rb  = dyn.rigid_body;
    f = @(x, v) rb(x, v(1), v(2), v(3), v(4), P, constants);

    % Termo afim: taxa de estado no ponto de operação
    c = f(x0, v0);

    % A = ∂f/∂x (diferenças centrais)
    A = zeros(nx, nx);
    for j = 1:nx
        xp = x0; xp(j) = x0(j) + h_x(j);
        xm = x0; xm(j) = x0(j) - h_x(j);
        A(:, j) = (f(xp, v0) - f(xm, v0)) / (2 * h_x(j));
    end

    % B = ∂f/∂v (diferenças centrais — exato, pois f é linear em v)
    B = zeros(nx, nv);
    for j = 1:nv
        vp = v0; vp(j) = v0(j) + h_v(j);
        vm = v0; vm(j) = v0(j) - h_v(j);
        B(:, j) = (f(x0, vp) - f(x0, vm)) / (2 * h_v(j));
    end
end
