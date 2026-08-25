%TEST_TM_REFACTOR  Testa o refactor PWM → [T,Mx,My,Mz] (alocação + corpo rígido).
%
% Roda checagens PASS/FAIL em duas frentes:
%   (a) REGRESSÃO  — o split não mudou a dinâmica (forces→rigid_body == vtol_dynamics);
%                    sim_window roda nos 3 modos.
%   (b) MODELO LINEAR [T,M] — B analítico/constante, A com polos de amortecimento,
%                    c=0 no hover, mixer inversível, .slx ≈ .m.
%
% Pré-requisito: rode linearize.m antes (gera outputs/linear_model.mat atualizado).
%
% Uso:  >> test_tm_refactor

clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = setup_paths();

proj = parameters();
constants.m = proj.m;  constants.g = proj.g;
P_J = load(fullfile(paths.outputs,'P_identified.mat')).P_final;
[fT, fQ] = motor_models();
dyn = vtol_dynamics('get_handles');

np = 0; nf = 0;
check = @(name, cond) fprintf('  [%s] %s\n', tern(cond,'PASS','FAIL'), name);
% (contadores via subfunção no fim)

fprintf('\n=== (a) REGRESSÃO: split alocação/corpo-rígido ===\n');

% T1: composição forces→rigid_body reproduz vtol_dynamics (3 e 9 estados)
tt = (0:0.1:2)'; pwm = [1700+30*sin(tt), 1690-20*tt, 1660+10*cos(tt), 1680+5*tt];
y9 = [0.2;-0.1;0.05; 0.1;-0.05;0.3; 0.4;-0.2;0.15];
d_full = vtol_dynamics(1.3, y9, P_J, tt, pwm, fT, fQ, constants);
pwm_t  = interp1(tt, pwm, 1.3, 'linear', 'extrap');
[Tt,Mx,My,Mz] = dyn.forces(pwm_t, P_J, fT, fQ);
d_split = dyn.rigid_body(y9, Tt, Mx, My, Mz, P_J, constants);
e1 = max(abs(d_full - d_split));
check(sprintf('forces->rigid_body == vtol_dynamics (max|Δ|=%.1e)', e1), e1 < 1e-12);
[np,nf] = tally(e1<1e-12, np, nf);

% T2: forces vetorizado (Nx4) bate com loop escalar
[Tv,Mxv,Myv,Mzv] = dyn.forces(pwm, P_J, fT, fQ);
ok2 = isequal(size(Tv),[size(pwm,1) 1]) && abs(Tv(14)-sum(P_J(5:8).*fT(pwm(14,:))'))<1e-9;
check('forces vetorizado (Nx4 -> Nx1)', ok2);  [np,nf]=tally(ok2,np,nf);

fprintf('\n=== (b) MODELO LINEAR [T,M] ===\n');
lm = load_linear_model(paths);

% T3: B analítico (G3, G4, 1/Jy, G8, -1/m)
Jx=P_J(1);Jy=P_J(2);Jz=P_J(3);Jxz=P_J(4); g0=Jx*Jz-Jxz^2;
Bana=zeros(9,4);
Bana(1,2)=Jz/g0;  Bana(1,4)=Jxz/g0;      % p: Mx,Mz
Bana(2,3)=1/Jy;                           % q: My
Bana(3,2)=Jxz/g0; Bana(3,4)=Jx/g0;        % r: Mx,Mz
Bana(9,1)=-1/constants.m;                 % w: T
e3 = max(abs(lm.B(:)-Bana(:)));
check(sprintf('B == analítico (max|Δ|=%.1e)', e3), e3 < 1e-9);  [np,nf]=tally(e3<1e-9,np,nf);

% T4: -diag(A(1:3)) == [Dp Dq Dr]
e4 = max(abs(-diag(lm.A(1:3,1:3)) - P_J(13:15)));
check(sprintf('A: -diag(1:3) == [Dp Dq Dr] (max|Δ|=%.1e)', e4), e4 < 1e-6);  [np,nf]=tally(e4<1e-6,np,nf);

% T5: c=0 no hover (trim de força é equilíbrio)
e5 = norm(lm.trim_residual);
check(sprintf('c = f(x0,v0) ~ 0  (|c|=%.1e)', e5), e5 < 1e-6);  [np,nf]=tally(e5<1e-6,np,nf);

% T6: u0 = [mg,0,0,0]
e6 = norm(lm.u0 - [constants.m*constants.g;0;0;0]);
check(sprintf('u0 == [m·g;0;0;0] (Δ=%.1e)', e6), e6 < 1e-6);  [np,nf]=tally(e6<1e-6,np,nf);

% T7: mixer inversível  M_alloc * M_mix == I
e7 = max(abs(reshape(lm.M_alloc*lm.M_mix - eye(4),[],1)));
check(sprintf('M_alloc*M_mix == I (max|Δ|=%.1e, cond=%.1f)', e7, cond(lm.M_alloc)), e7 < 1e-6);
[np,nf]=tally(e7<1e-6,np,nf);

% T8: controlável
ok8 = (rank(ctrb(lm.A,lm.B)) == 9);
check('sistema controlável (rank ctrb = 9)', ok8);  [np,nf]=tally(ok8,np,nf);

fprintf('\n=== RESULTADO: %d PASS, %d FAIL ===\n', np, nf);
if nf==0
    fprintf('TUDO VERDE. Pra o teste de integração (sim_window 3 modos + .slx vs .m),\n');
    fprintf('rode: validate_params  |  compare_v4_vs_validate  |  setup_quad_linear+sim\n');
else
    fprintf('>>> %d falha(s) — investigar antes de prosseguir. <<<\n', nf);
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
function [np,nf] = tally(ok,np,nf), if ok, np=np+1; else, nf=nf+1; end, end
