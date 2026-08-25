function lm = load_linear_model(paths)
%LOAD_LINEAR_MODEL  Carrega o modelo linear de hover salvo (outputs/linear_model.mat)
%   e AVISA se ele estiver desatualizado em relação ao P_identified.mat atual.
%
%   lm = load_linear_model()        % usa setup_paths() pra achar os caminhos
%   lm = load_linear_model(paths)   % reusa um struct paths já resolvido
%
%   Fonte ÚNICA de carregamento do modelo linear — usado por setup_quad_linear.m
%   e compare_v4_vs_validate.m. Garante que ambos usam EXATAMENTE o mesmo A, B,
%   x0, u0 (gerados uma vez por linearize.m), sem recalcular e sem divergir.
%
%   lm tem os campos salvos por linearize.m: A B C D x0 u0 sys eig_A P
%   pwm_trim M_alloc M_mix trim_residual rank_Co.
%   Entrada do modelo: v=[T,Mx,My,Mz] (forças); u0=[m·g;0;0;0]. A alocação
%   PWM⇄[T,M] vem em M_alloc (forward) e M_mix=inv(M_alloc) (controle→PWM).
%
%   Se o P usado na linearização != P_identified.mat atual → WARNING pedindo
%   pra rodar linearize.m de novo (evita usar um modelo linear stale por engano).

    if nargin < 1 || isempty(paths)
        paths = setup_paths();
    end

    lm_path = fullfile(paths.outputs, 'linear_model.mat');
    if ~exist(lm_path, 'file')
        error('load_linear_model:missing', ...
            ['linear_model.mat não existe em %s.\n' ...
             'Rode linearize.m primeiro pra gerar o modelo linear.'], lm_path);
    end
    lm = load(lm_path);

    % --- Checagem de stale: P salvo vs P_identified.mat atual ---
    pf = fullfile(paths.outputs, 'P_identified.mat');
    if isfield(lm, 'P') && exist(pf, 'file')
        Pcur = load(pf).P_final;
        if numel(Pcur) == numel(lm.P)
            dP = max(abs(Pcur(:) - lm.P(:)));
            if dP > 1e-9
                warning('load_linear_model:stale', ...
                    ['linear_model.mat foi gerado com um P DIFERENTE do ' ...
                     'P_identified.mat atual (max|ΔP|=%.3g).\n' ...
                     '         >>> Rode linearize.m pra atualizar o modelo linear. <<<'], dP);
            end
        else
            warning('load_linear_model:sizeMismatch', ...
                'P salvo (%d) e P_identified atual (%d) têm tamanhos diferentes.', ...
                numel(lm.P), numel(Pcur));
        end
    end
end
