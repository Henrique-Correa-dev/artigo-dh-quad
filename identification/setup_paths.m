function paths = setup_paths()
%SETUP_PATHS  Resolve caminhos absolutos do projeto e adiciona subpastas ao path.
%
% USO:
%   paths = setup_paths();
%   load(fullfile(paths.data, 'log_data.mat'));
%   save(fullfile(paths.outputs, 'P_identified.mat'), 'P_final');
%
% Chame no início de qualquer script. Idempotente (pode chamar várias vezes).

    here = fileparts(mfilename('fullpath'));   % raiz do projeto identification/

    paths.root           = here;
    paths.data           = fullfile(here, '1_data');
    paths.model          = fullfile(here, '2_model');
    paths.identification = fullfile(here, '3_identification');
    paths.simulink       = fullfile(here, '4_simulink');
    paths.validation     = fullfile(here, '5_validation');
    paths.linear         = fullfile(here, '6_linear');
    paths.outputs        = fullfile(here, 'outputs');
    paths.images         = fullfile(paths.outputs, 'images');
    paths.reference      = fullfile(here, 'reference');
    paths.docs           = fullfile(here, 'docs');
    paths.legacy         = fullfile(here, 'legacy');

    % Adicionar todas as subpastas ATIVAS ao MATLAB path
    % (legacy/, docs/, reference/ NÃO entram — código antigo / docs / refs externas)
    active_folders = {paths.root, paths.data, paths.model, paths.identification, ...
                      paths.simulink, paths.validation, paths.linear, ...
                      paths.outputs};

    for k = 1:numel(active_folders)
        if exist(active_folders{k}, 'dir')
            addpath(active_folders{k});
        end
    end
end
