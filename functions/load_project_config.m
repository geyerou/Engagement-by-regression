function cfg = load_project_config()
code_dir = fileparts(fileparts(mfilename('fullpath')));
project_dir = fileparts(code_dir);
% s000 writes the active analysis configuration here. This allows separate
% derivative trees for different smoothing/preprocessing variants.
config_file = fullfile(project_dir,'current_config.mat');
if ~exist(config_file, 'file')
    % Backward-compatible fallback for projects created before analysis
    % branches were introduced.
    config_file = fullfile(project_dir, 'derivatives', ...
        'result_s000_config_project', 'config.mat');
end
if ~exist(config_file, 'file')
    error('Configuration not found. Run s000_config_project.m first.');
end
tmp = load(config_file, 'cfg');
cfg = tmp.cfg;
addpath(cfg.functions_dir);
end
