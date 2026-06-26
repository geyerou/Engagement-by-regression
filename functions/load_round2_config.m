function cfg = load_round2_config()
project_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
f=fullfile(project_dir,'current_config_round2.mat');
if ~exist(f,'file'), error('Run s100_config_round2 first.'); end
S=load(f,'cfg'); cfg=S.cfg; addpath(cfg.functions_dir);
end
