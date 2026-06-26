function RunAllFigures()
% Generate publication panels for all main and supplementary figures.
code_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(code_dir,'functions'));
scripts={'Figure01','Figure02','Figure03','Figure04','Figure05', ...
    'FigureS01','FigureS02','FigureS03','FigureS04','FigureS05', ...
    'FigureS06','FigureS07','FigureS08'};
for i=1:numel(scripts)
    fprintf('\n=== Running %s (%d/%d) ===\n',scripts{i},i,numel(scripts));
    evalin('base',sprintf("run('%s')",fullfile(code_dir,[scripts{i} '.m'])));
end
fprintf('\nAll figure panels complete. Output root: %s\n', ...
    fullfile(fig_project_root(),'figures'));
end
