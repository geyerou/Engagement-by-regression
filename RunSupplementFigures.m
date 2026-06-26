function RunSupplementFigures()
code_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(code_dir,'functions'));
scripts={'FigureS01','FigureS02','FigureS03','FigureS04', ...
    'FigureS05','FigureS06','FigureS07','FigureS08'};
for i=1:numel(scripts)
    fprintf('\n=== Running %s (%d/%d) ===\n',scripts{i},i,numel(scripts));
    evalin('base',sprintf("run('%s')",fullfile(code_dir,[scripts{i} '.m'])));
end
end
