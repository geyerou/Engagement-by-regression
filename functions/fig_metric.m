function value=fig_metric(csv_file,name)
T=readtable(csv_file,'Delimiter',',','TextType','string', ...
    'VariableNamingRule','preserve');
m=string(T{:,1}); v=T{:,2};
if ~isnumeric(v), v=str2double(string(v)); end
idx=find(m==string(name),1);
if isempty(idx), value=NaN; else, value=double(v(idx)); end
end
