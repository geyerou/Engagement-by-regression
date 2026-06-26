function T = read_jhu_label_names(xml_file)
xml_text=fileread(xml_file);
tokens=regexp(xml_text,'<label index="(\d+)"[^>]*>(.*?)</label>','tokens');
id=[]; name=strings(0,1);
for i=1:numel(tokens)
    k=str2double(tokens{i}{1});
    if k==0, continue; end
    id(end+1,1)=k; %#ok<AGROW>
    name(end+1,1)=strtrim(string(tokens{i}{2})); %#ok<AGROW>
end
T=table(id,name,'VariableNames',{'label_id','label_name'});
end
