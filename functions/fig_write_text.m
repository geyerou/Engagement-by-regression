function fig_write_text(file_path,txt)
fid=fopen(file_path,'w');
if fid<0, error('Cannot write %s',file_path); end
cleanup=onCleanup(@()fclose(fid));
txt=char(txt);
txt=strrep(txt,'\n',newline);
fprintf(fid,'%s',txt);
end
