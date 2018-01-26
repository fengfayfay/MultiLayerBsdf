function gm2pbrtinput(dir,obj)
filename = [dir,'weights.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.ComponentProportion);
fclose(fileID);

filename = [dir,'means.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.mu);
fclose(fileID);

filename = [dir,'covars.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.Sigma);
fclose(fileID);
end