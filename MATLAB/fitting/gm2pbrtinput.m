function gm2pbrtinput(dir,obj,reflect)
if reflect
    add = '_reflect';
else
    add ='';
end
filename = [dir,'weights',add,'.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.ComponentProportion);
fclose(fileID);

filename = [dir,'means',add,'.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.mu);
fclose(fileID);

filename = [dir,'covars',add,'.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.Sigma);
fclose(fileID);
end
