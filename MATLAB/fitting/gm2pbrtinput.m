function gm2pbrtinput(dir,obj)
filename = [dir,'weights.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.ComponentProportion);
fclose(fileID);


filename = [dir,'means.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.mu);
fclose(fileID);
disp(size(obj.mu));
disp(obj.mu);

filename = [dir,'covars.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%4.6f\n',obj.Sigma);
fclose(fileID);
disp(size(obj.Sigma));
disp(obj.Sigma);
end
