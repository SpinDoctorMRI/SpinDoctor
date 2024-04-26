function b_tensors = read_b_tensors(file)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(file,'r');
fgetl(fid);
b_tensors = fscanf(fid,'%f',[5,Inf]);
fclose(fid);
end