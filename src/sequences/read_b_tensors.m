function b_tensors = read_b_tensors(file)
%READ_B_TENSORS Reads axis-symmetric b tensors from a file.
fid = fopen(file,'r');
fgetl(fid);
b_tensors = fscanf(fid,'%f',[5,Inf]);
fclose(fid);
end