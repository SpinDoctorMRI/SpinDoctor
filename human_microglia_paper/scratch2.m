addpath(genpath('setups'));
addpath(genpath('src'));
fid = fopen("cells_human_separated.txt",'r');
tline = fgetl(fid);
i=1;
clear names;
while ischar(tline)
    names(i) = string(tline);
    tline = fgetl(fid);
    i = i + 1;
end
ncells = length(names);
clear type;
for i = 1:ncells
    file = names(i);
    sep_file = split(file,"/");
    type(i) = sep_file(3);
end
for i = 1:ncells
    output_signals(names(i),sprintf('human_signals/%s',type(i)),sprintf('final_signals_matlab/%s',type(i)));
end