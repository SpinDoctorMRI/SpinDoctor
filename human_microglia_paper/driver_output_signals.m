clear type;
for i = 1:ncells
    file = names(i);
    sep_file = split(file,"/");
    type(i) = sep_file(3);
end

for i = 1:ncells
    output_signals(names(i),setup_file,sprintf('human_signals/%s',type(i)),sprintf('final_signals_matlab/%s',type(i)));
end