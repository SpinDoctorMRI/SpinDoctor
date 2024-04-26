function sequences = read_scheme(filename)

[~,name,~] = fileparts(filename);
fid = fopen(filename,'r');
tline = fgetl(fid);
s = string(tline);
if s ~= "VERSION: GRADIENT_WAVEFORM"
    error(sprintf("%s not formatted correctly. Must have %s as first line",filename,"VERSION: GRADIENT_WAVEFORM"));
end
tline = fgetl(fid);
sequences = {};
nseqeunce = 0;
while ischar(tline)
    vec = str2num(tline);
    nseqeunce = nseqeunce + 1;
    K = vec(1); dt = 1e6*vec(2);
    g_x = vec(3*(1:K)); g_y = vec(1+3*(1:K)); g_z = vec(2+3*(1:K)); 
    g = 1e3*[g_x;g_y;g_z];


    sequence = SequenceCamino(K,dt,g,sprintf("%s_%d",name,nseqeunce));
    sequences{nseqeunce} = sequence;
    tline = fgetl(fid);

end

fclose(fid);

