function cells = read_cells_inputfile(cellfilename)
%READ_CELLS_INPUTFILE Create geometrical configuration and write to file.
%
%   cellfilename: string


% Read geometry from file
disp("Loading cell configuration from " + cellfilename);
fid = fopen(cellfilename, "r");
fgetl(fid);
tline = fgetl(fid);
ncell = sscanf(tline, "%d", 1);
fgetl(fid);
tline = fgetl(fid);
shape = sscanf(tline, "%s", 1);
if shape == "sphere"
    d = 3;
else
    d = 2;
end
centers = zeros(d, ncell);
radii = zeros(1, ncell);
fgetl(fid);
for icell = 1:ncell
    tline = fgetl(fid);
    vec = sscanf(tline, "%f");
    centers(:, icell) = vec(2:1+d);
    radii(icell) = vec(2+d);
end
fclose(fid);

% Create output structure
cells.centers = centers;
cells.radii = radii;
