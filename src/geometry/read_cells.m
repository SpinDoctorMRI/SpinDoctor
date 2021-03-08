function cells = read_cells(cellfilename)
%READ_CELLS_INPUTFILE Read cell configuration (centers and radii).
%
%   cellfilename: string
%
%   cells: struct


% Read geometry from file
disp("Loading cell configuration from " + cellfilename);
fid = fopen(cellfilename, "r");
fgetl(fid);
tline = fgetl(fid);
ncell = sscanf(tline, "%d", 1);
fgetl(fid);
tline = fgetl(fid);
d = sscanf(tline, "%d", 1);
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
