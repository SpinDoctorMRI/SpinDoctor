function save_cells(cells, cellfilename)
%SAVE_CELLS Write geometrical cell configuration to file.
%   The file contains the spatial dimension: 2 for cylinders, 3 for spheres
%
%   cells: struct
%   cellfilename: string


% Extract cell parameters
centers = cells.centers;
radii = cells.radii;

% Sizes
ncell = length(radii);
dim = size(centers, 1);

% Save geometry to file
disp("Writing geometry to " + cellfilename);
fid = fopen(cellfilename, "w");
fprintf(fid, "Number of cells:\n");
fprintf(fid, "%d\n", ncell);
fprintf(fid, "Dimension:\n");
fprintf(fid, dim + "\n");
c = "x,y";
if dim == 3
    c = c + ",z";
end
fprintf(fid, "Number n, Center (%s), Radius r\n", c);
for icell = 1:ncell
    if fid ~= -1
        fprintf(fid, "%d " + join(repmat("%g", 1, dim)) + " %g\n", icell, centers(:, icell), radii(icell));
    end
end
fclose(fid);
