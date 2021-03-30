function save_surfaces(filename, surfaces)
%SAVE_SURFACES Save surface triangulation.
%   The files may be passed to Tetgen with `filename`.node.
%
%   filename: string
%   surfaces: struct


default_refinement = 0.1;

% Extract surface triangulation
points = surfaces.points;
facets = surfaces.facets;
facetmarkers = surfaces.facetmarkers;
regions = surfaces.regions;

npoint = size(points, 2);
nfacet = size(facets, 2);
nregion = size(regions, 2);

% Write list of nodes in .node file
fid = fopen(filename + ".node", "w");
fprintf(fid, "# Part 1 - node list\n");
fprintf(fid, "# node count, 3D, no attribute, no boundary marker\n");
fprintf(fid, "%d %d %d %d\n", npoint, 3, 0, 0);
fprintf(fid, "# Node index, node coordinates\n");
if fid ~= -1
    fprintf(fid, "%d %26.16f %26.16f %26.16f\n", [1:npoint; points]);
end
fclose(fid);

% Write .poly file
fid = fopen(filename + ".poly", "w");

% Write list of holes (refer to separate file)
fprintf(fid, "# Part 1 - node list\n");
fprintf(fid, "#  0 indicates the node list is stored in file .node\n");
fprintf(fid, "0\n");

% Write list of facets
fprintf(fid, "# Part 2 - facet list\n");
fprintf(fid, "# facet count, yes boundary marker\n");
fprintf(fid, "%d %d\n", nfacet, 1);
fprintf(fid, "# Node index, node coordinates\n");
if fid ~= -1
    for ifacet = 1:nfacet
        fprintf(fid, "%d %d %d\n", 1, 0, facetmarkers(ifacet));
        fprintf(fid, "%d %d %d %d \n", 3, facets(:, ifacet));
    end
end

% Write list of holes (empty)
fprintf(fid, "# Part 3 - hole list\n");
fprintf(fid, "0\n");

% Write list of interior points with their corresponding compartment
fprintf(fid, "# Part 4 - region list\n");
fprintf(fid, "%d\n", nregion);
for iregion = 1:nregion
    fprintf(fid, "%d %f %f %f %d %f\n", iregion, regions(:, iregion), iregion, default_refinement);
end

fclose(fid);
