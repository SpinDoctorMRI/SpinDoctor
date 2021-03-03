function surfaces = read_surfaces(filename)
%READ_SURFACES Load surface triangulation.
%
%   filename: string
%   surfaces: struct


% Get list of nodes from .node file
fid = fopen(filename + ".node", "r");
fgetl(fid);
fgetl(fid);
tline = fgetl(fid);
npoint = sscanf(tline, "%d", 1);
fgetl(fid);
points = zeros(3, npoint);
for ipoint = 1:npoint
    tline = fgetl(fid);
    vec = sscanf(tline, "%f");
    points(:, ipoint) = vec(2:4);
end
fclose(fid);

% Read .poly file
fid = fopen(filename + ".poly", "r");

% Read list of holes (refer to separate file)
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Read list of facets
fgetl(fid);
fgetl(fid);
tline = fgetl(fid);
nfacet = sscanf(tline, "%d", 1);
fgetl(fid);
facets = zeros(3, nfacet);
facetmarkers = zeros(1, nfacet);
for ifacet = 1:nfacet
    tline = fgetl(fid);
    vec = sscanf(tline, "%d");
    facetmarkers(ifacet) = vec(3);
    tline = fgetl(fid);
    vec = sscanf(tline, "%d");
    facets(:, ifacet) = vec(2:4);
end

% Read list of holes (empty)
fgetl(fid);
fgetl(fid);

% Read list of interior points with their corresponding compartment
fgetl(fid);
tline = fgetl(fid);
nregion = sscanf(tline, "%d", 1);
regions = zeros(3, nregion);
for iregion = 1:nregion
    tline = fgetl(fid);
    vec = sscanf(tline, "%f");
    regions(:, iregion) = vec(2:4);
end

fclose(fid);

% Extract surface triangulation
surfaces.points = points;
surfaces.facets = facets;
surfaces.facetmarkers = facetmarkers;
surfaces.regions = regions;
