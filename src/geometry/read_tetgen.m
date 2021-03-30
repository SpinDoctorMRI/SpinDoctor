function femesh_all = read_tetgen(filename)
%READ_TETGEN Read mesh from Tetgen.
%
%   filename: string
%
%   femesh_all: struct with fields
%       points_all: [3 x npoint]
%       facets_all: [3 x nfacet]
%       elements_all: [4 x nelement]
%       facetmarkers: [1 x nfacet]
%       elementmarkers: [1 x nelement]


disp("Reading from Tetgen FE mesh from " + filename);

% Read nodes
fid = fopen(filename + ".node", "r");
if fid ~= -1
    fscanf(fid, "%d", [4, 1]);
    vec = fscanf(fid, "%f", [4, Inf]);
end
fclose(fid);
points = vec(2:end, :);

% Read facets and their associated boundaries
fid = fopen(filename + ".face", "r");
if fid ~= -1
    fscanf(fid, "%d", [2, 1]);
    vec = fscanf(fid, "%f", [5, inf]);
end
fclose(fid);
facets = vec(2:4, :);
facetmarkers = vec(5, :);

% Read elements and their associated compartments
fid = fopen(filename + ".ele", "r");
if fid ~= -1
    fscanf(fid, "%d", [3, 1]);
    vec = fscanf(fid, "%f", [6, inf]);
end
fclose(fid);
elements = vec(2:5, :);
elementmarkers = vec(6, :);

% Create output structure
femesh_all.points = points;
femesh_all.facets = facets;
femesh_all.elements = elements;
femesh_all.facetmarkers = facetmarkers;
femesh_all.elementmarkers = elementmarkers;
