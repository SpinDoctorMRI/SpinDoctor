function [volumes, surface_areas] = get_vol_sa(femesh)
%GET_VOL_SA Get volume and surface area quantities.
%
%   femesh: struct
%   direction: double(3, 1)
%
%   volumes: double(1, ncompartment)
%  	surface_areas: double(1, ncompartment)


% Number of compartments
ncompartment = femesh.ncompartment;

% Initialize output arguments
volumes = zeros(1, ncompartment);
surface_areas = zeros(1, ncompartment);

for icmpt = 1:ncompartment
    % FE Mesh
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    facets = [femesh.facets{icmpt, :}];

    % Volumes
    volumes(icmpt) = get_volume_mesh(points, elements);

    % Surface area
    surface_areas(icmpt) = get_surface_mesh(points, facets);
end
