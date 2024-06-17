function femesh_soma = initialise_soma(points,e,soma_element_map)
%%INITIALISE_SOMA given the ids of elements of the global mesh which forms
%%the soma, the finite element mesh for the soma is found.
% 
%   points: (3,N) array of x,y,z positions of points in the full cell
%       finite element mesh
%   e: (4,N_) array of tetrahedra in full cell finite elemen mesh
%   soma_element_map: (N__,1) array of ids of the elements in e which should 
%       be taken to form the soma finite element mesh
% 
%   femesh_soma: finite element mesh of soma

    % Extract soma femesh from a triangulation with a mask
    mask = false([size(points,2),1]);
    mask(unique(e(:,soma_element_map))) = true;
    %Storing information within the femesh
    soma_points = find(  mask );
    dendrite_points = find(~mask);
    femesh.soma_points = {soma_points};
    femesh.dendrite_points = {dendrite_points};
    indices_map = cumsum(mask);
    % Create the suface mesh and geometrical information of the soma.
    e =e(:,soma_element_map);
    e = indices_map(e);
    p = points(:,soma_points);
    warning('off','all');
    soma = triangulation(e',p');
    warning('on','all');
    femesh_soma.points ={soma.Points'};
    femesh_soma.elements = {soma.ConnectivityList'}; 
    femesh_soma.facets = {freeBoundary(soma)'};
    femesh_soma.ncompartment = 1;
    femesh_soma.nboundary = 1;
    femesh_soma.element_map = soma_element_map;

    [volumes, areas] = get_vol_sa(femesh_soma);
    femesh_soma.volumes = volumes;
    femesh_soma.total_volume = sum(volumes, 'all');
    femesh_soma.areas = areas;
    femesh_soma.total_area = sum(areas, 'all');

    femesh_soma.point_map = {[1:size(femesh_soma.points{1},2)]};
