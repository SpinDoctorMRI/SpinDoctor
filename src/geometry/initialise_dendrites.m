function [femesh_dendrites,element_markers]= initialise_dendrites(p,e,bins,dendrite_elements)
%%INITIALISE_DENDRITES creates femesh structures for the dendrites, after they have been sorted into branches
% 
%   p: (3,N) array of points from full cell finite element mesh
%   e: (3,N_) array of elements from full cell finite element mesh
%   bins: (N__,1) array with each entry corresponding to the branch id of
%       the corresponding element in dendrite_elements
%   dendrite_elements: (N__,1) array of the elements from e which are to be
%       sorted into dendrite branches
    dendrite_ids = unique(bins);
    femesh_dendrites = cell(length(dendrite_ids),1);
    element_markers = zeros(1,size(e,2));
    for i = 1:length(dendrite_ids)
        femesh_dendrite = struct;
        femesh_dendrite.ncompartment =1;
        femesh_dendrite.nboundary =1;
        element_map = dendrite_elements(bins == dendrite_ids(i));   
        element_markers(element_map) = i + 1;
        elements = e(:,element_map);
        point_map = unique(elements);
        points = p(:,point_map);
        femesh_dendrite.points = {points};
        mask = false([size(points,2),1]);
        mask(point_map) = true;
        indices_map = cumsum(mask);
        elements = indices_map(elements);
    
        femesh_dendrite.elements = {elements};
        femesh_dendrite.point_map = {point_map};
        femesh_dendrite.element_map = element_map;
        mesh = triangulation(elements',points');
        femesh_dendrite.facets = {freeBoundary(mesh)'};
        [volumes, areas] = get_vol_sa(femesh_dendrite);
        femesh_dendrite.volumes = volumes;
        femesh_dendrite.total_volume = sum(volumes, 'all');
        femesh_dendrite.areas = areas;
        femesh_dendrite.total_area = sum(areas, 'all');
        femesh_dendrites{dendrite_ids(i)} = femesh_dendrite;
    end