function [femesh,femesh_soma,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path,soma_mesh_path)
%%SEGMENT_FEMESH creates finite element meshes corresponding to the soma
%%and dendrite branches.
%   
%   femesh: finite element mesh of the full cell
%   swc_file: swc file of the full cell.
%   tetgen_path: location and name of tetgen files corresponding to femesh
%   
%   femesh_soma: finite element mesh of soma
%   femesh_dendrites: cell array of finite element meshes for dendrites

    
    % Check for pre-existing decomposition
    if isfile(sprintf('%s_soma_elements',tetgen_path))
        disp('Reading from previous segmentation')
        p = [femesh.points{:}];
        e = [femesh.elements{:}];
        npoint = size(p,2); nelement = size(e,2);
        fid = fopen(sprintf('%s_soma_elements',tetgen_path),'r');
        soma_element_map=fscanf(fid,'%d',[1,Inf]);
        fclose(fid);
        
        femesh_soma= initialise_soma(p,e,soma_element_map);

        i = 1;
        bins = zeros([nelement,1]);
        while isfile(sprintf('%s_dendrite_%d_elements',tetgen_path,i))
            fid = fopen(sprintf('%s_dendrite_%d_elements',tetgen_path,i),'r');
            element_map = fscanf(fid,'%d',[1,Inf]);
            fclose(fid);
            bins(element_map) = i;
            i = i+1;
        end
        mask = true([nelement,1]);
        mask(femesh_soma.element_map) = false;
        dendrite_elements = find(mask);
        bins = bins(bins>0);
        [femesh_dendrites,elementmarkers] = initialise_dendrites(p,e,bins,dendrite_elements);
        
        % Previous method, using comparts in one finite element mesh to
        % handle segmentation.
        % elementmarkers(soma_element_map) = 1;
        % femesh.elementmarkers = elementmarkers;
        % femesh.points = [femesh.points{:}]; femesh.facets = [femesh.facets{:}]; femesh.elements = [femesh.elements{:}];
        % femesh.facetmarkers = ones(1,size(femesh.facets,2));
        % femesh = split_mesh(femesh);
    else
        disp('Segmenting finite element mesh')
        % Extract swc information and tetrahedra adjacency
        swc = Swc(swc_file);neighbours =read_tetgen_neigh(string(tetgen_path));
        % Create soma mesh
        if nargin ==4
             [femesh_soma,soma_element_map] = separate_soma(femesh,swc,neighbours,soma_mesh_path);
        else
             [femesh_soma,soma_element_map] = separate_soma(femesh,swc,neighbours);
        end
    
        % Create dendrites
        [femesh_dendrites,elementmarkers] = find_dendrites_femesh(femesh,femesh_soma,neighbours);
        ndendrites = length(femesh_dendrites);
        
        % Previous method, using comparts in one finite element mesh to
        % handle segmentation.
        % elementmarkers(soma_element_map) = 1;
        % femesh.elementmarkers = elementmarkers;
        % femesh.facetmarkers = ones(1,size(femesh.facets,2));
        % femesh = split_mesh(femesh);
        % Save soma mesh
        fid = fopen(sprintf('%s_soma_elements',tetgen_path),'w');
        fprintf(fid,'%d\n',femesh_soma.element_map);
        fclose(fid);
    
        % Save dendrites
        for i =1:ndendrites
            fid = fopen(sprintf('%s_dendrite_%d_elements',tetgen_path,i),'w');
            fprintf(fid,'%d\n',femesh_dendrites{i}.element_map);
            fclose(fid);
        end


    end