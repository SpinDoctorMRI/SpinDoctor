classdef Soma 
    properties
        points % point cloud on the surface
        normals % normals on surface
        mask % mask of points to keep
        centre 
        radius 
        volume % volume of segment
        bbox % bounding box of segment
        index % index of node
        parent_index % index of parent node
        type % type of node
        area % surface_area of the soma
      end
    methods
        function soma = Soma(centre,radius,connectivity,type)
            soma.centre = centre;
            soma.radius = radius;
            soma.volume = (4*pi*radius^3)/3;
            soma.area = 4*pi*radius^2;
            soma.bbox = soma.make_bbox();
            soma.index =connectivity(1);
            soma.parent_index = connectivity(2);
            soma.type = type;
        end
        function bbox = make_bbox(soma)
            R = soma.radius;
            bbox = R*[-1,1,-1,-1,-1,1,1,1;-1,1,1,-1,1,-1,-1,1;-1,1,1,1,-1,-1,1,-1]+ soma.centre';
        end
        function soma =  make_point_cloud(soma,density)
            % makes solid 3d point cloud
            if nargin < 2
                density = 1.0;
            end
            n = floor(soma.area*density);
            [x,y,z] = unitsphere(n);
            p = soma.radius*[x;y;z]'+ soma.centre;
            soma.points = p;
            soma.normals = [x;y;z]';
            soma.mask = true([n,1]);
        end
        function [in, on, outer, out_near] =intersect(soma,seg)
            p = seg.points;
            dist = vecnorm(p - soma.centre, 2,2) - soma.radius;
            
            m = true([length(p),1]);
            [in, on, outer, out_near] = create_mask(m ,dist,soma.radius);
        end


        function new_mask = intersect_old(soma,seg)
            % detects all points in seg which intersect obj and updates
            % seg mask
            outside = vecnorm( seg.points - soma.centre,2,2) > soma.radius;
            new_mask = outside & seg.mask;
        end
        function [p,n] = output_filtered_points(soma)
            % updates filtered points according to the mask
            p = soma.points(soma.mask,:);
            n = soma.normals(soma.mask,:);
        end
       
        function plot(soma,color)
            p = output_filtered_points(soma);
            if nargin <2
                color = soma.type/8;
            end
            x = p(:,1);
            y = p(:,2);
            z = p(:,3);
            scatter3(x,y,z,color);
        end

    end
end

