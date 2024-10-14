classdef Segment 
    properties
        points % point cloud on the surface
        normals % normals to point cloud
        mask % mask of points to keep
        centres % 2,3 array of centres
        radii % 2 array of radii
        length % length of segment
        volume % volume of segment
        area % surface area of segment
        axis % direction of segment
        bbox % bounding box of segment
        index % index of node
        parent_index % index of parent node
        type % type of node
        rotation_matrix % matrix to align z axis from centre 1 to centre 2
        include_caps
    end

    methods
        function obj = Segment(centres,radii,connectivity,type,include_caps)
            obj.centres = centres;
            obj.radii = radii;
            obj.length = norm(centres(1,:) - centres(2,:));
            obj.axis = (centres(2,:) - centres(1,:))/obj.length;
            obj.rotation_matrix = make_rotation_matrix(obj.axis);
            [volume,area] = obj.get_geometric_measures();
            obj.volume = volume;
            obj.area = area;
            obj.bbox = obj.make_bbox();
            obj.index =connectivity(1);
            obj.parent_index = connectivity(2);
            obj.type = type;
            if nargin == 4
                obj.include_caps = true;
            end
            obj.include_caps = include_caps;
        end
        function bbox = make_bbox(obj)
            R = obj.length/2 + max(obj.radii);
            bbox = R*[-1,1,-1,-1,-1,1,1,1;-1,1,1,-1,1,-1,-1,1;-1,1,1,1,-1,-1,1,-1]+ mean(obj.centres,1)';
        end

        function [volume,area]  = get_geometric_measures(obj)
            R = max(obj.radii);
            r = min(obj.radii);
            if R== r
                volume = obj.length*pi*R^2;
                area = 2*pi*R*obj.length + 4*pi*R^2;
            else
                h = R*obj.length/(R - r);
                volume = (h*pi*R^2 - (h - obj.length)*pi*r^2)/3;
                area  = pi*R*sqrt(h^2+R^2) - pi*r*sqrt((h-obj.length)^2+r^2) + 2*pi*(R^2 + r^2);
            end
        end
        
        function obj = make_point_cloud(obj,density)
            if nargin < 2
                density = 1;
            end
            n = floor(obj.area*density);
            if obj.include_caps
                n1 =  floor(2*pi*obj.radii(1)^2*density);
                n2 =  floor(2*pi*obj.radii(2)^2*density);
                [points,theta,normals] = local_frustrum(n - n1 -n2,obj.radii(1),obj.radii(2),obj.length);
                [x,y,z] = unitsphere(2*n1);
                
                points1 = [x;y;z];
                I1 = [0 0 1]*points1 <= 0;
                bottom_cap = obj.radii(1)*( obj.rotation_matrix*points1(:,I1) ) + obj.centres(1,:)' ;
                [x,y,z] = unitsphere(2*n2);
                points2 = [x;y;z];
                I2 = [0 0 1]*points2 >= 0;
                top_cap = obj.radii(2)*( obj.rotation_matrix*points2(:,I2) ) + obj.centres(2,:)' ;
                surface =  obj.rotation_matrix*points  + obj.centres(1,:)';
                p = [surface,bottom_cap,top_cap];
                
                obj.points = p';
                obj.normals = [obj.rotation_matrix*normals, obj.rotation_matrix*points1(:,I1),obj.rotation_matrix*points2(:,I2)]';
                obj.mask = true([length(p),1]);
            else
                [points,theta,normals] = local_frustrum(n,obj.radii(1),obj.radii(2),obj.length);
                 % obj.points = (obj.rotation_matrix*points)'  + obj.centres(1,:);
                 % obj.normals = (obj.rotation_matrix*normals)';
                 % obj.mask = true([length(points),1]);
                    % Changed on 21/02
                n1 =  floor(2*pi*obj.radii(1)^2*density);
                [x,y,z] = unitsphere(2*n1);
                points1 = [x;y;z];
                I = [0 0 1]*points1 <= 0;
                bottom_cap = obj.radii(1)*( obj.rotation_matrix*points1(:,I) ) + obj.centres(1,:)' ;
                frustrum = (obj.rotation_matrix*points)  + obj.centres(1,:)';
                p = [frustrum,bottom_cap];
                
                obj.points = p';
                obj.normals = [obj.rotation_matrix*normals, obj.rotation_matrix*points1(:,I)]';
                obj.mask = true([length(p),1]);                
            end
        end

        function [p,n] = output_filtered_points(obj)
            % updates filtered points according to the mask
            p = obj.points(obj.mask,:);
            n = obj.normals(obj.mask,:);
        end
       
        function plot(obj,color)
            % p = output_filtered_points(obj);
            p = obj.points;
            if nargin <2
                color = obj.type/8;
            end
            x = p(:,1);
            y = p(:,2);
            z = p(:,3);
            scatter3(x,y,z,color);
        end
        
        function [in, on, outer, out_near] =intersect(obj,p,eps)
            [N,M] = size(p);
            r_min = min([obj.radii,obj.length]);
            if M ~=3
                p = p';
                N = size(p,1);
            end

         
            z = (p -obj.centres(1,:))*obj.axis';            
            lateral_mask = (z>=0) & (z<= obj.length);
            q = p - obj.centres(1,:);
            q = q -z*obj.axis;
            l = z/obj.length;
            r = (1-l)*obj.radii(1) + l*obj.radii(2);
            dist = vecnorm(q,2,2) - r ;
            [lateral_in,lateral_on, lateral_outer, lateral_out_near] = create_mask(lateral_mask,dist,r_min,eps);
    
            if obj.include_caps
                bottom_mask = z < 0;
                dist = (vecnorm(p - obj.centres(1,:),2,2)  )- obj.radii(1);
                [bottom_in,bottom_on, bottom_outer, bottom_out_near] = create_mask(bottom_mask,dist,r_min,eps);
                
                top_mask = z > obj.length;
                dist = (vecnorm(p - obj.centres(2,:),2,2)  ) - obj.radii(2);
                [top_in,top_on, top_outer, top_out_near] = create_mask(top_mask,dist,r_min,eps);
    
                in = bottom_in | top_in | lateral_in;
                on = bottom_on | top_on | lateral_on;
                outer = bottom_outer | top_outer | lateral_outer;
                out_near = bottom_out_near | top_out_near | lateral_out_near;
            else
                % in = lateral_in;
                % on =  lateral_on;
                % outer = lateral_outer | not(lateral_mask); 
                % out_near =lateral_out_near;
                % Changed 21/02.
                bottom_mask = z < 0;
                dist = (vecnorm(p - obj.centres(1,:),2,2)  )- obj.radii(1);
                [bottom_in,bottom_on, bottom_outer, bottom_out_near] = create_mask(bottom_mask,dist,r_min,eps);
                top_mask = z > obj.length;
                in = bottom_in  | lateral_in;
                on = bottom_on | lateral_on;
                outer = bottom_outer  | lateral_outer | top_mask;
                out_near = bottom_out_near | lateral_out_near;
            end



            

           
            
        end
        
        
        
    end
        


end