classdef Swc_point_cloud
    %SWC_POINT_CLOUD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        swc % swc file to store point and connectivity data
        somas % array of soma objects (point cloud and geometrical data)
        segments % array of segment objects (point cloud and geometrical data)
        density % density of point cloud data
        intersection_matrix
    end
    
    methods
        function swc_p = Swc_point_cloud(swc)
            swc_p.swc  = Swc(swc);
            [somas,segments] =  make_segment_list(swc_p.swc.point_data,swc_p.swc.connectivity_data,swc_p.swc.type_data,1:size(swc_p.swc.point_data,1) ,true); % Changed true to false 20/02
            swc_p.somas = somas;
            swc_p.segments = segments;            
            
        end
        % function  swc_p = create_all_point_clouds(swc_p,density)
        %     if  nargin == 1
        %         swc_p.density = 1;      
        %     end
        %     [somas,segments] =  make_pointclouds(swc_p.somas,swc_p.segments,swc_p.density,false);                     
        %     swc_p.somas = somas;
        %     swc_p.segments = segments;
        % end
        % function [points,normals] = extract_point_cloud(swc_p)
        %     somas = swc_p.somas; segments = swc_p.segments;
        %     N = size(somas,1);M = size(segments,1); total = N+M;
        %     points = [];normals = [];
        %     for i =1:N
        %         [p,n] = somas(i).output_filtered_points();
        %         if length(p) > 0
        %             points= [points;p];
        %             normals = [normals;n];
        %         end
        %     end
        %     for i =1:M
        %         [p,n] = segments(i).output_filtered_points();
        %         if length(p) > 0
        %             points= [points;p];
        %             normals = [normals;n];
        %         end
        %     end
        % end
        
        % function plot(swc_p)
        %     [p,~] = swc_p.extract_point_cloud();
        %     scatter3(p(:,1),p(:,2),p(:,3)),
        % end
        function N = total_points(swc_p)
            N = 0;
            M = length(swc_p.somas);
            for i = 1:M
                N = N + length(swc_p.somas(i).points);
            end
            M = length(swc_p.segments);
            for i = 1:M
                N = N + length(swc_p.segments(i).points);
            end
        end
        function area = find_min_area(swc_p)
            M = length(swc_p.somas);
            area = Inf;
            for i = 1:M
                area = min(swc_p.somas(i).area,area);
            end
            M = length(swc_p.segments);
            for i = 1:M
                area = min(swc_p.segments(i).area,area);            
            end
        end    
    end
end

