function [somas,segments]=  make_segment_list(point_data,radius_data,connectivity_data,type_data,indices,include_caps)
    if nargin ==3
        indices = 1:size(point_data,1);
        include_caps = true;
    elseif nargin ==4
        include_caps = true;
    end
    segments =[];
    somas = [];
    N = size(point_data,1);
    for i = 1:N
        centre = point_data(i,1:3);
        parent = connectivity_data(i,2);
        radius = radius_data(i);

        if (ismember(i,indices) || ismember(connectivity_data(i,2),indices) ) && (type_data(i) ~= 1 && connectivity_data(i,2) ~= -1)
            
            centres = [centre;point_data(parent,1:3)];
            radii= [radius,radius_data(parent)];
            
            if type_data(parent) == 1
                radii(2) = mean(radii);
            end

            segments = [segments; Segment(centres,radii,connectivity_data(i,:),type_data(i),include_caps)];            
        elseif ismember(i,indices) && (type_data(i) == 1 || connectivity_data(i,2) == -1)
            somas = [somas; Soma(centre,radius,connectivity_data(i,:),type_data(i))];
        end
    end
end