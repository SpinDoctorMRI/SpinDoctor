function h_values = computeDiameters_dir(points ,elements, velocity)
    % elements: Nx4 matrix where each row contains indices of the tetrahedron's vertices
    % points: Mx3 matrix where each row is a point (x, y, z)
    % velocity: 1x3 vector representing the velocity
    num_elements = size(elements, 1);
    h_values = zeros(num_elements, 1);

%     for e = 1:num_elements
%         % Extract the vertices of the current tetrahedron
%         vertices = points(elements(e,:), :);
% 
%         % Compute the volume of the tetrahedron
%         K = det([vertices(2,:) - vertices(1,:); vertices(3,:) - vertices(1,:); vertices(4,:) - vertices(1,:)]) / 6;
% 
%         % Compute the areas and dot products
%         areas = zeros(4,1);
%         dots = zeros(4,1);
%         faces = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
%         for i = 1:4
%             face = vertices(faces(i,:), :);
%             normal = cross(face(2,:) - face(1,:), face(3,:) - face(1,:));
%             areas(i) = 0.5 * norm(normal);
%             normal = normal / norm(normal);
%             v_norm = velocity / norm(velocity);
%             dots(i) = abs(dot(normal, v_norm));
%         end
% 
%         % Find the area of the face most orthogonal to the velocity
%         [~, idx] = min(dots);
%         area_Fv = areas(idx);
% 
%         % Compute h for this tetrahedron
%         h_values(e) = K / area_Fv;
%     end
    
    % ... (as before)
    v_norm = velocity / norm(velocity);  % normalized velocity

    for e = 1:num_elements
        % Extract the vertices of the current tetrahedron
        vertices = points(elements(e,:), :);
        % ... (as before)
        
        % Compute the volume of the tetrahedron
        K = det([vertices(2,:) - vertices(1,:); vertices(3,:) - vertices(1,:); vertices(4,:) - vertices(1,:)]) / 6;


        % Compute projected areas and dot products
        proj_areas = zeros(4,1);
        faces = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
%         dots = zeros(4,1);
        for i = 1:4
            face = vertices(faces(i,:), :);
            normal = cross(face(2,:) - face(1,:), face(3,:) - face(1,:));
            area = 0.5 * norm(normal);
            normal = normal / norm(normal);
            cosine_theta = abs(dot(normal, v_norm));
            proj_areas(i) = area * cosine_theta;
%             dots(i) = abs(dot(normal, v_norm));
        end

        % Find the projected area of the face most orthogonal to the velocity
%         [~, idx] = min(dots);
%         proj_area_Fv = proj_areas(idx);
        proj_area_Fv = 2*mean(proj_areas);

        % Compute h for this tetrahedron
        h_values(e) = 3*K / proj_area_Fv;
    end
end


