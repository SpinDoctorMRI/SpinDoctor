function inradii = computeDiameters_undir(points, elements)
    % Number of elements
    numElements = size(elements, 1);
    
    % Initialize radii to zeros
    inradii = zeros(numElements, 1);
    
    % Iterate over each element
    for i = 1:numElements
        % Extract the four vertices of the tetrahedron
        v1 = points(elements(i, 1), :);
        v2 = points(elements(i, 2), :);
        v3 = points(elements(i, 3), :);
        v4 = points(elements(i, 4), :);
        
        % Compute volume V of the tetrahedron
        V = abs(dot(v1 - v4, cross(v2 - v4, v3 - v4))) / 6;
        
        % Compute surface area S of the tetrahedron
        S = 0.5 * (norm(cross(v2-v1, v3-v1)) + norm(cross(v2-v1, v4-v1)) + ...
                   norm(cross(v2-v3, v4-v3)) + norm(cross(v1-v3, v4-v3)));
        
        % Compute Inradius
        inradii(i) = 6* V / S;
        
        % Compute the edge lengths
%         a = norm(v2 - v1);
%         b = norm(v3 - v1);
%         c = norm(v4 - v1);
%         D = norm(v2 - v3); % opposite edge
%         e = norm(v3 - v4);
%         f = norm(v4 - v2);
        
        % Cayley-Menger Determinant
%         CM = det([
%             0, 1, a^2, b^2, c^2, 1;
%             1, 0, D^2, e^2, f^2, 1;
%             a^2, D^2, 0, f^2, e^2, 1;
%             b^2, e^2, f^2, 0, D^2, 1;
%             c^2, f^2, e^2, D^2, 0, 1;
%             1, 1, 1, 1, 1, 0
%         ]);
        
        % Compute Circumradius
%         circumradii(i) = abs(a*b*c*D / (4 * sqrt(abs(CM))));
    end
end

