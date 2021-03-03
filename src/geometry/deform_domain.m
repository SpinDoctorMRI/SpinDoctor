function nodes = deform_domain(nodes, deformation)
%DEFORM_DOMAIN Deform domain.


height =  max(nodes(3, :)) -  min(nodes(3, :));
width = max(nodes(1, :)) -  min(nodes(1, :));

bend = deformation(1);
twist = deformation(2);

thvec = nodes(3, :) / height * twist;


nodes(1, :) = nodes(1, :) + bend * 30 * width * (nodes(3, :) / height).^2;

center = 1 / 2 * (max(nodes(1:2, :), [], 2) + min(nodes(1:2, :), [], 2));


% rmax = max(vecnorm(nodes(1:2, :)));

% nodes(1, :) = nodes(1, :) + rmax;
nodes(1:2, :) = nodes(1:2, :) - center;

nodes(1:2, :) = [cos(thvec) .* nodes(1, :) - sin(thvec) .* nodes(2, :);
    sin(thvec) .* nodes(1, :) + cos(thvec) .* nodes(2, :)];

% nodes(1, :) = nodes(1, :) - rmax;
nodes(1:2, :) = nodes(1:2, :) + center;
