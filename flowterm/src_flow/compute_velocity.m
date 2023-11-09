function velocity_new = compute_velocity(femesh, velocity, setup)
    % Create velocity field defined on each element
    ncompartment = femesh.ncompartment;
    velocity_new = cell(1, ncompartment);
    for icmpt = 1:ncompartment
        % Finite elements
        elements = femesh.elements{icmpt};
    
        velocity_new{icmpt} = zeros(3,size(elements,2));
        velocity_new{icmpt}(3,:) = ones(1,size(elements,2) );
        velocity_new{icmpt} = velocity{icmpt}*velocity_new{icmpt};
    end
end

