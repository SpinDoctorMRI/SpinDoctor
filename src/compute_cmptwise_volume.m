function volumes = compute_cmptwise_volume(volumes, compartments)
%COMPUTE_CMPTWISE_VOLUME Compute total volumes for in, out and ecs
%compartments.
%
%   volumes: double(1, ncompartment)
%   compartments: string(1, ncompartment)


volume_in = sum(volumes(compartments=="in"));
volume_out = sum(volumes(compartments=="out"));
volume_ecs = sum(volumes(compartments=="ecs"));

if volume_in == 0 && volume_ecs == 0
    volumes = volume_out;
elseif volume_in > 0 && volume_ecs == 0
    volumes = [volume_in volume_out];
elseif volume_in == 0 && volume_ecs > 0
    volumes = [volume_out volume_ecs];
elseif volume_in > 0 && volume_ecs > 0
    volumes = [volume_in volume_out volume_ecs];
end

end
