function mean_diffusivity = compute_mean_diffusivity(diffusivity, femesh)
%COMPUTE_MEAN_DIFFUSIVITY Compute volume weighted mean of diffusivities over compartments.
%   Take trace of each diffusion tensor, divide by 3.
%
%   diffusivity: double(3, 3, ncompartment)
%   femesh: struct
%
%   mean_diffusivity: double

mean_diffusivity = trace(sum(diffusivity .* ...
    shiftdim(femesh.volumes, -1), 3)) / (3 * sum(femesh.total_volume));
end
