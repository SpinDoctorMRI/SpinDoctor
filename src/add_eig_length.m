function lap_eig = add_eig_length(lap_eig, diffusivity)
%ADD_EIG_LENGTH Add eigen length scale to lap_eig.
%
%   lap_eig: struct
%   diffusivity: double
%
%   lap_eig: struct

for icmpt = 1:length(lap_eig)
    lap_eig(icmpt).length_scales = eig2length(lap_eig(icmpt).values, diffusivity);
end
