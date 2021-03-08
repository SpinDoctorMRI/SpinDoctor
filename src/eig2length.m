function length_scales = eig2length(eigenvalues, diffusivity)
%EIG2LENGTH Convert Laplace eigenvalues into length scale.
%
%   eigenvalues: double(neig, 1)
%   diffusivity: double
%
%   length_scales: double(neig, 1)


length_scales = pi * sqrt(diffusivity ./ max(0, eigenvalues));
