function eigenvalues = length2eig(length_scales, diffusivity)
%LENGTH2EIG Convert length scales into Laplace eigenvalues.
%
%   length_scales: double(neig, 1)
%   diffusivity: double
%
%   eigenvalues: double(neig, 1)


eigenvalues = diffusivity * pi^2 ./ length_scales^2;
