function tensors = fit_tensor(directions, adc)
%FIT_TENSOR Fit effective diffusion tensors to directionalized ADCs.
%   The six components of the symmetric diffusion tensors are fitted by least
%   squares to the gradient directions and resulting ADCs.
%
%   directions: [3 x ndirection]
%   adc: [ncompartment x ndirection]
%
%   tensors: [3 x 3 x ncompartment]


% Sizes
[ncompartment, ndirection] = size(adc);

% Get vectors of gradient interactions (basis vectors for symmetrical diffusion
% tensors)
g = directions;
G = zeros(ndirection, 6);
G(:, 1) = g(1, :).^2;
G(:, 2) = g(2, :).^2;
G(:, 3) = g(3, :).^2;
G(:, 4) = g(1, :) .*  g(2, :);
G(:, 5) = g(1, :) .*  g(3, :);
G(:, 6) = g(2, :) .*  g(3, :);

% Deduce vector of effective diffusion tensor components in each compartment by
% least squares fit of the directions to the computed ADCs
Dvec = G \ adc';

% Tensorize component vectors (symmetric)
tensors = zeros(3, 3, ncompartment);
tensors(1, 1, :) = Dvec(1, :);
tensors(2, 2, :) = Dvec(2, :);
tensors(3, 3, :) = Dvec(3, :);
tensors(1, 2, :) = Dvec(4, :);
tensors(2, 1, :) = Dvec(4, :);
tensors(1, 3, :) = Dvec(5, :);
tensors(3, 1, :) = Dvec(5, :);
tensors(2, 3, :) = Dvec(6, :);
tensors(3, 2, :) = Dvec(6, :);

