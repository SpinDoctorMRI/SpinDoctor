function signal = solve_analytical(setup, volumes)
%SOLVE_ANALYTICAL Compute signal for analytical cylinders or spheres.
%
%   setup: struct
%
%   signal: [namplitude x nsequence]
%
% Units in analytical module - all SI-units:
%   permeability: m / s
%   water diffusion coefficient: m^2 / s
%   sphere radius: m
%   gyromagnetic ratio: rad / T / s
%   gradient pulse duration: s
%   time between pulses: s


% Sizes
namplitude = length(setup.gradient.values);
nsequence = length(setup.gradient.sequences);

signal = zeros(namplitude, nsequence);

warning("solve_analytical is not available, returning zero");
