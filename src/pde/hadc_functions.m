function [ode_func, Jacobian] = hadc_functions(K, surfint, sequence)
%HADC_FUNCTIONS Generate HADC ODE function and its state-Jacobian.
%
% The ODE is defined by
%
%     M * dy/dt = ode_func(t, y),
%
% and the Jacobian
%
%     d(ode_func)/dy = Jacobian(t, y).
%
% The matrix M refers to the mass matrix.
%
%   K: double(nnode, nnode)
%       Stiffness matrix
%   surfint: double(nnode, 1)
%       Boundary integral of normal component of gradient direction
%	sequence
%
%   ode_func: ODE function at time t and state y
%  	Jacobian: Jacobian of ODE function with respect to the state y


ode_func = @(t, y) -(K * y) + surfint .* sequence.integral(t);
Jacobian = -K;
