function [ode_func, Jacobian] = btpde_functions(K, Q, J, qvalue, sequence)
%BTPDE_FUNCTIONS Generate BTPDE ODE function and its state-Jacobian.
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
%   K: Stiffness matrix
% 	Q: Flux matrix
% 	J: Mass matrix weighted by dot(g, x), where g is the gradient direction
%           and x the spatial variable
%  	qvalue
%	sequence
%
%   ode_func: ODE function at time t and state y
%  	Jacobian: Jacobian of ODE function with respect to the state y


ode_func = @(t, y) -(K * y + Q * y + 1i * qvalue * sequence.call(t) * (J * y));
Jacobian = @(t, y) -(K + Q + 1i * qvalue * sequence.call(t) * J);
