function [ode_func, Jacobian] = btpde_functions(K, Q, R, J, q, seq)
%BTPDE_FUNCTIONS Generate BTPDE ODE function and its state-Jacobian.
%
% The ODE is defined by
%
%   M * dy/dt = ode_func(t, y),
%
% and the Jacobian
%
%   d(ode_func)/dy = Jacobian(t, y).
%
% The matrix M refers to the mass matrix.
%
% Args:
%   K: [npoint x npoint] - Stiffness matrix
%   Q: [npoint x npoint] - Flux matrix
%   R: [npoint x npoint] - Relaxation matrix
%   J: [npoint x npoint] - Moment matrix in gradient direction
%   q: [1 x 1]           - Q-value
%   seq: Sequence        - Gradient sequence
%
%   ode_func: ODE function at time t and state y
%   Jacobian: Jacobian of ODE function with respect to the state y


ode_func = @(t, y) -(K * y + Q * y + R * y + 1i * q * seq.call(t) * (J * y));
Jacobian = @(t, y) -(K + Q + R + 1i * q * seq.call(t) * J);
