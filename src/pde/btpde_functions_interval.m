function [ode_func, Jacobian] = btpde_functions_interval(K, Q, J, q, sequence, interval_midpoint)
%BTPDE_FUNCTIONS_INTERVAL Generate BTPDE ODE function and its state-Jacobian.
%   This function gives the ODE function for a given interval only, determined
%   by its midpoint.
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
%     Q: Flux matrix
%     J: Mass matrix weighted by dot(g, x), where g is the gradient direction
%           and x the spatial variable
%      q: q-value
%    sequence
%     interval_midpoint
%
%   ode_func: ODE function at time t and state y
%      Jacobian: Jacobian of ODE function with respect to the state y


% Check for possible simplifications
if isa(sequence, "PGSE") || isa(sequence, "DoublePGSE")
        % For (double) PGSE, the time profile is constant on each interval
        constant_timeprofile = true;
elseif isa(sequence, "CosOGSE") || isa(sequence, "SinOGSE")
        if interval_midpoint < sequence.delta || sequence.Delta <= interval_midpoint
            % The time profile is time dependent during the pulses
            % Use full form
            constant_timeprofile = false;
        else
            % The time profile is zero on this interval
            constant_timeprofile = true;
        end
else
    % Use full form
    constant_timeprofile = false;
end

if constant_timeprofile
    % The time profile is constant on this interval. The Jacobian function is
    % thus a constant matrix, and the ODE function is linear in the state y
    Jacobian = -(K + Q + 1i * q * sequence.call(interval_midpoint) * J);
    ode_func = @(t, y) Jacobian * y;
else
    % The time profile is not constant on this interval, use full form
    ode_func = @(t, y) -(K * y + Q * y + 1i * q * sequence.call(t) * (J * y));
    Jacobian = @(t, y) -(K + Q + 1i * q * sequence.call(t) * J);
end
