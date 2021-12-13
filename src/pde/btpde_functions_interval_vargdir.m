function [ode_func, Jacobian] = btpde_functions_interval_vargdir(K, Q, R, J, q, seq, interval_midpoint)
%BTPDE_FUNCTIONS_INTERVAL Generate BTPDE ODE function and its state-Jacobian.
%   This function gives the ODE function for a given interval only, determined
%   by its midpoint.
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
%   K: [npoint x npoint]      - Stiffness matrix
%   Q: [npoint x npoint]      - Flux matrix
%   R: [npoint x npoint]      - Relaxation matrix
%   J: [npoint x npoint]      - Moment matrix in gradient direction
%   q: [1 x 1]                - Q-value
%   seq: Sequence             - Gradient sequence
%   interval_midpoint [1 x 1] - Midpoint of interval
%
%   ode_func: ODE function at time t and state y
%   Jacobian: Jacobian of ODE function with respect to the state y


% Check for possible simplifications
if isa(seq, "PGSE") || isa(seq, "DoublePGSE")
        % For (double) PGSE, the time profile is constant on each interval
        constant_timeprofile = true;
elseif isa(seq, "CosOGSE") || isa(seq, "SinOGSE")
        if interval_midpoint < seq.delta || seq.Delta <= interval_midpoint
            % The time profile is time dependent during the pulses
            % Use full form
            constant_timeprofile = false;
        else
            % The time profile is zero on this interval
            constant_timeprofile = true;
        end
elseif isa(seq, "VarGdir")
    constant_timeprofile = true;
else
    % Use full form
    
    %%%%%%%%%% For CoreyBaron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set this to true because within each defined time step, the
    %%% time profile is constant.
    %  constant_timeprofile = false;
    constant_timeprofile = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

if constant_timeprofile
    % The time profile is constant on this interval. The Jacobian function is
    % thus a constant matrix, and the ODE function is linear in the state y
    Jacobian = -(K + Q + R + 1i * q * seq.call(interval_midpoint) * J);
    ode_func = @(t, y) Jacobian * y;
else
    % The time profile is not constant on this interval, use full form
    ode_func = @(t, y) -(K * y + Q * y + R * y + 1i * q * seq.call(t) * (J * y));
    Jacobian = @(t, y) -(K + Q + R + 1i * q * seq.call(t) * J);
end
