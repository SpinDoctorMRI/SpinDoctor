function mag = btpde_midpoint_seq_const(seq,ug,q,b,g,rho,dt,theta,M,Jx,K,Q,R,idir,iseq,iamp,ndirection,nsequence,namplitude,npoint_cmpts,solver_str)
%BTPDE_SEQ Solves the Bloch-Torrey partial differential equation for a
%single constant direction sequence
%   
%  This solves the Bloch-Torrey partial differential equation by a finite
%  element method for a single sequence which is of the form f(t)*ug
%
%   seq : (AbsSequence) a single sequence
%   ug: (array) unit direction vector
%   q : (float) q-value
%   b: (float) b-value    
%   g: (float) g-value
%   rho: (array) initial magnetisation
%   J: (matrix) gradient direction matrix
%   K : (matrix) stiffness matrix
%   Q : (matrix) flux matrix
%   R : (matrix) relaxation matrix
%   solve_ode: (func) ode solver
%   options_template: (odeset) parameters for ODE solver
%   idir: (int) direction index
%   iseq: (int) sequence index 
%   iamp: (int) amplitude index
%   ndirection: (int) number of sequences
%   nsequence: (int) number of sequences
%   namplitude: (int) number of amplitudes

% Get intervals based on the properties of the time profile
[timelist, interval_str, timeprofile_str] = seq.intervals;

% Number of intervals
ninterval = length(timelist) - 1;

% Assemble gradient direction dependent finite element matrix
A = ug(1) * Jx{1} + ug(2) * Jx{2} + ug(3) * Jx{3};

% Initial conditions
t = 0;
y = rho;

% Solve for each interval consecutively
for iint = 1:ninterval

    % Display state of iterations
    fprintf( ...
        join([
            "Solving BTPDE of size %d using %s"
            "  Direction %d of %d: ug = [%.2f; %.2f; %.2f]"
            "  Sequence  %d of %d: f = %s"
            "  Amplitude %d of %d: g = %g, q = %g, b = %g"
            "  Interval  %d of %d: I = %s, %s\n"
        ], newline), ...
        sum(npoint_cmpts), solver_str, ...
        idir, ndirection, ug, ...
        iseq, nsequence, seq, ...
        iamp, namplitude, g, q, b, ...
        iint, ninterval, interval_str(iint), timeprofile_str(iint) ...
    );
    
    % Midpoint of interval
    interval_midpoint = (timelist(iint) + timelist(iint + 1)) / 2;

    % Right hand side Jacobian
    J = -(K + Q + R + 1i * seq.call(interval_midpoint) * q * A);

    % Factorize left hand side matrix
    [L, U, P, QQ, D] = lu(M - dt * theta * J);

    % Right hand side matrix
    E = M + dt * (1 - theta) * J;

    % Time step loop
    while t + dt < timelist(iint + 1)
        % Advance by dt
        % t
        y = QQ * (U \ (L \ (P * (D \ (E * y)))));
        t = t + dt;
    end

    % Adapt time step to advance to end of interval. This requires a new
    % LU decomposition
    dt_last = timelist(iint + 1) - t;
    E = M + dt_last * (1 - theta) * J;
    [L, U, P, QQ, D] = lu(M - dt_last * theta * J);
    y = QQ * (U \ (L \ (P * (D \ (E * y)))));
    t = t + dt_last;
end
mag = y;
end