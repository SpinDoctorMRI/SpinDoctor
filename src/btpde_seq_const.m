function mag = btpde_seq_const(seq,ug,q,b,g,rho,Jx,K,Q,R,solve_ode,options_template,idir,iseq,iamp,ndirection,nsequence,namplitude,npoint_cmpts,solver_str)
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
    J = ug(1) * Jx{1} + ug(2) * Jx{2} + ug(3) * Jx{3};
    
    % Initial magnetization
    mag = rho;
    
    % Solve for each interval consecutively
    for iint = 1:ninterval
    
        % Add a third point to the interval, so that the ODE solver does not
        % store the magnetization for all time steps during the solve. If
        % there were only two points in the interval, the ODE solver would
        % store all time steps. This would require a lot of memory,
        % especially during parfor iterations
        interval_midpoint = (timelist(iint) + timelist(iint + 1)) / 2;
        time_list_interval = [timelist(iint), interval_midpoint, timelist(iint + 1)];
    
        % Display state of iterations
        fprintf( ...
            join([
                "Solving BTPDE of size %d using %s:"
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
    
        % Create new ODE functions on given interval
        [ode_function, Jacobian] = btpde_functions_interval( ...
            K, Q, R, J, q, seq, interval_midpoint);
    
        % Update options with new Jacobian, which is either a
        % function handle or a constant matrix, depending on the
        % time profile
        options = odeset(options_template, "Jacobian", Jacobian);
    
        % Solve ODE on domain, starting from the magnetization at
        % the end of the previous interval (mag)
        [~, y] = solve_ode(ode_function, time_list_interval, mag, options);
    
        % Magnetization at end of interval
        mag = y(end, :).';
    end
end