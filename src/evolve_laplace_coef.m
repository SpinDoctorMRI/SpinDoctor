function nu = evolve_laplace_coef(nu0, seq, iqA, LT2, ninterval)
%EVOLVE_LAPLACE_COEF Compute Laplace coefficients of final magnetization.
%   If PGSE, use three intervals; 
%   if DoublePGSE, use seven intervals; 
%   otherwise many intervals. TO DO: optimization for OGSE
%
%   nu0: double(neig, 1)
%   seq: Sequence
%   iqA: complex(neig, neig)
%   LT2: double(neig, neig)
%   ninterval: int   
%
%   nu: double(neig, 1)


% BT operator in Laplace basis
K = LT2 + iqA;

if isa(seq, "PGSE")
    % Laplace coefficients of final magnetization
    nu = expmv(-seq.delta, K, nu0);
    nu = expm(-(seq.Delta - seq.delta) * LT2) * nu;
    nu = expmv(-seq.delta, K', nu);

    % edK = expm(-seq.delta * K);
    % edL = expm(-(seq.Delta - seq.delta) * LT2);
    % nu = edK' * (edL * (edK * nu0));
elseif isa(seq, "DoublePGSE")
    % Laplace coefficients of final magnetization
    edK = expm(-seq.delta * K);
    edL = expm(-(seq.Delta - seq.delta) * LT2);
    etL = expm(-seq.tpause * LT2);
    nu = edK' * (edL * (edK * (etL * (edK' * (edL * (edK * nu0))))));

    % nu = expmv(-seq.delta, K, nu0);
    % nu = expm(-(seq.Delta - seq.delta) * LT2) * nu;
    % nu = expmv(-seq.delta, K', nu);
    % nu = expm(-seq.tpause * LT2) * nu;
    % nu = expmv(-seq.delta, K, nu);
    % nu = expm(-(seq.Delta - seq.delta) * LT2) * nu;
    % nu = expmv(-seq.delta, K', nu);
elseif isa(seq, "SinOGSE") || isa(seq, "CosOGSE")
    % Transform Laplace coefficients using piecewise constant
    % approximation of time profile    
    
    % ninterval is the number of intervals in one period
    n_quater_periods = 4 * seq.nperiod;
    quater_period = seq.delta / n_quater_periods;
    nsubinterval = round(ninterval / n_quater_periods);
    % at least 3 intervals in a quater period
    nsubinterval = max([3, nsubinterval]);

    % Create time intervals for a quater period
    time = linspace(0, quater_period, nsubinterval+1);
    % BT operator in Laplace basis for a given time profile value
    K = @(ft) LT2 + ft * iqA;

    nu = nu0;
    % first pulse
    for iquater = 1:n_quater_periods
        itime = time + (iquater-1) * quater_period;
        nu = riemann_approximation(nu, K, itime, seq);
    end
    % pause
    nu = expm(-(seq.Delta - seq.delta) * LT2) * nu;
    % second pulse
    for iquater = 1:n_quater_periods
        itime = time + seq.Delta + (iquater-1) * quater_period;
        nu = riemann_approximation(nu, K, itime, seq);
    end
else
    % Transform Laplace coefficients using piecewise constant
    % approximation of time profile
    
    % Create time intervals for time profile approximation
    time = linspace(0, seq.echotime, ninterval + 1);
    % BT operator in Laplace basis for a given time profile value
    K = @(ft) LT2 + ft * iqA;
    nu = riemann_approximation(nu0, K, time, seq);
end
end


function nu = riemann_approximation(nu, K, time, seq)
    ninterval = length(time) - 1;
    for i = 1:ninterval
        % Time step and time profile on given interval
        dt = time(i + 1) - time(i);
        ft = (seq.call(time(i + 1)) + seq.call(time(i))) / 2;

        % Laplace coefficients of magnetization at end of interval
        nu = expmv(-dt, K(ft), nu);
        % nu = expm(-dt * K(ft)) * nu;
    end
end
