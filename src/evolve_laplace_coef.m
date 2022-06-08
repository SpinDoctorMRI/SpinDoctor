function nu = evolve_laplace_coef(nu, seq, K, LT2, ninterval)
%EVOLVE_LAPLACE_COEF Compute Laplace coefficients of final magnetization.
%   If PGSE, use three intervals; 
%   if DoublePGSE, use seven intervals; 
%   otherwise many intervals.
%
%   nu: double(neig, 1)
%   seq: Sequence
%   K: complex(neig, neig) or function handle
%   LT2: double(neig, neig)
%   ninterval: int
%
%   nu: double(neig, 1)


LT2_is_matrix = size(LT2, 1) ~= 1;

if isa(seq, "PGSE")
    % Assemble matrix K
    if LT2_is_matrix
        K = K + LT2;
    else
        s = size(K);
        diag_ind = 1:s(1)+1:s(1)*s(2);
        K(diag_ind) = K(diag_ind) + LT2;
    end
    K = -seq.delta*K;

    % Compute Laplace coefficients of final magnetization
    [cost, K_shift, s, m, mu, tol] = expmv_cost(K, nu);
    if cost > 1e4
        clear K_shift;
        % first pulse
        K = expm(K);
        nu = K * nu;
        % between pulses
        t = seq.Delta - seq.delta;
        nu = between_pulses(t, nu, LT2, LT2_is_matrix);
        % second pulse
        K = conj(K);
        nu = K * nu;
    else
        clear K;
        % first pulse
        nu = expmv(K_shift, nu, s, m, mu, tol);
        % between pulses
        t = seq.Delta - seq.delta;
        nu = between_pulses(t, nu, LT2, LT2_is_matrix);
        % second pulse
        K_shift = conj(K_shift);
        nu = expmv(K_shift, nu, s, m, conj(mu), tol);
    end

elseif isa(seq, "DoublePGSE")
    if LT2_is_matrix
        K = K + LT2;
    else
        s = size(K);
        diag_ind = 1:s(1)+1:s(1)*s(2);
        K(diag_ind) = K(diag_ind) + LT2;
    end
    
    % Compute Laplace coefficients of final magnetization
    % % first pulse
    % K = expm(-seq.delta * K);
    % nu = K * nu;
    % % between pulses
    % t = seq.Delta - seq.delta;
    % nu = between_pulses(t, nu, LT2, LT2_is_matrix);
    % % second pulse
    % K = conj(K);
    % nu = K * nu;
    % % between pulses
    % t = seq.tpause;
    % nu = between_pulses(t, nu, LT2, LT2_is_matrix);
    % % third pulse
    % K = conj(K);
    % nu = K * nu;
    % % between pulses
    % t = seq.Delta - seq.delta;
    % nu = between_pulses(t, nu, LT2, LT2_is_matrix);
    % % fourth pulse
    % K = conj(K);
    % nu = K * nu;

    % first pulse
    nu = expmv(-seq.delta, K, nu);
    % between pulses
    t = seq.Delta - seq.delta;
    nu = between_pulses(t, nu, LT2, LT2_is_matrix);
    % second pulse
    K = conj(K);
    nu = expmv(-seq.delta, K, nu);
    % between pulses
    t = seq.tpause;
    nu = between_pulses(t, nu, LT2, LT2_is_matrix);
    % third pulse
    K = conj(K);
    nu = expmv(-seq.delta, K, nu);
    % between pulses
    t = seq.Delta - seq.delta;
    nu = between_pulses(t, nu, LT2, LT2_is_matrix);
    % fourth pulse
    K = conj(K);
    nu = expmv(-seq.delta, K, nu);

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

    % first pulse
    for iquater = 1:n_quater_periods
        itime = time + (iquater-1) * quater_period;
        nu = riemann_approximation(nu, @get_K, itime, seq);
    end
    % between pulses
    t = seq.Delta - seq.delta;
    nu = between_pulses(t, nu, LT2, LT2_is_matrix);
    % second pulse
    for iquater = 1:n_quater_periods
        itime = time + seq.Delta + (iquater-1) * quater_period;
        nu = riemann_approximation(nu, @get_K, itime, seq);
    end

else
    % Transform Laplace coefficients using piecewise constant
    % approximation of time profile

    % Create time intervals for time profile approximation
    time = linspace(0, seq.echotime, ninterval + 1);
    nu = riemann_approximation(nu, @get_K, time, seq);
end

function K_prime = get_K(ft)
    if LT2_is_matrix
        K_prime = K*ft + LT2;
    else
        K_prime = K*ft;
        s = size(K_prime);
        diag_ind = 1:s(1)+1:s(1)*s(2);
        K_prime(diag_ind) = K_prime(diag_ind) + LT2;
    end
end

end


function nu = between_pulses(t, nu, LT2, LT2_is_matrix)
    if abs(t) > 1e-16
        if LT2_is_matrix
            LT2 = -t*LT2;

            [cost, LT2_shift, s, m, mu, tol] = expmv_cost(LT2, nu);
            if cost > 1e4
                clear LT2_shift;
                nu = expm(LT2) * nu;
            else
                clear LT2;
                nu = expmv(LT2_shift, nu, s, m, mu, tol);
            end
        else
            nu = exp(-t * LT2') .* nu;
        end
    end
end


function nu = riemann_approximation(nu, get_K, time, seq)
    ninterval = length(time) - 1;
    for i = 1:ninterval
        % Time step and time profile on given interval
        dt = time(i + 1) - time(i);
        ft = (seq.call(time(i + 1)) + seq.call(time(i))) / 2;

        % Laplace coefficients of magnetization at end of interval
        % nu = expm(-dt * get_K(ft)) * nu;
        nu = expmv(-dt, get_K(ft), nu);
    end
end
