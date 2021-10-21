function mf_jn = compute_mf_jn(eigvals, setup)
%COMPUTE_MF_JN Compute the quantity J(lambda_n, f).
%
%   eigvals
%   setup: struct
%
%   mf_jn


disp("Computing the quantity J(lambda_n, f)");

% Extract experiment parameters
sequences = setup.gradient.sequences;
ninterval = setup.mf.ninterval;

% Sizes
nsequence = length(sequences);
neig = length(eigvals);

% Initialize output argument
mf_jn = zeros(nsequence, neig);

% Compute Jn
for iseq = 1:nsequence
    fprintf("  Experiment %d of %d\n", iseq, nsequence);

    % Gradient sequences
    seq = sequences{iseq};

    % Time
    ntime = ninterval + 1;
    echotime = seq.echotime;
    time = linspace(0, echotime, ntime);
    dtime = echotime / ntime;

    for ieig = 1:neig
        lambda = eigvals(ieig);
        if abs(lambda) < 1e-16
            tmp = 0;
        elseif isa(seq, "PGSE")
            tmp = J_PGSE(lambda, seq);
        elseif isa(seq, "DoublePGSE")
            tmp = J_DoublePGSE(lambda, seq);
        elseif isa(seq, "SinOGSE")
            tmp = J_SinOGSE(lambda, seq);
        elseif isa(seq, "CosOGSE")
            tmp = J_CosOGSE(lambda, seq);
        else
            fprintf("Eigenvalue %d of %d: Numerical integration\n", ieig, neig);
            tmp = @(t) integral(@(s) exp(lambda * (s - t)) .* seq.call(s), 0, t, ...
                "Waypoints", [seq.delta, seq.Delta]);% , "AbsTol", 1e-6, "RelTol", 1e-3);
            tmp = lambda * seq.integral(time) .* arrayfun(tmp, time);
            tmp = dtime * trapz(tmp) / seq.bvalue_no_q;
        end
        mf_jn(iseq, ieig) = tmp;
    end
end
