function mf_jn = compute_mf_jn(eigvals, diffusivity, setup)
%COMPUTE_MF_JN Compute the quantity J(lambda_n, f).
%
%   eigvals
%   diffusivity
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
    d = seq.delta;
    D = seq.Delta;

    % Time
    ntime = ninterval + 1;
    echotime = seq.echotime;
    time = linspace(0, echotime, ntime);
    dtime = echotime / ntime;

    for ieig = 1:neig
        % c = diffusivity * eigvals(ieig);
        c = eigvals(ieig);
        if abs(c) < 1e-15
            tmp = 0;
        elseif isa(seq, "PGSE")
            tmp = -1 / c^2 * ( ...
                + exp(-c * (D + d)) ...
                + exp(-c * (D - d)) ...
                - 2 * exp(-c * d)...
                - 2 * exp(-c * D)...
                + 2 * (1 - c * d));
        elseif isa(seq, "SinOGSE")
            n = seq.nperiod;
            a = 2 * pi * n / d;
            tmp = a^2 / (a^2+c^2)^2 * ( ...
                + 2 ...
                + d * c...
                + d * c^3 / a^2 ...
                + exp(-c * (D - d)) ...
                - 2 * exp(-c * D) ...
                - 2 * exp(-c * d) ...
                + exp(-c * (D + d)));
        elseif isa(seq, "CosOGSE")
            n = seq.nperiod;
            a = 2 * pi * n / d;
            ed = exp(-c * d);
            eD = exp(-c * D);
            em = exp(-c * (D - d));
            ep = exp(-c * (D + d));
            w2 = a * c^2;
            tmp = (w2 * (2 * eD - em + 2 * ed - ep) ...
                + c * a * (a^2 * d + c^2 * d - 2 * c)) ...
                / (a * (a^2 + c^2)^2);
            % elseif isa(seq, "DoublePGSE")
            %     tmp = TODO;
        else
            fprintf("Eigenvalue %d of %d: Numerical integration\n", ieig, neig);
            tmp = @(t) integral(@(s) exp(c * (s - t)) .* seq.call(s), 0, t, ...
                "Waypoints", [seq.delta, seq.Delta]);% , "AbsTol", 1e-6, "RelTol", 1e-3);
            tmp = c * seq.integral(time) .* arrayfun(tmp, time);
            tmp = dtime * trapz(tmp);
        end
        mf_jn(iseq, ieig) = tmp / seq.bvalue_no_q / diffusivity;
    end
end
