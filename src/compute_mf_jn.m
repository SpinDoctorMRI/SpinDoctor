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

    for ieig = 1:neig
        lambda = eigvals(ieig);
        if abs(lambda) < 1e-16
            tmp = 0;
        elseif ismember(class(seq), ["PGSE", "DoublePGSE", "CosOGSE", "SinOGSE"])
            tmp = seq.J(lambda);
        else
            tmp = seq.J(lambda, ninterval);
        end
        mf_jn(iseq, ieig) = tmp;
    end
end
