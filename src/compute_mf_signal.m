function [signal_allcmpts, itertimes] = compute_mf_signal(setup, initial_signal, lap_eig, neig_lap)
%COMPUTE_MF_SIGNAL Compute the Matrix Formalism total signal.
%   This function only works under the assumption of a constant intial spin
%   density.
%
%   setup: struct
%   initial_signal
%   lap_eig
%   neig_lap (optional)
%
%   signal_allcmpts
%   itertimes


% Time function evaluation
starttime = tic;

do_reduce = nargin == nargin(@compute_mf_signal);

% Extract experiment parameters
sequences = setup.gradient.sequences;

% Extract gradient directions
dir_points = setup.gradient.directions.points;
dir_inds = setup.gradient.directions.indices;
dir_opposite = setup.gradient.directions.opposite;

lap_eig_vals = lap_eig.values;

% Sizes
nsequence = length(sequences);
namplitude = length(setup.gradient.values);
ndirection = setup.gradient.ndirection;
ninterval = setup.mf.ninterval;
neig = length(lap_eig_vals);

% Initialize output arguments
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

% Compute signal

if do_reduce && neig_lap < neig
    inds = 1:neig_lap;
else
    inds = 1:neig;
end
lap_eig_vals = lap_eig_vals(inds);

L_mat = diag(lap_eig_vals);
A_mat = lap_eig.moments(inds, inds, :);
for idir = dir_inds
    g = dir_points(:, idir);

    for iseq = 1:nsequence
        % Experiment parameters
        seq = sequences{iseq};

        % Create time intervals for time profile approximation
        time = linspace(0, seq.echotime, ninterval + 1);

        for iamp = 1:namplitude
            q = setup.gradient.qvalues(iamp, iseq);
            b = setup.gradient.bvalues(iamp, iseq);

            fprintf("Computing MF signal for:\n" ...
                + "  Direction %d of %d: g = [%.2f; %.2f; %.2f]\n" ...
                + "  Sequence  %d of %d: f = %s\n" ...
                + "  Amplitude %d of %d: q = %g, b = %g\n", ...
                idir, ndirection, g, ...
                iseq, nsequence, seq, ...
                iamp, namplitude, q, b);

            % Measure time
            itertime = tic;

            % Compute gdir contribution to BT operator
            W_mat = sum(A_mat .* shiftdim(g, -2), 3);

            % If PGSE, use three intervals, otherwise all intervals in time
            if isa(seq, "PGSE")
                K = L_mat + 1i * q * W_mat;
                edK = expm(-seq.delta * K);
                edL = exp(-(seq.Delta - seq.delta) * lap_eig_vals);
                H11 = edK(:, 1)' * (edL .* edK(:, 1));
            else
                % FE discretized BT operator for a given time profile
                K = @(ft) L_mat + sqrt(-1) * q * ft * W_mat;

                % Compute the first column of the matrix H using a
                % piecewise constant approximation of the time profile
                H1 = eye(neig, 1);
                for iint = 1:ninterval
                    dt = time(iint + 1) - time(iint);
                    ft = 1 / 2 * (seq.call(time(iint + 1)) + seq.call(time(iint)));
                    H1 = expm(-dt * K(ft)) * H1;
                end
                H11 = H1(1);
            end

            % Signal is obtained by multiplying the initial signal with H11
            signal_allcmpts(iamp, iseq, idir) = H11 * initial_signal;

            % Store computational time
            itertimes(iamp, iseq, idir) = toc(itertime);
        end
    end
end


for idir = dir_inds
    if ~isempty(dir_opposite{idir})
        % Opposite direction
        signal_allcmpts(:, :, dir_opposite{idir}) = conj(signal_allcmpts(:, :, idir));
    end
end

% Display elapsed time of function evaluation
toc(starttime);
