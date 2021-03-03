function [signal_allcmpts, itertimes] = compute_mf_signal_sub_bt(experiment, initial_signal, lap_eig, bt_eig, neig_bt)
%COMPUTE_MF_SIGNAL_SUB_BT Compute the explicit Matrix Formalism.
%   This only works for PGSE, and requires a decomposed BT-matrix
%
%   experiment: struct with fields
%       ndirection
%     	directions
%    	qvalues
%    	bvalues
%     	sequences
% 	initial_signal
% 	lap_eig
%  	bt_eig (optional)
%
%   signal_allcmpts
%  	itertimes


starttime = tic;

do_reduce = nargin == nargin(@compute_mf_signal_sub_bt);

% Extract experiment parameters
sequences = experiment.sequences;

% Extract gradient directions
dir_points = experiment.directions.points;
dir_inds = experiment.directions.indices;
dir_opposite = experiment.directions.opposite;

% Sizes
nsequence = length(sequences);
namplitude = size(experiment.bvalues, 2);
ndirection = size(dir_points, 1);
neig = length(lap_eig.values);

% Initialize output arguments
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

for iseq = 1:nsequence
    % Experiment parameters
    seq = sequences{iseq};

    for iamp = 1:namplitude
        q = experiment.qvalues(iamp, iseq);
        b = experiment.bvalues(iamp, iseq);
        for idir = dir_inds
            fprintf("Computing MF signal for:\n" ...
                + "  Direction %d of %d: [%.2f; %.2f; %.2f]\n" ...
                + "  Time profile %d of %d: %s\n" ...
                + "  Amplitude %d of %d: q-value %g, b-value %g\n", ...
                idir, ndirection, dir_points(:, idir), ...
                iseq, nsequence, seq, ...
                iamp, namplitude, q, b);

            % Measure time
            itertime = tic;

            if do_reduce && neig_bt < neig
                ind = 1:neig_bt;
            else
                ind = 1:neig;
            end

            % inds_sort = 1:neig;
            % [~, inds_sort] = sort(abs(bt_eig.invVsortC1{iamp, iseq, idir}.'), "descend");
            % [~, inds_sort] = sort(abs(bt_eig.Vsort{iamp, iseq, idir}(1, :)), "descend");
            % ind = inds_sort(ind);

            % Load BT eigendecomposition
            V = bt_eig.Vsort{iamp, iseq, idir}(:, ind);
            D = bt_eig.Dsort{iamp, iseq, idir}(ind);
            % invV = bt_eig.invVsort{iamp, iseq, idir}(ind, :);
            invVC1 = bt_eig.invVsortC1{iamp, iseq, idir}(ind);

            H11 = invVC1' * diag(exp(-seq.delta * D))' * V'...
                * diag(exp(-(seq.Delta - seq.delta) * lap_eig.values))...
                * V * diag(exp(-seq.delta * D)) * invVC1;


            % Signal is obtained by multiplying the initial signal with H11
            signal_allcmpts(iamp, iseq, idir) = H11 * initial_signal;
            if ~isempty(dir_opposite{idir})
                % Opposite direction
                signal_allcmpts(iamp, iseq, dir_opposite{idir}) ...
                    = conj(signal_allcmpts(iamp, iseq, idir));
            end

            % Store computational time
            itertimes(iamp, iseq, idir) = toc(itertime);
        end
    end
end

toc(starttime);
