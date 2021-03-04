%% BT eigendecomposition

if isfield(setup, "mf")

    disp("Computing or loading the Bloch-Torrey eigenfunctions");

    % Create structure for storing results
    bt_eig.Vsort = cell(namplitude, nsequence, ndirection);
    bt_eig.Dsort = cell(namplitude, nsequence, ndirection);
    bt_eig.invVsort = cell(namplitude, nsequence, ndirection);
    bt_eig.invVsortC1 = cell(nsequence, namplitude, ndirection);
    bt_eig.totaltime = zeros(namplitude, nsequence);

    % Compute or load the Bloch-Torrey eigenfunctions
    for iseq = 1:nsequence
        for iamp = 1:namplitude
            q = setup.gradient.qvalues(iamp, iseq);

            % Make file name
            fname = sprintf("bt_eig_lengthscale%g_q%g.mat", setup.mf.length_scale, q);
            fname = save_dir_path_spindoctor + "/" + fname;
            
            % Load existing or compute new decomposition
            if isfile(fname)
                % Load Bloch-Torrey eigendecomposition
                disp("load " + fname);
                load(fname);
            else
                % Compute Bloch-Torrey eigendecomposition
                bt_eig_tmp = compute_blochtorrey_eig(q, lap_eig, setup.gradient.directions);

                % Extract results
                Vsort = bt_eig_tmp.Vsort;
                Dsort = bt_eig_tmp.Dsort;
                invVsort = bt_eig_tmp.invVsort;
                invVsortC1 = bt_eig_tmp.invVsortC1;
                totaltime = bt_eig_tmp.totaltime;

                % Save BT eigendecomposition
                disp("save " + fname + " -v7.3 -struct bt_eig_tmp");
                save(fname, "-v7.3", "-struct", "bt_eig_tmp");
            end

            % Store results
            bt_eig.Vsort(iamp, iseq, :) = Vsort;
            bt_eig.Dsort(iamp, iseq, :) = Dsort;
            bt_eig.invVsort(iamp, iseq, :) = invVsort;
            bt_eig.invVsortC1(iamp, iseq, :) = invVsortC1;
            bt_eig.totaltime(iamp, iseq) = totaltime;
        end
    end

    % Clear temporary variables
    clear bt_eig_tmp
    clear Vsort Dsort invVsort invVsortC1 totaltime
    clear fname
end


% % Copy identical matrices to avoid excess memory usage (if values_type = "q")
% for iseq = 2:nsequence
%     for iamp = 1:namplitude
%         for idir = 1:ndirection
%             bt_eig.Vsort{iamp, iseq, idir} = bt_eig.Vsort{iamp, 1, idir};
%             bt_eig.Dsort{iamp, iseq, idir} = bt_eig.Dsort{iamp, 1, idir};
%             bt_eig.invVsort{iamp, iseq, idir} = bt_eig.invVsort{iamp, 1, idir};
%             bt_eig.invVsortC1{iamp, iseq, idir} = bt_eig.invVsortC1{iamp, 1, idir};
%         end
%     end
% end
