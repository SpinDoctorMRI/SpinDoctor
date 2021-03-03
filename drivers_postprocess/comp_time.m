%% MF signal
% [~, mf_signal_allcmpts, mf_ctime] = compute_mf_signal_sub_bt(experiment, initial_signal_domains, lap_eig, directions, bt_eig);
[~, mf_signal_allcmpts, mf_ctime] = compute_mf_signal(experiment, initial_signal_domains, lap_eig);

%% Display computational time
disp("Laplace eigendecomposition:");
disp(lap_eig.ctime);
% disp("BT eigendecomposition (per direction):");
% disp(bt_eig.ctime / ndirection);
disp("MF computational time (per direction):");
disp(shiftdim(mean(mf_ctime, 4), 1));

%%
disp("BTPDE computational time (per direction):");
disp(mean(btpde.ctime, 3));
