%% MF signal
% [~, mf_signal_allcmpts, mf_ctime] = compute_mf_signal_sub_bt(experiment, initial_signal_domains, lap_eig, directions, bt_eig);
[mf_signal_allcmpts, mf_itertimes] = compute_mf_signal(setup, initial_signal, lap_eig);

%% Display computational time
disp("Laplace eigendecomposition:");
disp(lap_eig.totaltime);
% disp("BT eigendecomposition (per direction):");
% disp(bt_eig.ctime / ndirection);
disp("MF computational time (per direction):");
disp(shiftdim(mean(mf_itertimes, 3), 1));

%%
disp("BTPDE computational time (per direction):");
disp(mean(btpde.itertimes, 3));
