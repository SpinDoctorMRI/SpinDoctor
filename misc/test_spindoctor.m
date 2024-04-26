%TEST_spindoctor
% Run some tests to guarantee SpinDoctor's functionality.

% import SpinDoctor
clear
cd("../");
addpath(genpath("src"));

test_path = "setups";
test_list = [
    "setup_cylinders_mini_settings"
    "setup_spheres_mini_settings"
    "setup_cylinders_full_settings"
    "setup_spheres_full_settings"
    "setup_neuron_mini_settings"
    "setup_neuron_full_settings"
    "setup_1axon_analytical"
    "setup_1sphere_analytical"
    "setup_15spheres"
    "setup_2axons_deform"
    "setup_5axons_myelin_relax"
    "setup_4axons_flat"
    "setup_30axons_flat"
    "setup_30axons"
    "setup_200axons"
    "setup_neuron"
    
];
% Choose whether to save magnetization field
magnetization_flag = true;
for itest = 1:length(test_list)
    try
        run(fullfile(test_path, test_list(itest) + ".m"));

        % Prepare simulation
        [setup, femesh, surfaces]  = prepare_simulation(setup);

        %% Solve BTPDE
        if isfield(setup, "btpde")
            disp("Computing or loading the BTPDE signals");
            savepath = create_savepath(setup, "btpde");
            btpde = solve_btpde(femesh, setup, savepath, magnetization_flag);
            btpde2 = load_btpde(setup, savepath, magnetization_flag);
            
            % test
            btpde = rmfield(btpde, 'totaltime');
            btpde2 = rmfield(btpde2, 'totaltime');
            assert(all(DataHash(btpde) ==  DataHash(btpde2)));
        end

        %% Solve BTPDE using midpoint method
        if isfield(setup, "btpde_midpoint")
            disp("Computing or loading the BTPDE midpoint signals");
            savepath = create_savepath(setup, "btpde_midpoint", 'saved_simul');
            btpde_mp = solve_btpde_midpoint( ...
                femesh, setup, savepath, magnetization_flag ...
            );
            btpde_mp2 = load_btpde_midpoint(setup, savepath, magnetization_flag);

            % test
            btpde_mp = rmfield(btpde_mp, 'totaltime');
            btpde_mp2 = rmfield(btpde_mp2, 'totaltime');
            assert(all(DataHash(btpde_mp) ==  DataHash(btpde_mp2)));
        end

        %% Solve HADC model
        if isfield(setup, "hadc")
            disp("Computing or loading the homogenized apparent diffusion coefficient");
            savepath = create_savepath(setup, "hadc");
            hadc = solve_hadc(femesh, setup, savepath);
            hadc2 = load_hadc(femesh, setup, savepath);

            % test
            hadc = rmfield(hadc, 'totaltime');
            hadc2 = rmfield(hadc2, 'totaltime');
            assert(all(DataHash(hadc) ==  DataHash(hadc2)));
        end

        %% Laplace eigendecomposition
        if isfield(setup, "mf")
            % Laplace eigendecomposition
            eigenpath = create_savepath(setup, "lap_eig");
            lap_eig = compute_laplace_eig(femesh, setup.pde, setup.mf, eigenpath);
            lap_eig2 = load_laplace_eig(eigenpath, setup.mf, setup.pde.mean_diffusivity);
            assert(all(DataHash(lap_eig) ==  DataHash(lap_eig2)));

            % Compute MF magnetization and signal
            savepath = create_savepath(setup, "mf");
            mf = solve_mf(femesh, setup, lap_eig, savepath, magnetization_flag);
            mf2 = load_mf(setup, savepath, magnetization_flag);
            % test
            mf = rmfield(mf, 'totaltime');
            mf2 = rmfield(mf2, 'totaltime');
            assert(all(DataHash(mf) ==  DataHash(mf2)));

            mf_hadc = solve_mf_hadc(femesh, setup, lap_eig);
        end

        assert(max(abs(btpde.signal - btpde_mp.signal)/setup.pde.initial_signal, [], 'all') < 0.1);
        assert(max(abs(btpde.signal - mf.signal)/setup.pde.initial_signal, [], 'all') < 0.1);
        assert(max(abs(hadc.adc_allcmpts - mf_hadc.adc_allcmpts), [], 'all') < 0.1);

        % Solve Karger model
        if isfield(setup, "karger")
            % Solve analytical analytical model
            karger = solve_karger(femesh, setup);
        end
    catch e
        disp(e)
        fprintf('Error in %s', test_list(itest));
    end
end
