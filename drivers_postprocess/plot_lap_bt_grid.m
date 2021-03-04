disp(["Plotting the of grid of LAP and EIG significant eigenfunctions"]);
disp("mf_signal_reduce uses a subset of Lap eigs and a subset of BTPDE eigs only");
disp("squeeze(real(mf.signal_allcmpts))/volume_allcmpts uses the full set of Lap eigs and BTPDE eigs");

color_vec{1} = "k";
color_vec{2} = "b";
color_vec{3} = "r";
color_vec{4} = "c";
color_vec{5} = "g";
color_vec{6} = "y";
ncolor = length(color_vec);

marker_vec{1} = "d";
marker_vec{2} = "s";
marker_vec{3} = "+";
marker_vec{4} = "d";
marker_vec{5} = "s";
marker_vec{6} = "+";
nmarker = length(marker_vec);

units.t = "ms";
units.q = "\mu m^{-1}";
units.eigval = "ms^{-1}";


maxkappa = max(setup.pde.permeability);

neig = length(lap_eig.length_scales);
neig_lap = find(lap_eig.length_scales >= 1.0, 1, "last");
eigindex_use_btpde = 1:min(neig_lap, neig);
eigindex_use_lap = 1:min(neig_lap, neig);
eiglength_vec = lap_eig.length_scales(eigindex_use_lap);
eiglength_vec_full = lap_eig.length_scales;

experi_vec = 1:nsequence;

b_vec = 1:namplitude;

sequences = setup.gradient.sequences;

%% Plot BT eigenvalues
clear mf_signal_reduce mf_signal_reduce_fromgrid
clear bt_eig_values bt_eig_significant
clear legend_vec

dir_vec = 1;%:ndirection;

figure;
hold on;

% a = get(gcf, "position");
% set(gcf, "position", 0.7*a);

legend_vec = {};
iplot = 0;
for iseq = 1:nsequence
    for iamp = flip(b_vec)
        for idir = dir_vec
            Dsort = bt_eig.Dsort{iamp, iseq, idir};
            Vsort = bt_eig.Vsort{iamp, iseq, idir};
            invVsortC1 = bt_eig.invVsortC1{iamp, iseq, idir};
            bt_eig_values = Dsort(eigindex_use_btpde);
            bt_eig_significant = 1:50;
            % bt_eig_significant = abs(Vsort(1, eigindex_use_btpde)) >= 1e-2;
            % bt_eig_significant = abs(invVsortC1(eigindex_use_btpde)) >= 1e-2;
            % bt_eig_significant = 1:length(bt_eig_values);
            values = bt_eig_values(bt_eig_significant);
            % values = bt_eig_values;
            % values = 1000 * values;
            values = values / mean_diffusivity;
            h = plot(real(values), imag(values), color_vec{iamp} + marker_vec{iamp});
            set(h, "markersize", 20 - iamp * 3);
            set(h, "linewidth", 1 + iamp / 2);
            if length(dir_vec) > 1
                leg_dir_str = sprintf("dir %d, ", idir);
            else
                leg_dir_str = "";
            end
            gamma = 2.67513 * 1e-04;
            % legend_vec{iplot + 1} = sprintf("%sq=%g%s", leg_dir_str, setup.gradient.qvalues(iamp, iseq), units.q);
            % legend_vec{iplot + 1} = sprintf("%sq=%g", leg_dir_str, setup.gradient.qvalues(iamp, iseq));
            legend_vec{iplot + 1} = sprintf("%sg=%.3f", leg_dir_str, setup.gradient.qvalues(iamp, iseq) / gamma);
            iplot = iplot + 1;
        end
    end
end
% xlim([0 1]);
% ylim([-1, 1]);

xlim([0 0.3]);
grid on;
if length(dir_vec) > 1
    dir_str = sprintf("dir [" + join(repmat("%d", 1, length(dir_vec))) + "]", dir_vec);
else
    dir_str = sprintf("dir=[%.2f; %.2f; %.2f]", setup.gradient.directions.points(:, dir_vec));
end

% xlabel(sprintf("Real part (%s)", units.eigval));
% ylabel(sprintf("Imaginary part (%s)", units.eigval));
xlabel("Real part");
ylabel("Imaginary part");
legend(legend_vec);
set(gca, "fontsize", 14);
% title_str = sprintf("Bloch-Torrey eigenvalues, \\kappa=%g, %s", maxkappa, dir_str);
title_str = sprintf("\\kappa=%g, %s", maxkappa, dir_str);
title(title_str, "fontsize", 18);

% idir = dir_vec(1);
% exportgraphics(gca, sprintf("output/k%g_dir%d.png", maxkappa, idir));
% exportgraphics(gca, sprintf("output/k%g_dir%d_wrap.png", maxkappa, idir));


%% Plot real part of BT eigenvalues against invVsortC1
clear mf_signal_reduce mf_signal_reduce_fromgrid
clear bt_eig_values

figure;
hold on
legend_vec = {};
iplot = 0;
for iseq = experi_vec
    for iamp = b_vec
        for idir = dir_vec
            Dsort = bt_eig.Dsort{iamp, iseq, idir};
            invVsortC1 = bt_eig.invVsortC1{iamp, iseq, idir};
            bt_eig_values = real(Dsort(eigindex_use_btpde));
            values = 1000 * bt_eig_values;
            h = plot(log10(abs(invVsortC1(eigindex_use_btpde))), values, ...
                color_vec{iamp} + marker_vec{idir});
            set(h, "markersize", 20 - 3 * iamp);
            set(h, "linewidth", 1 + iamp / 2);
            legend_vec{iplot+1} = sprintf("dir %d, t=%g%s, q=%g%s", idir, ...
                0 / 1000, units.t, setup.gradient.qvalues(iamp, iseq), units.q);
            iplot = iplot + 1;
        end
    end
end
% xlim([-5, 1]);
% ylim([0, 1]);
xlabel("BT eigenfunction coefficient magnitude (log 10)");
ylabel(sprintf("BT eigenvalue, real part (%s)", units.eigval));
legend(legend_vec, "location", "northwest", "fontsize", 14); set(gca, "fontsize", 14);
grid on;


%% Plot Laplace eig significance
clear mf_signal_reduce mf_signal_reduce_fromgrid;
clear bt_eig_values;
for iamp = b_vec
    figure;
    hold on;
    legend_vec = {};
    iplot = 0;
    for idir = dir_vec
        for iseq = 1:nsequence%experi_vec
            seq = setup.gradient.sequences{iseq};
            sdelta = seq.delta;
            bdelta = seq.Delta;
            Dsort = bt_eig.Dsort{iamp, iseq, idir};
            Vsort = bt_eig.Vsort{iamp, iseq, idir};
            invVsortC1 = bt_eig.invVsortC1{iamp, iseq, idir};
            Vbtdecay = diag(exp(-sdelta*Dsort(eigindex_use_btpde))) * invVsortC1(eigindex_use_btpde);
            Vlapdecaysqrt = sqrt(exp(-(bdelta - sdelta) * lap_eig.values(eigindex_use_lap)));
            clap = Vsort(eigindex_use_lap, eigindex_use_btpde) * Vbtdecay;
            vvec = Vlapdecaysqrt .* Vsort(eigindex_use_lap, eigindex_use_btpde) * Vbtdecay;

            aaa = vvec' * vvec;
            aaa_noinf = vvec(2:end)' * vvec(2:end);
            mattmp = repmat(transpose(Vbtdecay), [length(eigindex_use_lap), 1]);
            mattmp2 = diag(Vlapdecaysqrt) * Vsort(eigindex_use_lap, eigindex_use_btpde) .* mattmp;
            mat = abs((mattmp2(:, :)));
            matmax = sqrt(aaa);
            mattmp2_keeplarge = mattmp2;
            mattmp2_keeplarge(abs(mattmp2) <= 5e-3 * matmax) = 0;
            vvec2 = sum(mattmp2_keeplarge, 2);
            aaa2 = sum(abs(vvec2).^2);

            lap_eig_significant = abs(clap) >= 0.01;
            clap_copy = nan(size(clap));
            clap_copy(lap_eig_significant) = clap(lap_eig_significant);
            h = plot(lap_eig.length_scales(eigindex_use_lap(2:end)), ...
                abs(clap_copy(2:end)), ...
                color_vec{mod(iseq - 1, ncolor) + 1} + marker_vec{mod(iseq - 1, nmarker) + 1});

            set(h, "markersize", 20 - 3 * iamp);
            set(h, "linewidth", 1 + iamp / 2);
            legend_vec{iplot + 1} = sprintf("dir %d, t=%g%s, q=%g%s", idir, sdelta / 1000, units.t, setup.gradient.qvalues(iamp, iseq), units.q);

            iplot = iplot + 1;
        end
    end
    set(gca, "ylim", [0, 1]);
    set(gca, "xlim", [0, 400]);
    xlabel("LAP eigenvalue length scale");
    ylabel("LAP eigenfunction coefficient magnitude (log 10)");
    legend(legend_vec, "location", "northwest", "fontsize", 14); set(gca, "fontsize", 14);
    grid on;
end



%% Plot 4
clear mf_signal_reduce mf_signal_reduce_fromgrid
clear bt_eig_values
lap_eig_significant_TEover2_size = zeros(max(experi_vec), max(b_vec), max(dir_vec));
lap_eig_significant_TEover2_minindex = zeros(max(experi_vec), max(b_vec), max(dir_vec));
lap_eig_significant_TEover2_maxindex = zeros(max(experi_vec), max(b_vec), max(dir_vec));
lap_eig_significant_TEover2_minlength = zeros(max(experi_vec), max(b_vec), max(dir_vec));

mf_signal_reduce = zeros(nsequence, namplitude, ndirection);
mf_signal_reduce_fromgrid = zeros(nsequence, namplitude, ndirection);
for iamp = b_vec
    figure;
    hold on;
    clear legend_vec
    iplot = 0;
    for idir = dir_vec
        for iseq = experi_vec
            seq = setup.gradient.sequences{iseq};
            sdelta = seq.delta;
            bdelta = seq.Delta;
            Dsort = bt_eig.Dsort{iamp, iseq, idir};
            Vsort = bt_eig.Vsort{iamp, iseq, idir};
            invVsortC1 = bt_eig.invVsortC1{iamp, iseq, idir};

            Vbtdecay = diag(exp(-sdelta*Dsort(eigindex_use_btpde))) * invVsortC1(eigindex_use_btpde);
            Vlapdecaysqrt = sqrt(exp(-(bdelta - sdelta) * lap_eig.values(eigindex_use_lap)));
            clap = Vsort(eigindex_use_lap, eigindex_use_btpde) * Vbtdecay;
            vvec = Vlapdecaysqrt .* Vsort(eigindex_use_lap, eigindex_use_btpde) * Vbtdecay;

            aaa = vvec' * vvec;
            aaa_noinf = vvec(2:end)' * vvec(2:end);

            mf_signal_reduce(iamp, iseq, idir) = aaa;

            vp = nan(size(vvec));
            ii = find(abs(vvec).^2>=0.0001);
            if ~isempty(ii) && aaa >= 0.01
                lap_eig_significant_TEover2{iseq}{iamp}{idir} = ii;
                vp(ii) = vvec(ii);

                lap_eig_significant_TEover2_size(iamp, iseq, idir) = length(ii);

                lap_eig_significant_TEover2_minindex(iamp, iseq, idir) = min(ii);
                lap_eig_significant_TEover2_maxindex(iamp, iseq, idir) = max(ii);
                lap_eig_significant_TEover2_minlength(iamp, iseq, idir) = eiglength_vec(max(ii));

                aaa2 = vvec(ii)' * vvec(ii);
                mf_signal_reduce_fromgrid(iamp, iseq, idir) = aaa2;

                h = plot(eiglength_vec(2:end), abs(vp(2:end)), ...
                    color_vec{mod(iseq - 1, ncolor) + 1} + marker_vec{mod(iseq - 1, nmarker)+1});
                set(h, "markersize", 8 + 3 * iseq);
                set(h, "linewidth", iseq);
                legend_vec{iplot + 1} = ("dir=" + num2str(idir) + ", t = " + num2str(sdelta / 1000) + "ms, q=" ...
                    + num2str(setup.gradient.qvalues(iamp, iseq)) + ", \lambda_0 contrib=" ...
                    + num2str(abs(vvec(1)), 4));
                iplot = iplot + 1;
            end
        end
    end
    xlim([0, 400]);
    ylim([0, 1.1]);
    xlabel("LAP eigenvalue length scale");
    ylabel("Contribution to signal");
    legend(legend_vec, "location", "northwest", "fontsize", 14);
    set(gca, "fontsize", 14);
    grid on;
end

disp("MF subset of BT and LAP eigs");
disp(squeeze(mf_signal_reduce(experi_vec, b_vec, dir_vec)));
disp("MF subset of BT and LAP eigs, remove the small ones");
disp(squeeze(mf_signal_reduce_fromgrid(experi_vec, b_vec, dir_vec)));

disp("Number of significant lap eigens");
disp(squeeze(lap_eig_significant_TEover2_size(experi_vec, b_vec, dir_vec)));
disp("Min index of significant lap eigens");
disp(squeeze(lap_eig_significant_TEover2_minindex(experi_vec, b_vec, dir_vec)));
disp("Max index of significant lap eigens");
disp(squeeze(lap_eig_significant_TEover2_maxindex(experi_vec, b_vec, dir_vec)));
disp("Min length of significant lap eigens");
disp(squeeze(lap_eig_significant_TEover2_minlength(experi_vec, b_vec, dir_vec)));



%% Plot Laplace agains BT indices
hh = gcf;
ifigure0 = hh.Number;
ifigure = 0;
clear mf_signal_reduce mf_signal_reduce_fromgrid
clear bt_eig_values
for idir = dir_vec
    ifigure = ifigure + 1;
    figure(ifigure0 + 3 * (ifigure-1) + 1);
    % figure(ifigure0 + 3 * (ifigure-1) + 2);
    hold on;
    iplot = 0;
    for iamp = b_vec
        for iseq = experi_vec
            seq = setup.gradient.sequences{iseq};
            sdelta = seq.delta;
            bdelta = seq.Delta;

            Dsort = bt_eig.Dsort{iamp, iseq, idir};
            Vsort = bt_eig.Vsort{iamp, iseq, idir};
            invVsortC1 = bt_eig.invVsortC1{iamp, iseq, idir};

            Vbtdecay = diag(exp(-sdelta * Dsort(eigindex_use_btpde))) * invVsortC1(eigindex_use_btpde);
            clap = Vsort(eigindex_use_lap, eigindex_use_btpde) * Vbtdecay;
            Vlapdecaysqrt = sqrt(exp(-(bdelta - sdelta) * lap_eig.values(eigindex_use_lap)));
            vvec = Vlapdecaysqrt .* clap;

            aaa = vvec' * vvec;
            aaa_noinf = vvec(2:end)' * vvec(2:end);

            mattmp = repmat(transpose(Vbtdecay), [length(eigindex_use_lap), 1]);
            mattmp2 = diag(Vlapdecaysqrt) * Vsort(eigindex_use_lap, eigindex_use_btpde) .* mattmp;
            mat = abs(mattmp2(:, :));
            matmax = sqrt(aaa);
            mattmp2_keeplarge = mattmp2;

            [ii, jj] = find(abs(mattmp2) > 0.0001);

            lap_eig_significant_maxindex(iamp, iseq, idir) = max(ii);
            lap_eig_significant_minlength(iamp, iseq, idir) = eiglength_vec(max(ii));
            bt_eig_significant_maxindex(iamp, iseq, idir) = max(jj);
            bt_eig_significant_maxvalue(iamp, iseq, idir) = bt_eig.Dsort{iamp, iseq, idir}(max(jj));
            lap_eig_significant_maxvalue(iamp, iseq, idir) = lap_eig.values(max(ii));

            mattmp2_keeplarge(abs(mattmp2) <= 0.0001) = 0;
            vvec2 = sum(mattmp2_keeplarge, 2);
            aaa2 = sum(abs(vvec2).^2);

            mf_signal_reduce(iamp, iseq, idir) = aaa;
            mf_signal_reduce_fromgrid(iamp, iseq, idir) = aaa2;

            bt_eig_values = real(Dsort(eigindex_use_btpde));

            iplot = iplot + 1;

            % figure(ifigure0 + 3 * (ifigure-1) + 1);
            % subplot(length(b_vec), length(experi_vec), iplot);
            % h=pcolor(lap_eig.length_scales(eigindex_use_lap(2:end)), ...
            %     1000 * bt_eig_values, ...
            %     abs(mattmp2_keeplarge(2:end, :).'));
            % set(h, "edgecolor", "none");
            % colormap(flipud(gray));
            % matmax2 = max(max(abs(mat(1:end, :))));
            % pmax = abs(vvec(1));
            % colorbar;
            % caxis([0, pmax / 100]);
            % % caxis([0, 0.01]);
            % % set(gca, "ylim", [0, 3]);
            % xlabel("LAP eig length scale");
            % % set(gca, "xlim", [0, 3]);
            % ylabel("- BT eig value (ms^{-1})");
            % title(sprintf("gdir %d, sdelta=%g, g=.2e",...
            %     idir, sdelta, setup.gradient.qvalues(1, iamp))
            % %axis equal;
            % %grid on;
            % axis square;

            figure(ifigure0 + 3 * (ifigure - 1) + 1);
            subplot(length(b_vec), length(experi_vec), iplot);
            spy(abs(mattmp2_keeplarge(1:end, 1:end).'));
            title("idir=" + num2str(idir) + ", sdelta=" + num2str(sdelta) ...
                + ", g=" + num2str(setup.gradient.qvalues(iamp, iseq)));
            xlabel("LAP eig index");
            ylabel("BT eig index");
        end
    end
end

% disp("BTPDE");
% disp(squeeze(real(btpde.signal_allcmpts(dir_vec, experi_vec, b_vec))) / volume_allcmpts);
% disp("MF full set of BT and LAP eigs");
% disp(squeeze(real(mf.signal_allcmpts(experi_vec, b_vec, dir_vec))) / volume_allcmpts);
disp("MF subset of BT and LAP eigs");
disp(squeeze(mf_signal_reduce(b_vec, experi_vec, dir_vec)));
disp("MF subset of BT and LAP eigs, remove the small ones");
disp(squeeze(mf_signal_reduce_fromgrid(b_vec, experi_vec, dir_vec)));

disp("Max index of significant LAP eigens");
disp(squeeze(lap_eig_significant_maxindex(b_vec, experi_vec, dir_vec)));
disp("Min length of significant LAP eigens");
disp(squeeze(lap_eig_significant_minlength(b_vec, experi_vec, dir_vec)));
disp("Max index of significant BT eigens");
disp(squeeze(bt_eig_significant_maxindex(b_vec, experi_vec, dir_vec)));




%%
clear mf_signal_reduce mf_signal_reduce_fromgrid
clear bt_eig_values
mf_signal_reduce = zeros(nsequence, namplitude, ndirection);
mf_signal_reduce_fromgrid = zeros(nsequence, namplitude, ndirection);
for idir = dir_vec
    figure;
    hold on;
    iplot = 0;
    gdir = setup.gradient.directions.points(:, idir);
    for iamp = b_vec
        W_mat = sum(lap_eig.moments .* shiftdim(gdir, -2), 3);
        L_mat = diag(lap_eig.values);
        for iseq = experi_vec
            qval = setup.gradient.qvalues(iamp, iseq);
            seq = setup.gradient.sequences{iseq};
            sdelta = seq.delta;
            bdelta = seq.Delta;

            Dsort = bt_eig.Dsort{iamp, iseq, idir};
            Vsort = bt_eig.Vsort{iamp, iseq, idir};
            invVsortC1 = bt_eig.invVsortC1{iamp, iseq, idir};

            V_mat = L_mat + sqrt(-1) * W_mat * qval;
            [Vsmall, Dsmall] = eig(V_mat(eigindex_use_lap, eigindex_use_lap));
            Dsmall = diag(Dsmall);
            c1small = zeros(size(Vsmall, 1), 1);
            c1small(1) = 1;
            invVsmallC1 = Vsmall \ c1small;

            Vbtdecay = diag(exp(-sdelta * Dsort)) * invVsortC1;
            clap = Vsort * Vbtdecay;
            Vlapdecaysqrt = sqrt(exp(-(bdelta - sdelta) * lap_eig.values));
            vvec = Vlapdecaysqrt .* clap;
            aaa = vvec' * vvec;

            Vbtdecay = diag(exp(-sdelta * Dsmall)) * invVsmallC1;
            clap = Vsmall * Vbtdecay;
            Vlapdecaysqrt = sqrt(exp(-(bdelta - sdelta) * lap_eig.values(eigindex_use_lap)));

            vvec = Vlapdecaysqrt .* clap;
            aaa = vvec'*vvec;
            mf_signal_reduce(iamp, iseq, idir) = aaa;
            aaa_noinf = vvec(2:end)' * vvec(2:end);

            mattmp = repmat(transpose(Vbtdecay), [length(eigindex_use_lap), 1]);
            mattmp2 = diag(Vlapdecaysqrt) * Vsmall .* mattmp;
            mat = abs(mattmp2(:, :));
            matmax = sqrt(aaa);
            mattmp2_keeplarge = mattmp2;

            [ii, jj] = find(abs(mattmp2) > 0.001);
            if ~isempty(ii)
                lap_eig_significant_maxindex(iamp, iseq, idir) = max(ii);
                lap_eig_significant_minlength(iamp, iseq, idir) = eiglength_vec(max(ii));
                bt_eig_significant_maxindex(iamp, iseq, idir) = max(jj);
                bt_eig_significant_maxvalue(iamp, iseq, idir) = bt_eig.Dsort{iamp, iseq, idir}(max(jj));
                lap_eig_significant_maxvalue(iamp, iseq, idir) = lap_eig.values(max(ii));
            else
                lap_eig_significant_maxindex(iamp, iseq, idir) = nan;
                lap_eig_significant_minlength(iamp, iseq, idir) = nan;
                bt_eig_significant_maxindex(iamp, iseq, idir) = nan;
                bt_eig_significant_maxvalue(iamp, iseq, idir) = nan;
                lap_eig_significant_maxvalue(iamp, iseq, idir) = nan;
            end

            mattmp2_keeplarge(abs(mattmp2) <= 0.001) = 0;
            vvec2 = sum(mattmp2_keeplarge, 2);
            aaa2 = sum(abs(vvec2).^2);

            mf_signal_reduce_fromgrid(iamp, iseq, idir) = aaa2;
            bt_eig_values = real(Dsort(eigindex_use_btpde));

            iplot = iplot + 1;

            subplot(length(b_vec), length(experi_vec), iplot);
            spy(abs(mattmp2_keeplarge(1:end, 1:end).'));
            title("idir=" + num2str(idir) + ", sdelta=" + num2str(sdelta) ...
                + ", g=" + num2str(setup.gradient.qvalues(iamp, iseq)));
            xlabel("LAP eig index");
            ylabel("BT eig index");
            if (~isempty(ii))
                set(gca, "xlim", [1, 2 * lap_eig_significant_maxindex(iamp, iseq, idir)]);
                set(gca, "ylim", [1, 2 * bt_eig_significant_maxindex(iamp, iseq, idir)]);
            end
        end
    end
end

disp(["BTPDE signal"]);
signal_btpde_signal = real((btpde.signal_allcmpts)) / sum(volumes);
disp(squeeze(signal_btpde_signal));
disp(["MF signal"]);
signal_mf_signal = real((mf.signal_allcmpts)) / sum(volumes);
disp(squeeze(signal_mf_signal));

disp(["REL diff between BTPDE and MF (percent)"]);
signal_mf_btpde_relerror = 100 * abs(signal_btpde_signal - signal_mf_signal) ./ abs(signal_btpde_signal);
disp(squeeze(signal_mf_btpde_relerror));
disp(["REL diff between BTPDE and MF: mean rel error (percent)"]);
signal_mf_btpde_meanrelerror = (mean(signal_mf_btpde_relerror, 1));
disp(squeeze(signal_mf_btpde_meanrelerror));

clear bmat;
for iseq = experi_vec
    for iamp = b_vec
        seq = setup.gradient.sequences{iseq};
        qv = setup.gradient.qvalues(iamp, iseq);
        bmat(iamp, iseq) = qv^2 * seq.bvalue_no_q;
    end
end

signal_mf_MFreduce_relerr = 100 * abs(signal_mf_signal - mf_signal_reduce) ./ signal_mf_signal;
signal_mf_meanvalue = squeeze(mean(signal_mf_signal, 1));
signal_mf_MFreduce_meanrelerror = squeeze(mean(signal_mf_MFreduce_relerr, 1));

signal_mfreduce_fromgrid_relerr = 100 * abs(mf_signal_reduce - mf_signal_reduce_fromgrid) ./ mf_signal_reduce;
signal_mfreduce_fromgrid_meanrelerror = squeeze(mean(signal_mfreduce_fromgrid_relerr, 1));

signal_btpde_MFreduce_relerr = 100 * abs(signal_btpde_signal - mf_signal_reduce)./signal_btpde_signal;
signal_btpde_meanvalue = squeeze(mean(signal_btpde_signal, 1));
signal_btpde_MFreduce_meanrelerror = squeeze(mean(signal_btpde_MFreduce_relerr(dir_vec, experi_vec, b_vec, 1)));

figure; hold on;
for iamp = b_vec
    plot(bmat(iamp, :), signal_mf_MFreduce_meanrelerror(iamp, :), color_vec{iamp} + marker_vec{iamp});
end
figure; hold on;
for iamp = b_vec
    plot(signal_mf_meanvalue(iamp, :), signal_mf_MFreduce_meanrelerror(iamp, :), color_vec{iamp} + marker_vec{iamp});
end

disp("MF full set of eigs: meanvalue");
disp(squeeze(signal_mf_meanvalue));
disp("BTPDE: meanvalue");
disp(squeeze(signal_btpde_meanvalue));

disp("MF full - MF subset of eigs: rel error (percent)");
disp(squeeze(signal_mf_MFreduce_meanrelerror));
disp("BTPDE - MF subset of eigs: rel error (percent)");
disp(squeeze(signal_btpde_MFreduce_meanrelerror));

disp("MF subset - MF subset remove small: rel error (percent)");
disp(squeeze(signal_mfreduce_fromgrid_meanrelerror));

load BTPDErefinement0pt1.mat
signal_btpde_h0pt1 = signal_btpde;

disp("BTPDE - BTPDE Ref: rel error (percent)");
disp(squeeze(100 * mean(abs(signal_btpde_h0pt1 - signal_btpde_signal) ./ abs(signal_btpde_h0pt1), 1)));

disp("MF full - BTPDE Ref: rel error (percent)");
disp(squeeze(100 * mean(abs(signal_btpde_h0pt1 - signal_mf_signal) ./ abs(signal_btpde_h0pt1), 1)));

disp("MF reduce - BTPDE Ref: rel error (percent)");
disp(squeeze(100 * mean(abs(signal_btpde_h0pt1 - mf_signal_reduce) ./ abs(signal_btpde_h0pt1), 1)));
