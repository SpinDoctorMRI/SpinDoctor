% Plot Laplace and Bloch-Torrey eigenfunctions.
%
% This plotting script may only be run after running the script
% "driver_mf" or "driver_mf_save_load".


%% Plot Laplace eigenfunctions
disp("Plotting the LAP eigenfunctions");
eigindex_use_lap = 2:5;
plot_diffdir = 0;

maxkappa = max(setup.pde.permeability);
neig = length(lap_eig.values);

% Split Laplace eigenfunctions into compartments
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);
lap_eig_funcs_sep = mat2cell(lap_eig.funcs, npoint_cmpts);

for ieig = eigindex_use_lap
    % Plot
    lscale_str = sprintf("L=%.1f", lap_eig.length_scales(ieig));
    title_str = sprintf("\\kappa=%g, eig %d/%d (%s)", maxkappa, ieig, neig, lscale_str);
    
    plot_field(femesh, lap_eig_funcs_sep, setup.pde.compartments, title_str, ieig);
    % plot_field_everywhere(femesh, lap_eig_funcs_sep, title_str, ieig);
    % set(gca, "fontsize", 14);
    
    % Plot diffusion direction
    if plot_diffdir
        color = "r";
        p = [24 29 0];
        diffdir = squeeze(lap_eig.moments(1, ieig, :));
        diffdir = 10 * diffdir / norm(diffdir, 2);
        % q = quiver3(p(1), p(2), p(3), diffdir(1), diffdir(2), diffdir(3), 0);
        q = quiver3(0, 0, 0.5, diffdir(1), diffdir(2), diffdir(3), 0);
        q.LineWidth = 1.5;
        q.Color = color;
        q.MaxHeadSize = 1;
        % text(p(1), p(2)-3, p(3), "d", "fontsize", 30, "color", color, ...
        %     "VerticalAlignment", "baseline");
    end
    % view(2);
    % caxis([-0.03, 0.05]);
    % title(title_str, "fontsize", 18);
    % exportgraphics(gca, sprintf("output/lap_k%g_%d.png", maxkappa, ieig));
end

% Clear temporary variables
clear lap_eig_funcs_sep


%% Plot BT eigenfunctions
disp("Plotting the BTPDE eigenfunctions");
neig = length(lap_eig.values);
eigindex_use_lap = 1:neig;
eigindex_use_btpde = 1:5;

% Number of points in each compartment
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

for iseq = 1%:nsequence
    seq = setup.gradient.sequences{iseq};
    for iamp = 1:namplitude
        qval = setup.gradient.qvalues(iamp, iseq);
        for idir = 1%[1 17 33 49]
            bt_eig_funcs_sep = cell(1, ncompartment);
            % Clear results from previous iterations
            clear bt_eig_funcs;
            
            % Convert Laplace eigenfunctions to BT eigenfunctions
            bt_eig_funcs(:, eigindex_use_btpde) = lap_eig.funcs(:, eigindex_use_lap)...
                * bt_eig.invVsort{iamp, iseq, idir}(eigindex_use_btpde, eigindex_use_lap).';
            
            % Split BT eigenfunctions into compartments
            bt_eig_funcs_sep = mat2cell(bt_eig_funcs, npoint_cmpts);
            
            bvalue_str = sprintf("q-value=%g", qval);
            dir_str = sprintf("dir=[%.2f; %.2f; %.2f]", setup.gradient.directions.points(:, idir));
            
            % Plot selected eigenfunctions
            for ieig = eigindex_use_btpde
                title_str = sprintf("BT eigenfunctions, %s, %s, %s, eig %d", ...
                    bvalue_str, dir_str, ieig);
                
                % Plot
                % plot_field(femesh, bt_eig_funcs_sep, cmpts_in, cmpts_out, cmpts_ecs, title_str, ieig);
                plot_field_everywhere(femesh, bt_eig_funcs_sep, title_str, ieig);
            end
        end
    end
end

%% Plot BT eigenfunctions on supports
disp("Plotting the BTPDE eigenfunctions and their supports");
neig = length(lap_eig.values);

experi_use = 1;%:nsequence;
b_use = namplitude;%1:namplitude
% dir_use = [1 17 33 49 64];
% dir_use = [33 49];
dir_use = 1;

eigindex_use_lap = 1:neig;
eigindex_use_btpde = [1:10];

support_tol = 0.1;
plot_supports = false;
plot_gdir = true;

maxkappa = max(setup.pde.permeability);

% Number of points in each compartment
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

nodes = [femesh.points{:}];

% Add number of nodes in all previous compartments to each node number in
% each compartment element collection
tmp = num2cell(cumsum([0 npoint_cmpts(1:end-1)]));
elements = cellfun(@plus, femesh.elements, tmp, "UniformOutput", false);

% Concatenate compartment element arrays
elements = [elements{:}];

for iseq = experi_use
    seq = setup.gradient.sequences{iseq};
    for iamp = b_use
        qval = setup.gradient.qvalues(iamp, iseq);
        for idir = dir_use
            clear bt_eig_func bt_eig_func_sum
            
            [~, inds_sort] = sort(abs(bt_eig.Vsort{iamp, iseq, idir}(1, :)), "descend");
            % [~, inds_sort] = sort(abs(bt_eig.invVsortC1{iamp, iseq, idir}.'), "descend");
            % eigindex_use_btpde = inds_sort(eigindex_use_btpde);
            
            bt_eig_funcs(:, eigindex_use_btpde) = lap_eig.funcs(:, eigindex_use_lap)...
                * bt_eig.invVsort{iamp, iseq, idir}(eigindex_use_btpde, eigindex_use_lap).';
            
            % % Take absolute value
            % bt_eig_funcs = abs(bt_eig_funcs);
            
            % Switch sign of eigenfunctions if more negative than positive
            bt_eig_funcs = sign(sum(real(bt_eig_funcs))) .* bt_eig_funcs;
            
            % Sum eigenfunctions (to plot in same figure, assuming disjoint
            % supports)
            bt_eig_funcs_sum = zeros(size(bt_eig_funcs, 1), 1);
            
            support_centers = zeros(3, neig);
            eigindex_str_cell = cell(1, neig);
            for ieig = eigindex_use_btpde
                % Max value (infty norm)
                maxvalue = max(abs(bt_eig_funcs(:, ieig)));
                
                % Identify support
                inodes_support{ieig} = abs(bt_eig_funcs(:, ieig)) > support_tol * maxvalue;
                
                % Approximate center of support
                % support_centers(:, ieig) = mean(nodes(:, inodes_support{ieig}), 2);
                support_centers(:, ieig) = nodes * abs(bt_eig_funcs(:, ieig)) ...
                    / norm(bt_eig_funcs(:, ieig), 1);
                support{ieig} = nodes(:, inodes_support{ieig});
                eigindex_str_cell{ieig} = num2str(ieig);
                
                % % Normalize
                % bt_eig_funcs(:, ieig) = bt_eig_funcs(:, ieig) / maxvalue;
                
                % % Set value to eigenvalue number
                % bt_eig_funcs(:, ieig) = ieig; % (max(eigindex_use_btpde) - ieig) + 1;
                % bt_eig_funcs(:, ieig) = 1;
                
                % % Set values below threshold to zero
                % bt_eig_funcs(~inodes_support{ieig}, ieig) = 0;
                
                % Sum selected eigenfunctions
                bt_eig_funcs_sum = bt_eig_funcs_sum + bt_eig_funcs(:, ieig);
            end
            
            % Split BT eigenfunctions into compartments
            bt_eig_funcs_sum_sep = mat2cell(bt_eig_funcs_sum, npoint_cmpts);
            
            % Title
            % gdir_str = sprintf("dir=[%.2f; %.2f; %.2f]", directions.points(:, idir));
            gdir_str = sprintf("dir %d/%d", idir, ndirection);
            eig_str = sprintf("eigenindex [" ...
                + join(repmat("%d", 1, length(eigindex_use_btpde))) + "]", ...
                eigindex_use_btpde);
            bvalue_str = sprintf("q=%g", qval);
            % title_str = sprintf("BT eigenfunctions, \\kappa=%g, %s,\n%s, %s", maxkappa, bvalue_str, gdir_str, eig_str);
            title_str = sprintf("\\kappa=%g, %s, %s", maxkappa, bvalue_str, gdir_str);
            
            % Plot eigenfunctions
            plot_field_everywhere(femesh, bt_eig_funcs_sum_sep, title_str);
            % plot_field(femesh, bt_eig_funcs_sum_sep,  cmpts_in, cmpts_out, cmpts_ecs, title_str);
            set(gca, "fontsize", 14);
            
            % Plot gradient direction
            if plot_gdir
                color = "r";
                p = [24 29 0];
                dir = 10 * setup.gradient.directions.points(:, idir);
                q = quiver3(p(1), p(2), p(3), dir(1), dir(2), dir(3), 0);
                q.LineWidth = 1.5;
                q.Color = color;
                q.MaxHeadSize = 1;
                text(p(1) - 1, p(2) - 5, p(3), "g", "fontsize", 30, "color", color, ...
                    "VerticalAlignment", "baseline");
            end
            
            % Place text slightly above domain
            hoverheight = 1;
            support_centers(3, eigindex_use_btpde)...
                = max(support_centers(3, eigindex_use_btpde)) + hoverheight;
            
            % Create gradient colormap
            colors = linspace(0, 1, length(eigindex_use_btpde));
            color1 = [1; 0; 0];
            color2 = [1; 0.5; 1];
            colors = (1 - colors) .* color1 + colors .* color2;
            colors = colors(:, randperm(length(eigindex_use_btpde)));
            icol = 1;
            
            % Plot supports (only for cylinders)
            if setup.geometry.cell_shape == "cylinder"
                for ieig = eigindex_use_btpde
                    % % Make rectangle, convex hull or exact outline
                    % points = box_points_top(nodes(:, inodes_support{ieig}));
                    points = convhull_points_top(support{ieig});
                    % points = contour_points_top(nodes, elements, inodes_support{ieig});
                    
                    % Style
                    color = "r";
                    % color = colors(:, icol);
                    linewidth = 1.5;
                    % linewidth = 3 - 2.5 * find(ieig == eigindex_use_btpde, 1) / length(eigindex_use_btpde);
                    
                    % Plot closure of support
                    if plot_supports
                        h = line(points(1, :), points(2, :), points(3, :));
                        h.Color = color;
                        h.LineWidth = linewidth;
                        h.LineStyle = ":";
                    end
                    
                    % Plot eigenvalue index above support center
                    t = text(support_centers(1, ieig), ...
                        support_centers(2, ieig), ...
                        support_centers(3, ieig), ...
                        eigindex_str_cell(ieig));% num2str(icol));
                    t.HorizontalAlignment = "center";
                    t.VerticalAlignment = "baseline";
                    t.FontSize = 20;
                    t.Color = color;
                    
                    icol = icol + 1;
                end
            end
            
            title(title_str, "fontsize", 18);
            view(2);
            
            % exportgraphics(gca, sprintf("output/k%g_q%g_dir%d.png", maxkappa, qval, idir));
            % End of figure
        end
    end
end

%% Functions for marking supports
function points = box_points_top(points)
% Identify corners of rectangle around support points in the x-y plane
pmin = min(points, [], 2);
pmax = max(points, [], 2);
points = [...
    pmin(1) pmin(1) pmax(1) pmax(1) pmin(1);
    pmin(2) pmax(2) pmax(2) pmin(2) pmin(2);
    pmax(3) pmax(3) pmax(3) pmax(3) pmax(3)];
end

function points = convhull_points_top(points)
% Keep support points defining a convex hull around all the support points, and
% arange them sequentially
points(3, :) = max(points(3, :));
k = convhull(points(1:2, :)');
points = points(:, k);
end

function points = contour_points_top(nodes, elements, support_inds)
% Create sequence of points defining the contour of the support in the plane on
% the top of the cylinders

% Remove nodes that are not in the plane
z = max(nodes(3, support_inds));
support_inds = support_inds & nodes(3, :)' > z - 1e-10;

% Find indices of support nodes
support_inds = find(support_inds);

% Join ghost nodes
for iself = 1:length(support_inds) - 1
    inode = support_inds(iself);
    
    % Ghost node if difference between nodes smaller than tolerance
    isghost = vecnorm(nodes(:, inode) - nodes(:, support_inds(iself+1:end))) < 1e-10;
    ighost = iself + find(isghost, 1); % at most one ghost node
    
    % Check if ghost node has been found
    if all(size(ighost))
        ighost = support_inds(ighost);
        
        % Replace all instances of node with ghost node in elements
        elements(elements == inode) = ighost;
        
        % Mark node to be removed in support_inds (with 0)
        support_inds(iself) = 0;
    end
end

% Remove nodes marked by 0 from support inds
support_inds = support_inds(support_inds ~= 0);

% Find elements containing support nodes
elements = shiftdim(elements, -1);
matches = support_inds == elements;

% Keep elements with one facet laying in the plane of nodes; that is - elements
% with three out of four nodes lying in the plane (no flat elements)
ie_keep = sum(matches, [1 2]) == 3;
elements = elements(:, :, ie_keep);
matches = matches(:, :, ie_keep);

% Remove the fourth point of each tetrahedron that is not lying in the plane,
% and obtain the facet defined by the three remaining points
facets = elements(any(matches, 1));
facets = reshape(facets, [3 sum(ie_keep)]);

% Extract the three edges from from each facet
edges = [facets([1 2], :) facets([1 3], :) facets([2 3], :)];

% Remove inner edges of the support (edges shared by two facets). Remaining edges
% will thus be outer edges of the (flat) support, as there are no facets outside
% the support. Equality is tested for both orientations of each edge
equal_edges = all(edges == permute(edges, [1 3 2]), 1);
equal_inverseedges = all(edges == flipud(permute(edges, [1 3 2])), 1);
equal_e_or_ie = equal_edges | equal_inverseedges;
double_edges = sum(equal_e_or_ie, 2) == 2;
keep = ~double_edges;
edges = edges(:, keep);
nedge = size(edges, 2);

% Sort edges by vicinity, starting from first node of first edge
edge_old = edges(:, 1);
begin_connex = 1;
end_connex = [];
for iedge = 2:nedge
    % Fine next edge (if any)
    [~, inext] = find(edges(:, iedge+1:end) == edge_old(2));
    inext = iedge + inext;
    
    % Switch found edge and misplaced edge
    edges(:, [iedge inext]) = edges(:, [inext iedge]);
    
    % Check that the new edge is oriented in the right way
    if edges(2, iedge) == edge_old(2)
        % Turn edge
        edges([1 2], iedge) = edges([2 1], iedge);
    elseif edge_old(2) == edges(1, begin_connex(end))
        % The previous edge has closed a connected contour, which means that the
        % current edge belongs to a new separated domain of the support. Mark
        % index, and start new connected domain.
        end_connex(end + 1) = iedge - 1;
        begin_connex(end + 1) = iedge;
    end
    
    % Prepare edge for next iteration
    edge_old = edges(:, iedge);
end
end_connex(end + 1) = nedge;

% Create list of points ordered according to the edges, and close contour by
% adding the first point at the end (for each connex)
inds = [];
for iconnex = 1:length(begin_connex)
    inds = [inds begin_connex(iconnex):end_connex(iconnex) begin_connex(iconnex)];
end
points = nodes(:, edges(1, inds));
end
