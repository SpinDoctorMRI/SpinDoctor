function plot_hardi(directions, signal, fig_title)
%PLOT_HARDI Plot signal in all gradient directions.
%
%   directions: struct
%   signal: double(ndirection, 1)
%   fig_title


% Check if all HARDI points lie in the plane (x-y plane)
is_flat = ~any(directions.points(3, :));

% Plot HARDI points in the plane or in 3D-space
if is_flat
    plot_hardi_2d(directions, signal, fig_title)
else
    plot_hardi_3d(directions, signal, fig_title)
end
