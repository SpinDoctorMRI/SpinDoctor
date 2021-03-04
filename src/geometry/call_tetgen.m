function call_tetgen(filename, refinement)
%CALL_TETGEN Summary of this function goes here
%   Detailed explanation goes here


% Tetgen command for corresponding operating system
if ispc
    tetgen_cmd = "src\tetgen\win64\tetgen";
elseif ismac
    tetgen_cmd = "src/tetgen/mac64/tetgen";
elseif isunix
    tetgen_cmd = "src/tetgen/lin64/tetgen";
else
    warning("Using Linux Tetgen command.")
    tetgen_cmd = "src/tetgen/tetGen/lin64/tetgen";
end


% Options for Tetgen command
if  nargin == nargin(@call_tetgen) && refinement > 0
    % Pass refinement to the 'a' flag of Tetgen. This gives a maximum
    % tetrahedron volume (not length, as in earlier versions)
    % tetgen_options = "-pq0.5AVa" + num2str(setup.pde.refinement);
    tetgen_options = "-pqAVa" + num2str(refinement);
else
    tetgen_options = "-pqAV";
end

% Call Tetgen
cmd = sprintf("%s %s %s", tetgen_cmd, tetgen_options, filename);
disp(cmd)
system(cmd);
