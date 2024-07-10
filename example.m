% Add SpinDoctor to path
addpath(genpath('setups'));
addpath(genpath('src'));

% Get paths of cells to be considered
fid = fopen("cells_human.txt",'r');                  % Run only on the full cells
% fid = fopen("cells_human_separated.txt",'r');      % Run only on individual somas and processes
                                                     

tline = fgetl(fid);
i=1;
clear names;
while ischar(tline)
    names(i) = string(tline);
    tline = fgetl(fid);
    i = i + 1;
end
ncells = length(names);


%% Matrix formalism experiments
% Parameters for experiments
tetgen_options = "-pq1.2a0.5VCn";
ls = '1.0';

% Setup file
setup_file='setup_morez';

for i =1:ncells
    tic
    driver_cell(names(i),setup_file,tetgen_options,ls);
    toc
end

%% Finite element ode15 experiments
% Parameters for experiments
tetgen_options = "-pq1.2a0.1VCn";
full_cell = '1'; % segmenting the full cell mesh not currently supported for human microglia
% Setup file
setup_file='setup_morez_ref_sol';

for i =1:ncells
    tic
    driver_btpde(names(i),setup_file,tetgen_options,full_cell);
    toc
end

%% Output error plots
% Parameters of tetgen and solve_mf to be tested.
lsval = [1.0,3.0];
hval = [0.5];

for i =1:ncells
    tic
    driver_get_morez_errors(names(i),lsval,hval,'0','figures_errors_morez_Q3');   
    close all;
    toc
end

%% Output matrix formalism signals, sorted into Amoeboid vs Ramified
clear type;
for i = 1:ncells
    file = names(i);
    sep_file = split(file,"/");
    type(i) = sep_file(3);
end
tetgen_options = "-pq1.2a0.5VCn";
ls = '1.0';
setup_file = 'setup_morez';
for i = 1:ncells
    output_signals(names(i),setup_file,tetgen_options,ls,sprintf('human_signals_morez_Q3/%s',type(i)),sprintf('final_signals_matlab/%s',type(i)));
end


