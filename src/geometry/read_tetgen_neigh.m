function neighbours =read_tetgen_neigh(file)
%%READ_TETGEN_NEIGH reads the neighbourhood adjacency information from
%%tetgen path
fid = fopen(sprintf('%s.neigh',file),"r");
fgetl(fid);
neighbours =fscanf(fid,'%d',[5,Inf]);