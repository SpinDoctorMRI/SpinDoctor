function [cell_shape,Rratio_nucleus,dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    ic_nucleus,ic_cytoplasm,ic_exterior,ic_llimit,ic_ulimit,kappa_nc,kappa_ce,include_box,box_gap,...
    create_geom,fname_geom,ncell,Hcyl,Rmean,Rmin,Rmax,Htetgen,para_deform] ...
    = read_simulation_parameters_domain(fname_domain)

ndim = 3;
fid=fopen(fname_domain);

tline = fgetl(fid);
cell_shape = sscanf(tline,'%f',1);

tline = fgetl(fid);
Rratio_nucleus = sscanf(tline,'%f',1);
  
if (Rratio_nucleus < 0 | Rratio_nucleus > 0.99)
    Rratio_nucleus = 0;
end
  
tline = fgetl(fid);
dcoeff_nucleus = sscanf(tline,'%f',1);

tline = fgetl(fid);
dcoeff_cytoplasm = sscanf(tline,'%f',1);

tline = fgetl(fid);
dcoeff_exterior = sscanf(tline,'%f',1);

% initial conditions
tline = fgetl(fid);
ic_nucleus = sscanf(tline,'%f',1);

tline = fgetl(fid);
ic_cytoplasm = sscanf(tline,'%f',1);

tline = fgetl(fid);
ic_exterior = sscanf(tline,'%f',1);

tline = fgetl(fid);
ic_llimit = sscanf(tline,'%f',ndim);

tline = fgetl(fid);
ic_ulimit = sscanf(tline,'%f',ndim);

tline = fgetl(fid);
kappa_nc = sscanf(tline,'%f',1);
tline = fgetl(fid);
kappa_ce = sscanf(tline,'%f',1);
tline = fgetl(fid);
include_box = sscanf(tline,'%f',1);
tline = fgetl(fid);
box_gap = sscanf(tline,'%f',1);

tline = fgetl(fid);
create_geom = sscanf(tline,'%f',1);

tline = fgetl(fid);
[strpos] = regexp(tline,"'");
fname_geom = tline(strpos(1)+1:strpos(2)-1);

tline = fgetl(fid);
ncell = sscanf(tline,'%f',1);
tline = fgetl(fid);
Hcyl= sscanf(tline,'%f',1);
tline = fgetl(fid);
Rmean = sscanf(tline,'%f',1);
tline = fgetl(fid);
Rmin = sscanf(tline,'%f',1);
tline = fgetl(fid);
Rmax = sscanf(tline,'%f',1);
tline = fgetl(fid);
Htetgen = sscanf(tline,'%f',1);

tline = fgetl(fid);
para_deform = sscanf(tline,'%f',2);

fclose(fid);
    
 