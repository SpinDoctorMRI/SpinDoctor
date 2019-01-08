function [domain_geom,domain_pde,domain_mesh] = read_simulation_domain_parameters(fname_domain)

ndim = 3;
fid=fopen(fname_domain);

tline = fgetl(fid);
domain_geom.Rratio_nucleus = sscanf(tline,'%f',1);
  
if (domain_geom.Rratio_nucleus < 0 | domain_geom.Rratio_nucleus > 0.99)
    domain_geom.Rratio_nucleus = 0;
end

tline = fgetl(fid);
domain_geom.include_box = sscanf(tline,'%f',1);

tline = fgetl(fid);
domain_geom.box_gap = sscanf(tline,'%f',1);
  
tline = fgetl(fid);
domain_pde.dcoeff_nucleus = sscanf(tline,'%f',1);

tline = fgetl(fid);
domain_pde.dcoeff_cytoplasm = sscanf(tline,'%f',1);

tline = fgetl(fid);
domain_pde.dcoeff_exterior = sscanf(tline,'%f',1);

% initial conditions
tline = fgetl(fid);
domain_pde.ic_nucleus = sscanf(tline,'%f',1);

tline = fgetl(fid);
domain_pde.ic_cytoplasm = sscanf(tline,'%f',1);

tline = fgetl(fid);
domain_pde.ic_exterior = sscanf(tline,'%f',1);

tline = fgetl(fid);
domain_pde.kappa_nc = sscanf(tline,'%f',1);
tline = fgetl(fid);
domain_pde.kappa_ce = sscanf(tline,'%f',1);

tline = fgetl(fid);
domain_mesh.Htetgen = sscanf(tline,'%f',1);


tline = fgetl(fid);
[strpos] = regexp(tline,"'");
domain_mesh.tetgen_cmd = tline(strpos(1)+1:strpos(2)-1);
fclose(fid);
    
 