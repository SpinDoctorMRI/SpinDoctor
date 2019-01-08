function [geom,fname_geom] = read_cells_parameters(fname_geometry)

fid=fopen(fname_geometry);

tline = fgetl(fid);
geom.cell_shape = sscanf(tline,'%f',1);

tline = fgetl(fid);
[strpos] = regexp(tline,"'");
fname_geom = tline(strpos(1)+1:strpos(2)-1);

tline = fgetl(fid);
geom.ncell = sscanf(tline,'%f',1);

tline = fgetl(fid);
geom.Rmin = sscanf(tline,'%f',1);
tline = fgetl(fid);
geom.Rmax = sscanf(tline,'%f',1);

tline = fgetl(fid);
geom.dmin = sscanf(tline,'%f',1);
tline = fgetl(fid);
geom.dmax = sscanf(tline,'%f',1);

tline = fgetl(fid);
geom.para_deform = sscanf(tline,'%f',2);

if (geom.cell_shape == 2)
	tline = fgetl(fid);
	geom.Hcyl= sscanf(tline,'%f',1);
end

fclose(fid);
    
 