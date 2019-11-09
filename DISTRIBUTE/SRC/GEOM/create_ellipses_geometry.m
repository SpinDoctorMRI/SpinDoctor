function [ncell,facets_cell,facets_labels_cell,nodes_cell,pt_in_cell,center,normal,Rcell] ...
    = create_ellipses_geometry(fname_cell,Rratio_nucleus)

ndim = 3;

fid=fopen(fname_cell);
tline = fgetl(fid);
tline = fgetl(fid);
ncell= sscanf(tline,'%f',1);
tline = fgetl(fid);
for icell = 1:ncell
    tline = fgetl(fid);
    tline = fgetl(fid);
    tline = fgetl(fid);
    atmp= sscanf(tline,'%f');
    center{icell}(1,1:3) = atmp(1:3);
    normal{icell}(1,1:3) = atmp(4:6);
    Rcell(icell,1:3) = atmp(7:9);    
end
fclose(fid);

if (Rratio_nucleus < 0 | Rratio_nucleus > 1)
    disp(['Rratio_nucleus out of range']);
end

if (Rratio_nucleus > 0)
    nregion = 2*ncell;
else
    nregion = ncell;
end

pt_in_cell = zeros(nregion,ndim);

for icell = 1:ncell
    pt_in_cell(icell,1:3) = center{icell}(1,1:3)+[Rcell(icell,1)*(Rratio_nucleus+(1-Rratio_nucleus)/2),0,0];
end

if (Rratio_nucleus > 0)
    for icell = 1:ncell
        pt_in_cell(icell+ncell,1:3) = center{icell}(1,1:3);
    end
end

Rmean = mean(mean(Rcell(:,:)));

for icell = 1:ncell
    Np = round(200*(Rcell(icell,1)/Rmean)^2);
    [nodes,facets,vol]=spheresurface_regularpoints(Rcell(icell,1),Np);    
    Nb = size(nodes,1);
	nodes(1:Nb,1) = nodes(1:Nb,1)+center{icell}(1,1);
    nodes(1:Nb,2) = nodes(1:Nb,2)+center{icell}(1,2);
    nodes(1:Nb,3) = nodes(1:Nb,3)+center{icell}(1,3);
	facets_cell{icell} = facets';
	nodes_cell{icell} = nodes;	
    
    facets_labels_cell{icell} = ones(1,size(facets_cell{icell},2));

    offset = Nb;
    
    if (Rratio_nucleus > 0)
        Np = round(200*(Rcell(icell,1)*Rratio_nucleus/Rmean)^2);
        [nodes,facets,vol]=spheresurface_regularpoints(Rcell(icell,1)*Rratio_nucleus,Np);
        Nb = size(nodes,1);
        nodes(1:Nb,1) = nodes(1:Nb,1)+center{icell}(1,1);
        nodes(1:Nb,2) = nodes(1:Nb,2)+center{icell}(1,2);
        nodes(1:Nb,3) = nodes(1:Nb,3)+center{icell}(1,3);
        facets_cell{icell} = [facets_cell{icell},facets'+offset];
        nodes_cell{icell} = [nodes_cell{icell};nodes];
        facets_labels_cell{icell} = [facets_labels_cell{icell},2*ones(1,size(facets',2))];
        
    end
    
end


