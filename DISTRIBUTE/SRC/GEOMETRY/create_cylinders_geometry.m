function [ncell,facets_cell,facets_labels_cell,nodes_cell,pt_in_cell,center,normal,Rcell,nslice_vec,...
     nodes_ind_bottomring,nodes_ind_topring] ...
    = create_cylinders_geometry(fname_cell,Rratio_nucleus)

ndim = 3;

fid=fopen(fname_cell);
  
tline = fgetl(fid);
tline = fgetl(fid);
ncell= sscanf(tline,'%f',1);
tline = fgetl(fid);
tline = fgetl(fid);
nslice_vec = sscanf(tline,'%f',ncell);

for icell = 1:ncell
	tline = fgetl(fid);
	tline = fgetl(fid);
	for islice = 1:nslice_vec(icell)+1
		tline = fgetl(fid);		
		atmp= sscanf(tline,'%f');
		ind = atmp(1);
		center{icell}(ind,1:3) = atmp(2:4);
		normal{icell}(ind,1:3) = atmp(5:7);
		Rcell{icell}(ind,1) = atmp(8);
	end
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
    pt_in_cell(icell,:) = center{icell}(2,:);
end

if (Rratio_nucleus > 0)
    for icell = 1:ncell
        pt_in_cell(icell+ncell,:) = center{icell}(2,:)+[Rcell{icell}(2)*(Rratio_nucleus+(1-Rratio_nucleus)/2),0,0];
    end  
end

for icell = 1:ncell
	Rvec(icell) = mean(Rcell{icell});
end

Rmean = mean(Rvec(1:ncell));

dr = 2*pi*Rmean/30;

for icell = 1:ncell

	clear nodes;
    
	nslice = nslice_vec(icell);
    Rcell_one = Rcell{icell};
    normal_one = normal{icell};
    center_one = center{icell};
   
    if (Rratio_nucleus == 0)        
        cell_closed = 1;
        [facets,nodes,ind_face_a,ind_face_b,facets_labels] = create_cylinders_wall(dr,nslice,Rcell_one,normal_one,center_one,cell_closed);
        ind_face_a_out = ind_face_a;
        ind_face_b_out = ind_face_b;
        
    else
        cell_closed = 1;
        [facets_in,nodes_in,ind_face_a_in,ind_face_b_in,facets_labels_in] ...
            = create_cylinders_wall(dr,nslice,Rratio_nucleus*Rcell_one,normal_one,center_one,cell_closed);
        Na_in = length(ind_face_a_in);
        Nb_in = length(ind_face_b_in);      
        offset = size(nodes_in,1);
        
        face_nlabel = max(facets_labels_in);

       	%unique(facets_labels_in)

        cell_closed = 0;
        [facets_out,nodes_out,ind_face_a_out,ind_face_b_out,facets_labels_out] ...
            = create_cylinders_wall(dr,nslice,Rcell_one,normal_one,center_one,cell_closed);
        Na_out = length(ind_face_a_out);
        Nb_out = length(ind_face_b_out);
        ind_face_a_out = ind_face_a_out+offset;
        ind_face_b_out = ind_face_b_out+offset;
        facets_out = facets_out+offset;
        
        facets_labels_out = facets_labels_out+face_nlabel;      
        
        %unique(facets_labels_out)
        face_nlabel = max(facets_labels_out);

        [Cmat_top,Pts] = cylinder_connectivity(Na_in,Na_out);
        
        
        
        Cmat_top_new = Cmat_top;
        for ii = 1:Na_in
          jj = find(Cmat_top == ii); 
          Cmat_top_new(jj) = ind_face_a_in(ii);
        end
        
        for ii = Na_in+1:Na_in+Na_out
            jj = find(Cmat_top == ii);
            Cmat_top_new(jj) = ind_face_a_out(ii-Na_in);
        end
        
    
        
        [Cmat_bottom,Pts] = cylinder_connectivity(Nb_in,Nb_out);
        
        Cmat_bottom_new = Cmat_bottom;
        for ii = 1:Nb_in
            jj = find(Cmat_bottom == ii);
            Cmat_bottom_new(jj) = ind_face_b_in(ii);
        end
        
        for ii = Nb_in+1:Nb_in+Nb_out
            jj = find(Cmat_bottom == ii);
            Cmat_bottom_new(jj) = ind_face_b_out(ii-Nb_in);
        end
        
        nodes = [nodes_in; nodes_out];
        facets = [facets_in,facets_out,Cmat_top_new,Cmat_bottom_new];
        
        facets_labels = [facets_labels_in,facets_labels_out,(face_nlabel+1)*ones(1,size([Cmat_top_new,Cmat_bottom_new],2))];
         
        %unique(facets_labels)
    end
	
    
	facets_cell{icell} = facets;
	nodes_cell{icell} = nodes;
    facets_labels_cell{icell} = facets_labels;
    nodes_ind_bottomring{icell} = ind_face_a_out;
    nodes_ind_topring{icell} = ind_face_b_out;
	
end
