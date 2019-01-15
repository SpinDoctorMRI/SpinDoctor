function [facets,nodes,ind_face_a,ind_face_b,facets_labels] = create_cylinders_wall(dr,nslice,Rcyl_one,normal_one,center_one,cyl_closed)

	islice = 0;
	Rb = Rcyl_one(islice+1,1);		
	Nb = round(2*pi*Rb/dr);	
	t = linspace(0,1,Nb+1)';	
	xb = Rb*cos(t(1:end-1)*2*pi);
	yb = Rb*sin(t(1:end-1)*2*pi);
	zb = 0*ones(size(xb));	
	nodes = [xb,yb,zb];
	nnodes = Nb; 	
	nor = normal_one(islice+1,:);
	nor = nor/norm(nor);
	[val,i]=min(abs(nor));
	if (i(1) == 1)
		v1 = cross([1,0,0],nor);
	elseif (i(1) == 2)
		v1 = cross([0,1,0],nor);
	elseif (i(1) == 3)
		v1 = cross([0,0,1],nor);
	end
	v1 = v1/norm(v1);
	v2 = cross(v1,nor);
	Mrot = [v1',v2',nor'];
	mat = repmat(center_one(islice+1,:)',[1,Nb])+Mrot*[xb';yb';zb'];	
	nodes(1:Nb,1:3) = [mat(1,:)',mat(2,:)',mat(3,:)'];		
	Ncirc_one(islice+1,1) = Nb;
	Na = Nb;

	%nslice = nslice_vec(icyl);
	
	for islice = 1:nslice	
		Rb = Rcyl_one(islice+1);
		Nb = round(2*pi*Rb/dr);
		Ncirc_one(islice+1,1) = Nb;
		t = linspace(0,1,Nb+1)';
		xb = Rb*cos(t(1:end-1)*2*pi);
		yb = Rb*sin(t(1:end-1)*2*pi);
%		zb = Hcirc(islice+1)*ones(size(xb));	
		zb = 0*ones(size(xb));
		nor = normal_one(islice+1,:);
		nor = nor/norm(nor);
		[val,i]=min(abs(nor));
		if (i(1) == 1)
			v1 = cross([1,0,0],nor);
		elseif (i(1) == 2)
			v1 = cross([0,1,0],nor);
		elseif (i(1) == 3)
			v1 = cross([0,0,1],nor);
		end
		v1 = v1/norm(v1);
		v2 = cross(v1,nor);
		Mrot = [v1',v2',nor'];
	  	mat = repmat(center_one(islice+1,:)',[1,Nb])+Mrot*[xb';yb';zb'];
		nodes(nnodes+1:nnodes+Nb,1:3) = [mat(1,:)',mat(2,:)',mat(3,:)'];
		
		nnodes = nnodes + Nb;	
		[Cmat,Pts] = cylinder_connectivity(Na,Nb);
		if (islice == 1) 
			ConnectMat = Cmat;
			offset = Na;
		else
			ConnectMat = [ConnectMat,Cmat+offset];
			offset = offset + Na;
		end	
		Na = Nb;
	end

	facets = ConnectMat;
   
    facets_labels = 1*ones(1,size(facets,2));
      
    ind_face_a = 1:Ncirc_one(1,1);
    ind_face_b = nnodes-Ncirc_one(end,1)+1:nnodes;
        
    if (cyl_closed == 0)
       
    else
        ind_face_a_center = nnodes+1;
        ind_face_b_center = nnodes+2;
        [face_a_Cmat,face_a_center] = face_connectivity(nodes(ind_face_a,:),ind_face_a,ind_face_a_center);
        [face_b_Cmat,face_b_center] = face_connectivity(nodes(ind_face_b,:),ind_face_b,ind_face_b_center);
        nodes(ind_face_a_center,1:3) = face_a_center;
        nodes(ind_face_b_center,1:3) = face_b_center;
        facets = [ConnectMat,face_a_Cmat,face_b_Cmat];
        facets_labels = [facets_labels,2*ones(1,size([face_a_Cmat,face_b_Cmat],2))];
    end