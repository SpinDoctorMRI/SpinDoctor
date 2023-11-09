function M=damping_matrix_3D(elements,coordinates,volumes, velocity)
%% 
NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension

elementtype = size(elements,2);

%particular part for a given element in a given dimension
if elementtype == 4
    NLB=4; %number of local basic functions, it must be known!
elseif elementtype ==10
    NLB=10;
end

coord=zeros(DIM,NLB,NE);
for d=1:DIM
    for i=1:NLB           
        coord(d,i,:)=coordinates(elements(:,i),d);       
    end   
end

IP = zeros(3,2^3);

[IP(1,:),IP(2,:),IP(3,:),weight]=tetraquad(2,[0,0,0;1,0,0;0,1,0;0,0,1]);

weight = 6* weight;

if elementtype==4
    phi = shape_func(IP,'P1',NE,NLB);
    dphi = phider(coord,IP,'P1');
elseif elementtype==10
    phi = shape_func(IP,'P2',NE,NLB);
    dphi = phider(coord,IP,'P2');
end

dphidx = dphi(1,:,:,:);
dphidy = dphi(2,:,:,:);
dphidz = dphi(3,:,:,:);

Z=0;
for poi = 1:length(weight)
    dphidxtemp = dphidx(:,:,poi,:);
    dphidx_temp(1,:,:) = squeeze(dphidxtemp);
    
    dphidytemp = dphidy(:,:,poi,:);
    dphidy_temp(1,:,:) = squeeze(dphidytemp);
    
    dphidztemp = dphidz(:,:,poi,:);
    dphidz_temp(1,:,:) = squeeze(dphidztemp);
    
    phitemp = phi(:,:,poi,:);
    phi_temp(1,:,:) = squeeze(phitemp);

    volumes_temp = volumes;
    
%     Z_temp1=astam(velocity(1).*volumes_temp*weight(poi),amtam(dphidx_temp,phi_temp));
%     Z_temp2=astam(velocity(2).*volumes_temp*weight(poi),amtam(dphidy_temp,phi_temp));
%     Z_temp3=astam(velocity(3).*volumes_temp*weight(poi),amtam(dphidz_temp,phi_temp));
% 
%     Z = Z + Z_temp1 + Z_temp2 + Z_temp3;
    
    Z_temp=astam(volumes_temp*weight(poi),amtam( astam(velocity(1,:),dphidx_temp)...
        +astam(velocity(2,:),dphidy_temp)+astam(velocity(3,:),dphidz_temp),phi_temp ) );
    Z = Z + Z_temp;
end

Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);
X=permute(Y,[2 1 3]); 
M=sparse(X(:),Y(:),Z(:)); 
end
