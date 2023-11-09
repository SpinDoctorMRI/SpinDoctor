function [vv,mv,Jx,Jy,Jz]=stabilized_matrix_3D(elements,points,volumes, velocity)

% Jx,Jy,Jz : x cdot phi cdot (velocity cdot nabla phi)
% mv : phi cdot (velocity cdot nabla phi)
% vv : (velocity cdot nabla phi) cdot (velocity cdot nabla phi)

% sigma is a 3*3 matrix
%% 
NE=size(elements,1); %number of elements
DIM=size(points,2); %problem dimension
% elementtype = size(elements,2); % p1 or p2
%particular part for a given element in a given dimension
NLB=4; %number of local basic functions, it must be known!
coord=zeros(DIM,NLB,NE);
for d=1:DIM
    for i=1:NLB           
        coord(d,i,:)=points(elements(:,i),d);       
    end   
end
IP = zeros(3,2^3);%zeros(3,4^3);
[IP(1,:),IP(2,:),IP(3,:),weight]=tetraquad(2,[0,0,0;1,0,0;0,1,0;0,0,1]);
weight = 6* weight;
phi = shape_func(IP,'P1',NE,NLB);
dphi = phider(coord,IP,'P1');

dphidx = dphi(1,:,:,:);
dphidy = dphi(2,:,:,:);
dphidz = dphi(3,:,:,:);

[x_global,y_global,z_global] = coord_func(phi,coord);

mv = 0;
vv = 0;
ctildeZ = 0;

Jx = 0;
Jy = 0;
Jz = 0;

for poi = 1:length(weight)
    dphidxtemp = dphidx(:,:,poi,:);
    dphidx_temp(1,:,:) = squeeze(dphidxtemp);
    
    dphidytemp = dphidy(:,:,poi,:);
    dphidy_temp(1,:,:) = squeeze(dphidytemp);
    
    dphidztemp = dphidz(:,:,poi,:);
    dphidz_temp(1,:,:) = squeeze(dphidztemp);
    
    phitemp = phi(:,:,poi,:);
    phi_temp(1,:,:) = squeeze(phitemp);
    
    
    velocity_times_nabla_phi = astam(velocity(1,:),dphidx_temp) ...
                            + astam(velocity(2,:),dphidy_temp) ...
                            + astam(velocity(3,:),dphidz_temp);
    
    vv_temp = amtam(velocity_times_nabla_phi,...
        velocity_times_nabla_phi);
    vv = vv + weight(poi)*vv_temp;

%     mv_temp = amtam(velocity_times_nabla_phi,phi_temp);   
    mv_temp = amtam(phi_temp,velocity_times_nabla_phi);
    mv = mv + weight(poi)*mv_temp;
    
    x_temp = x_global(:,poi);
    y_temp = y_global(:,poi);
    z_temp = z_global(:,poi);
    
    Jx_temp=astam(x_temp,amtam(phi_temp,...
        velocity_times_nabla_phi ));
    Jx = Jx + weight(poi)*Jx_temp;
    
    Jy_temp=astam(y_temp,amtam(phi_temp,...
        velocity_times_nabla_phi));
    Jy = Jy + weight(poi)*Jy_temp;
    
    Jz_temp=astam(z_temp,amtam(phi_temp,...
        velocity_times_nabla_phi));
    Jz = Jz + weight(poi)*Jz_temp;
end

vv = astam(volumes,vv);
mv = astam(volumes,mv);
Jx = astam(volumes,Jx);
Jy = astam(volumes,Jy);
Jz = astam(volumes,Jz);

end