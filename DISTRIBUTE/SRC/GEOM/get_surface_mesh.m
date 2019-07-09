function [SA,SAu] = get_surface_mesh(Pts,Fac,UG);

% define areas of boundary edges (EA) in FE mesh and out pointing
%normals (ENx,ENy) to boundary at the edges
% SA is the total surface area of boundary
% SAx is its projection onto the x-axis.
% define volumes of triangles in FE mesh
% VOL is the total volume of cmpt

x1 = Pts(1,Fac(1,:)); y1 = Pts(2,Fac(1,:)); z1 = Pts(3,Fac(1,:));
x2 = Pts(1,Fac(2,:)); y2 = Pts(2,Fac(2,:)); z2 = Pts(3,Fac(2,:));
x3 = Pts(1,Fac(3,:)); y3 = Pts(2,Fac(3,:)); z3 = Pts(3,Fac(3,:));


NFac = size(Fac,2);

FacA = zeros(1,NFac);
FacC = zeros(3,NFac);
FacN = zeros(3,NFac);

if (nargin >= 3)
    UG_vec = UG;
else
    UG_vec = eye(3);
end

tmp = zeros(size(UG_vec,1),1);

for iFac = 1:NFac
    v1 = [x1(iFac)-x2(iFac),y1(iFac)-y2(iFac),z1(iFac)-z2(iFac)];
    v2 = [x2(iFac)-x3(iFac),y2(iFac)-y3(iFac),z2(iFac)-z3(iFac)];
    vnor = -cross(v1,v2);
    FacA(iFac) = 1/2*norm(vnor);
    vnor = vnor/norm(vnor);
    FacN(1:3,iFac) = vnor;
    FacC(1,iFac) = mean([x1(iFac),x2(iFac),x3(iFac)]);
    FacC(2,iFac) = mean([y1(iFac),y2(iFac),y3(iFac)]);
    FacC(3,iFac) = mean([z1(iFac),z2(iFac),z3(iFac)]);
    
    vnor_vec = repmat(vnor,[size(UG_vec,1),1]);
    %UG_vec
    
    %w = dot(UG_vec,vnor_vec);
    
    %tmp = tmp+dot(UG_vec,vnor_vec)*w*FacA(iFac);
    
    tmp = tmp+(dot(UG_vec,vnor_vec,2)).^2*FacA(iFac);
end
%UG_vec
SAu = tmp';

%vnor = FacN{icmpt}(1:3,iFac);
%SAu_G{icmpt} = SAu_G{icmpt}+dot(UG,vnor)^2*FacA{icmpt}(iFac);

SA = sum(FacA(:));
%figure; h=quiver3(FacC{1}(1,:),FacC{1}(2,:),FacC{1}(3,:),FacN{1}(1,:),FacN{1}(2,:),FacN{1}(3,:));



