function [VOL,EleV,EleC] = get_volume_surface(Pts,Ele);
  
% define areas of boundary edges (EA) in FE mesh and out pointing 
%normals (ENx,ENy) to boundary at the edges
% SA is the total surface area of boundary
% SAx is its projection onto the x-axis.
  
    
  % define volumes of triangles in FE mesh
  % VOL is the total volume of cmpt
    
  x1 = Pts(1,Ele(1,:)); y1 = Pts(2,Ele(1,:)); z1 = Pts(3,Ele(1,:)); 
  x2 = Pts(1,Ele(2,:)); y2 = Pts(2,Ele(2,:)); z2 = Pts(3,Ele(2,:)); 
  x3 = Pts(1,Ele(3,:)); y3 = Pts(2,Ele(3,:)); z3 = Pts(3,Ele(3,:));
  x4 = Pts(1,Ele(4,:)); y4 = Pts(2,Ele(4,:)); z4 = Pts(3,Ele(4,:));
  
  NEle = size(Ele,2);
  EleV = zeros(1,NEle);
  EleC = zeros(3,NEle);
  for iEle = 1:NEle
    aa = [x1(iEle),y1(iEle),z1(iEle)];
    bb = [x2(iEle),y2(iEle),z2(iEle)];
    cc = [x3(iEle),y3(iEle),z3(iEle)];
    dd = [x4(iEle),y4(iEle),z4(iEle)]; 
    
    EleV(iEle) = 1/6*abs(dot((aa-dd),(cross(bb-dd,cc-dd))));
    EleC(1,iEle) = mean([x1(iEle),x2(iEle),x3(iEle),x4(iEle)]);
    EleC(2,iEle) = mean([y1(iEle),y2(iEle),y3(iEle),y4(iEle)]);
    EleC(3,iEle) = mean([z1(iEle),z2(iEle),z3(iEle),z4(iEle)]);
    
  end
  
  VOL = sum(EleV);



