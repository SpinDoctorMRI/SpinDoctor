function [shape] = shape_func(point,etype,Noele,NLB)
% Copyright (c) 2013, Talal Rahman, Jan Valdman

% SHAPEDER Returns the gradients of the shape functions with
%          respect to the reference coordinates (xi,eta,...).
%
%  point : point(nod,nop), the coordinates of the
%          points on the reference element.
% dshape : dshape(nod,nos,nop), the gradients of the
%          shape functions (second) at all points (third)
%          with respect to the reference cordinates.
%  etype : 'P0','P1','P2', etc., the element type.
%         
%          Note: 
%          nod - dimension of the element.
%          nop - number of points.
%          nos - number of shape functions.
%

nod = size(point,1);
nop = size(point,2);

switch nod
    
case{2}
    l1 = point(1,:);
    l2 = point(2,:);
    l3 = 1 - l1 - l2;
    
    switch etype
         
     case {'P1'}
         shape_parent_ele = [l3; l1; l2];
         
         shape_parent_ele = reshape(shape_parent_ele,1,3,nop);
         
     case {'P2'}
         %l1=x, l2=y, l3=1-x-y are barycentric coordinates
          %6 basis nodal functions
         %on point (0,0)        phi1=(2*l3-1)*l3         node 1
         %on point (1,0)        phi2=(2*l1-1)*l1         node 2
         %on point (0,1)        phi3=(2*l2-1)*l2         node 3

         %on point (0.5,0.5)    phi4=4*l1*l2             edge 23
         %on point (0,0.5)      phi5=4*l2*l3            edge 13
         %on point (0.5,0)      phi6=4*l1*l3             edge 12
         
         
         shape_parent_ele = [(2*l3-1).*l3;...
                            (2*l1-1).*l1;...
                            (2*l2-1).*l2;...
                            4*l1.*l2;...
                            4*l2.*l3;...
                            4*l1.*l3];
         
         shape_parent_ele = reshape(shape_parent_ele,1,6,nop);
         
     otherwise, error('Only P1 and P2 elements implemented.');
     
     end

case {3}
% 3-D elements.

     l1 = point(1,:);
     l2 = point(2,:);
     l3 = point(3,:);
     l4 = 1 - l1 - l2 -l3;
   
     switch etype
         
     case {'P1'}
         shape_parent_ele = [l4; l1; l2; l3];
         
         
         shape_parent_ele = reshape(shape_parent_ele,1,4,nop);
     
     case {'P2'}
     %Quadratic shape functions
     %l1=x, l2=y, l3=z, l4=1-x-y-z are barycentric coordinates
     %10 basis nodal functions
     %on point (0,0,0)        phi1=(2*l4-1)*l4         node 1
     %on point (1,0,0)        phi2=(2*l1-1)*l1         node 2
     %on point (0,1,0)        phi3=(2*l2-1)*l2         node 3      
     %on point (0,0,1)        phi4=(2*l3-1)*l3         node 4
     
     %on point (0,0.5,0.5)    phi5=4*l2*l3             edge 34
     %on point (0,0,0.5)      phi6=4*l3*l4             edge 14
     %on point (0.5,0,0)      phi7=4*l1*l4             edge 12
     %on point (0.5,0.5,0)    phi8=4*l1*l2             edge 23
     %on point (0.5,0,0.5)    phi9=4*l1*l3             edge 24
     %on point (0,0.5,0)      phi10=4*l2*l4            edge 13
     
     
     shape_parent_ele = [(2*l4-1).*l4; ...
            (2*l1-1).*l1; ...
            (2*l2-1).*l2; ...
            (2*l3-1).*l3; ...
            
%              4*l2.*l3;...%5
%              4*l3.*l4;...%6
%              4*l1.*l4;...%7
%              4*l1.*l2;...%8
%              4*l1.*l3;...%9
%              4*l2.*l4];%10
         
          
             4*l1.*l4;...%7
             4*l1.*l2;...%8
             4*l2.*l4;...%10
             4*l3.*l4;...%6
             4*l1.*l3;...%9
             4*l2.*l3];...%5
                                                    
      shape_parent_ele = reshape(shape_parent_ele,1,10,nop);

     otherwise, error('Only P1 and P2 elements implemented.');
     end

end

shape = repmat(shape_parent_ele,[1,1,1,Noele]);
end
