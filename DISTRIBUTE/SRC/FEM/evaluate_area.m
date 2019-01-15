function elements2area = evaluate_area(elements,coordinates)
% Copyright (c) 2018, Jan Valdman

v1=coordinates(elements(:,2),:)-coordinates(elements(:,1),:);
v2=coordinates(elements(:,3),:)-coordinates(elements(:,1),:);
   
matrix_3D=zeros(2,2,size(elements,1));
matrix_3D(1,1,:)=sum(v1.*v1,2);
matrix_3D(1,2,:)=sum(v1.*v2,2);
matrix_3D(2,1,:)=sum(v1.*v2,2);
matrix_3D(2,2,:)=sum(v2.*v2,2);
   
elements2area=(sqrt(amdet(matrix_3D))/2)';
end 
    
    
