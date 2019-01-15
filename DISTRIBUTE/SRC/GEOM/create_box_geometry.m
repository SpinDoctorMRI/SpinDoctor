function [nodes,facets,pt_in] = create_box_geometry(x1min,x1max,x2min,x2max,x3min,x3max)

	nodes = [
      0.0 0.0 0.0;
      1.0 0.0 0.0;
      1.0 1.0 0.0;
      0.0 1.0 0.0;
      0.0 0.0 1.0;
      1.0 0.0 1.0;
      1.0 1.0 1.0;
      0.0 1.0 1.0];
 
	facets = [
      1 2 3; 
      3 4 1;
      5 6 7; 
      7 8 5;
      1 2 6; 
      6 5 1;
      2 3 7; 
      7 6 2;
      3 4 8; 
      8 7 3;
      4 1 5; 
      5 8 4]';
  
	nodes(:,1) = nodes(:,1)*(x1max-x1min)+x1min;
	nodes(:,2) = nodes(:,2)*(x2max-x2min)+x2min;
	nodes(:,3) = nodes(:,3)*(x3max-x3min)+x3min; 
	pt_in(1,1) = (x1min+x1max)/2;
	pt_in(1,2) = (x2min+x2max)/2;
	pt_in(1,3) = (x3min+x3max)/2;