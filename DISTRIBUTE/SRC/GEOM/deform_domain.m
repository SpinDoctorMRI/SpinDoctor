function [nodes_new] = deform_domain(nodes,para_deform)

a_bend = para_deform(1);
a_twist = para_deform(2);

% a_bend = 0
% a_twist = 0

nodes_new = nodes';

nodes_new = [nodes_new(:,1)+a_bend*nodes_new(:,3).^2,nodes_new(:,2),nodes_new(:,3)];

rmax = max(sqrt(nodes_new(:,1).^2+nodes_new(:,2).^2));

nodes_new(:,1) = nodes_new(:,1)+rmax;

thvec = nodes_new(:,3)*a_twist;

nodes_new = [(cos(thvec).*nodes_new(:,1)-sin(thvec).*nodes_new(:,2)),...
    (sin(thvec).*nodes_new(:,1)+cos(thvec).*nodes_new(:,2)),nodes_new(:,3)];

nodes_new(:,1) = nodes_new(:,1)-rmax;

% nodes_new(:,3) = nodes_new(:,3)*5;
% 
% nodes_new(:,1) = nodes_new(:,1)*4;
% nodes_new(:,2) = nodes_new(:,2)*4;

nodes_new = nodes_new';

