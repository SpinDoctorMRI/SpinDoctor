function [nodes_new] = deform_domain(nodes,para_deform)

a_bend = para_deform(1);
a_twist = para_deform(2);

nodes_new = nodes';

thvec = nodes_new(:,3)*a_twist;
nodes_new = [(cos(thvec).*nodes_new(:,1)-sin(thvec).*nodes_new(:,2)),...
    (sin(thvec).*nodes_new(:,1)+cos(thvec).*nodes_new(:,2)),nodes_new(:,3)];
 
u1 = [0,0,1]';   
for ii = 1:size(nodes,2)
    u2 = [2*a_bend*nodes(3,ii),0,1]';
    u2 = u2/norm(u2);
    [Mrot] = rotate_nodes(u1,u2);
    nodes_new(ii,1:3) = nodes_new(ii,1:3)*Mrot;
end

nodes_new = [nodes_new(:,1)+a_bend*nodes_new(:,3).^2,nodes_new(:,2),nodes_new(:,3)];

nodes_new = nodes_new';

