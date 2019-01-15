function qmesh=mesh_quality( node_xyz, tetra_node) 
% Copyright (c) 2019, Van-Dang Nguyen

% % The mesh quality is based on the minimum, over all tetrahedrons, of 3 
% % times the radius of the insphere hin divided by the radius of the circumsphere hout.
% % Algorithm to compute hin and hout: https://en.wikipedia.org/wiki/Tetrahedron#Volume

tetra_xyz=cell(4,1);

for idv=1:4
    tetra_xyz{idv}=node_xyz(:,tetra_node(idv,:));
end;

% All vectors formed by 4 vertices
v21 = tetra_xyz{2} - tetra_xyz{1}; v31 = tetra_xyz{3} - tetra_xyz{1};
v41 = tetra_xyz{4} - tetra_xyz{1}; v32 = tetra_xyz{3} - tetra_xyz{2};
v42 = tetra_xyz{4} - tetra_xyz{2}; v43 = tetra_xyz{4} - tetra_xyz{3};


% volume of the tetrahedron
vol = 1/6*abs(dot(v41,cross(v42,v43)));

% Length of all edges
l21 = vecnorm(v21,2); l31 = vecnorm(v31,2);
l41 = vecnorm(v41,2); l32 = vecnorm(v32,2);
l42 = vecnorm(v42,2); l43 = vecnorm(v43,2);    

% Compute InRadius based on the volume and areas of 4 faces
% Each face area is computed by Heron's Formula
p123 = 0.5*(l21 + l31 + l32); p234 = 0.5*(l32 + l42 + l43);
p341 = 0.5*(l43 + l41 + l31); p124 = 0.5*(l21 + l41 + l42);

a123 = sqrt( abs(p123.*(p123-l21).*(p123-l31).*(p123-l32)) );
a234 = sqrt( abs(p234.*(p234-l32).*(p234-l42).*(p234-l43)) );
a341 = sqrt( abs(p341.*(p341-l43).*(p341-l41).*(p341-l31)) );
a124 = sqrt( abs(p124.*(p124-l21).*(p124-l41).*(p124-l42)) );

hin = 3*vol./(a123 + a234 + a341 + a124); 

% Compute Circumradius
r1 =  l21 .* l43 + l31 .* l42 + l41 .* l32;
r2 = -l21 .* l43 + l31 .* l42 + l41 .* l32;
r3 =  l21 .* l43 - l31 .* l42 + l41 .* l32;
r4 =  l21 .* l43 + l31 .* l42 - l41 .* l32;

hout = sqrt(r1.*r2.*r3.*r4)./(24*vol);

% Compute the mesh quality
tetrahedron_quality = 3.0 * hin ./ hout;
value_max  = max  ( tetrahedron_quality );
value_min  = min  ( tetrahedron_quality );
value_mean = mean ( tetrahedron_quality );
value_var  = var  ( tetrahedron_quality );
qmesh.quality=[value_min, value_mean, value_max, value_var];
qmesh.hout = hout;
qmesh.hin = hin;
qmesh.vol = vol;

