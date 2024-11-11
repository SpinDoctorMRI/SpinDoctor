function S3 = get_mitra_S3(femesh)
%GET_MITRA_SM Computes structural matrix for generalised Mitra formula.
%   Moutal, N., Maximov, I.I. and Grebenkov, D.S., 2019. 
%   Probing surface-to-volume ratio of an anisotropic medium by diffusion NMR with general gradient encoding. 
%   IEEE Transactions on Medical Imaging, 38(11), pp.2507-2522.

S = femesh.total_area;
[areas, facet_centers, normals] = get_surfacenormal_mesh(femesh.points{1}, femesh.elements{1}, femesh.facets{1});
nfaces = size(facet_centers,2);
n = reshape(normals,3,1,nfaces);
m = reshape(normals,1,3,nfaces);
k = pagemtimes(n,m);
areas = reshape(areas,1,1,nfaces);
S3 = sum(k.*areas,3)/S;
end

