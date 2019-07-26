function [Cmat,Pts] = cylinder_connectivity(Na,Nb)

tol = 1e-6;

if (Na > Nb) 
	Ra = 1;
	Rb = 0.9;
else		
	Ra = 0.9;
	Rb = 1;
end


t = linspace(0,1,Na+1);

xa = Ra*cos(t(1:end-1)*2*pi);
ya = Ra*sin(t(1:end-1)*2*pi);

t = linspace(0,1,Nb+1);

xb = Rb*cos(2*pi*t(1:end-1));
yb = Rb*sin(2*pi*t(1:end-1));

nodes(:,1) = [xa';xb'];
nodes(:,2) = [ya';yb'];

pgon = polyshape({xa,xb},{ya,yb});

tri = triangulation(pgon);

Dmat = tri.ConnectivityList';
Cmat = zeros(size(Dmat));
Pts = tri.Points;

npts = size(Pts,1);

for i = 1:npts
	ind = find(abs(Pts(i,1) - nodes(:,1)) <= tol & abs(Pts(i,2) - nodes(:,2)) <= tol );
	ii = find(Dmat == i);
	Cmat(ii) = ind; 
end



