% Write ellipses information input file for reading in pde solver
function create_ellipses_inputfile(nell,Rmean,Rmax,Rmin,fname)

% this is for rotating normals
Rotx = @(theta) [1,0,0; 0,cos(theta),-sin(theta);0,sin(theta),cos(theta)];
Roty = @(theta) [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
Rotz = @(theta) [cos(theta),-sin(theta),0; sin(theta),cos(theta),0;0,0,1];

center3d(1,1:3) = [0,0,0]; %[pos_x(1),pos_y(1),pos_z(1)];
R3d(1,1) = Rmean;

iell = 1;

if (nell >= 2)
    
    positions = (rand(nell*50000,3)-0.5)*Rmean*max(10,(nell)^(1/3));
    
    pos_x = positions(:,1);
    pos_y = positions(:,2);
    pos_z = positions(:,3);
    
    npos = size(positions,1);
    
    for ipos = 1:npos
        if (iell < nell)
            dist = sqrt((center3d(1:iell,1)-pos_x(ipos)).^2+(center3d(1:iell,2)-pos_y(ipos)).^2+(center3d(1:iell,3)-pos_z(ipos)).^2)-R3d(1:iell,1);
            Ruse = min(dist)*0.98;
            if (Ruse >= Rmin & Ruse <= Rmax)
                iell = iell+1;
                center3d(iell,1:3) = [pos_x(ipos),pos_y(ipos),pos_z(ipos)];
                R3d(iell,1) = Ruse;
            end
        end
    end
    if (iell < nell)
        disp(['did not find enough centers']);
        stop
    end
    
    xmin = min(center3d(1:nell,1)-R3d(1:nell,1));
    xmax = max(center3d(1:nell,1)+R3d(1:nell,1));
    ymin = min(center3d(1:nell,2)-R3d(1:nell,1));
    ymax = max(center3d(1:nell,2)+R3d(1:nell,1));
    zmin = min(center3d(1:nell,3)-R3d(1:nell,1));
    zmax = max(center3d(1:nell,3)+R3d(1:nell,1));
    
    x0 = (xmin+xmax)/2;
    y0 = (ymin+ymax)/2;
    z0 = (zmin+zmax)/2;
    center3d(:,1) = center3d(:,1)-x0;
    center3d(:,2) = center3d(:,2)-y0;
    center3d(:,3) = center3d(:,3)-z0;
    
    figure; hold on;
    
    for iell = 1:nell
        [X,Y,Z]=ellipsoid(center3d(iell,1),center3d(iell,2),center3d(iell,3),R3d(iell,1),R3d(iell,1),R3d(iell,1),20);
        surf(X,Y,Z);
        view(3);
        axis equal;
    end
    

end
for iell = 1:nell
    Rell{iell} = zeros(1,3);
    Rell{iell}(1,1:3) = R3d(iell,1);
    center{iell} = zeros(1,3);
    center{iell}(1,1:3) = center3d(iell,:);
    normal{iell} = zeros(1,3);
    normal{iell}(1,3) = 1;
    vec = Rotx(0.0)*normal{iell}(1,:)';
    normal{iell}(1,:) = vec';
end
fid = fopen(fname,'w');
disp(['Opening ',fname]);
fprintf(fid, '%s\n', 'Total number of ellipses: ');
fprintf(fid, '%d\n', nell);
fprintf(fid, '\n');

for iell = 1:nell
	fprintf(fid, '%s %d:\n', 'ellipse ', iell);
	fprintf(fid, '%s\n', 'Center, Normal, Radius');
	if (fid ~= -1) 
		fprintf(fid, '%26.16f %26.16f %26.16f %26.16f %26.16f %26.16f %26.16f %26.16f %26.16f\n', ...
            [center{iell}, normal{iell}, Rell{iell}]);
	end
end

fclose(fid);