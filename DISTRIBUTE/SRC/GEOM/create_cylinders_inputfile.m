% Write cylinder information input file for reading in pde solver
function create_cylinders_inputfile(ncyl,Rmin,Rmax,dmin,dmax,fname,Hmax)


% % this is for rotating cylinder slice normals
% Rotx = @(theta) [1,0,0; 0,cos(theta),-sin(theta);0,sin(theta),cos(theta)];
% Roty = @(theta) [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
% Rotz = @(theta) [cos(theta),-sin(theta),0; sin(theta),cos(theta),0;0,0,1];

Rmean = mean(Rmin,Rmax);

positions = (rand(ncyl*100000,2)-0.5)*Rmean*max(10,sqrt(ncyl))*2;

pos_x = positions(:,1);
pos_y = positions(:,2);
npos = size(positions,1);
  
center2d(1,1:2) = [pos_x(1),pos_y(1)];
R2d(1,1) = Rmean;

icyl = 1;
for ipos = 1:npos
    if (icyl < ncyl)
        dist = sqrt((center2d(1:icyl,1)-pos_x(ipos)).^2+(center2d(1:icyl,2)-pos_y(ipos)).^2)-R2d(1:icyl,1);
		
		Ruse = min(dist);
		Rmax1 = Ruse - dmin*Rmean;
		Rmin1 = Ruse - dmax*Rmean;
		
		lpt = max(Rmin,Rmin1);
		rpt = min(Rmax,Rmax1); 
			
		if (lpt <= rpt)
			Ruse = (lpt+rpt)/2;
			
        %if (Ruse >= Rmin*0.9 & Ruse <= Rmax*1.1)
            icyl = icyl+1;
            center2d(icyl,1:2) = [pos_x(ipos),pos_y(ipos)];
            R2d(icyl,1) = Ruse;
        end
    end
end


if (icyl < ncyl)
	disp(['did not find enough centers']);
	stop
end



xmin = min(center2d(1:ncyl,1)-R2d(1:ncyl,1));
xmax = max(center2d(1:ncyl,1)+R2d(1:ncyl,1));
ymin = min(center2d(1:ncyl,2)-R2d(1:ncyl,1));
ymax = max(center2d(1:ncyl,2)+R2d(1:ncyl,1));

x0 = (xmin+xmax)/2;
y0 = (ymin+ymax)/2;
center2d(:,1) = center2d(:,1)-x0;
center2d(:,2) = center2d(:,2)-y0;



rr = mean(R2d);

for icyl = 1:ncyl 

	nslice{icyl} = 2*round(Hmax);

	Rcyl{icyl} = zeros(nslice{icyl}+1,1);
	tvec = linspace(0,1,nslice{icyl}+1)';
    
	Rcyl{icyl} = ones(size(tvec))*R2d(icyl,1);
		
	Hcirc = linspace(0,Hmax,nslice{icyl}+1)'-Hmax/2;
	
	center{icyl} = zeros(nslice{icyl}+1,3);
	
	for islice = 1:nslice{icyl}+1
		center{icyl}(islice,1:2) = center2d(icyl,:);
	end
	center{icyl}(:,3) = Hcirc;

    center{icyl}(:,1:2) = center{icyl}(:,1:2);
    
	normal{icyl} = zeros(nslice{icyl}+1,3);
	normal{icyl}(:,3) = 1;


end


fid = fopen(fname,'w');
disp(['Writing geometry to ',fname]);
fprintf(fid, '%s\n', 'Total number of cylinders: ');
fprintf(fid, '%d\n', ncyl);
fprintf(fid, '%s\n', 'Number of slices in each cylinder: ');
fprintf(fid, '%d ', nslice{1:end});
fprintf(fid, '\n');
	
for icyl = 1:ncyl

	ind = [1:nslice{icyl}+1]';
	fprintf(fid, '%s %d:\n', 'Cylinder ', icyl);
	fprintf(fid, '%s\n', 'Index, Center, Normal, Radius');
	if (fid ~= -1) 
		fprintf(fid, '%d %26.16f %26.16f %26.16f %26.16f %26.16f %26.16f %26.16f\n', [ind, center{icyl}, normal{icyl}, Rcyl{icyl}]');
	end
end
fclose(fid);