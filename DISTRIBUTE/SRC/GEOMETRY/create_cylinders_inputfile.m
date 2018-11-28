% Write cylinder information input file for reading in pde solver
function create_cylinders_inputfile(ncyl,Rmean,Rmax,Rmin,Hmax,fname)


% this is for rotating cylinder slice normals
Rotx = @(theta) [1,0,0; 0,cos(theta),-sin(theta);0,sin(theta),cos(theta)];
Roty = @(theta) [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
Rotz = @(theta) [cos(theta),-sin(theta),0; sin(theta),cos(theta),0;0,0,1];

  
positions = (rand(ncyl*100000,2)-0.5)*Rmean*sqrt(ncyl)*2;

pos_x = positions(:,1);
pos_y = positions(:,2);
npos = size(positions,1);
  
center2d(1,1:2) = [pos_x(1),pos_y(1)];
R2d(1,1) = Rmean;

icyl = 1;
for ipos = 1:npos	
	dist = sqrt((center2d(1:icyl,1)-pos_x(ipos)).^2+(center2d(1:icyl,2)-pos_y(ipos)).^2)-R2d(1:icyl,1);
	Ruse = min(dist)*0.85;
	if (Ruse > Rmin)
		icyl = icyl+1;
		center2d(icyl,1:2) = [pos_x(ipos),pos_y(ipos)];
		R2d(icyl,1) = min(Rmax,Ruse);
	end
end
if (icyl < ncyl)
	disp(['did not find enough centers']);
	stop
end
t = linspace(0,1,10)*2*pi;

figure; hold on;

for icyl = 1:ncyl 
	plot(center2d(icyl,1)+R2d(icyl,1)*cos(t),center2d(icyl,2)+R2d(icyl,1)*sin(t),'-o');
end


rr = mean(R2d);

for icyl = 1:ncyl 

	nslice{icyl} = 50;

	Rcyl{icyl} = zeros(nslice{icyl}+1,1);
	tvec = linspace(0,1,nslice{icyl}+1)';
	
	Rcyl{icyl} = (1-0.05*(1-sin(tvec*2*pi*10))/2)*R2d(icyl,1);

	%Rcyl{icyl} = ones(size(tvec))*R2d(icyl,1);
		
	Hcirc = linspace(0,Hmax,nslice{icyl}+1)';
	
	center{icyl} = zeros(nslice{icyl}+1,3);
	
	for islice = 1:nslice{icyl}+1
		center{icyl}(islice,1:2) = center2d(icyl,:);
	end
	center{icyl}(:,3) = Hcirc;

	normal{icyl} = zeros(nslice{icyl}+1,3);
	normal{icyl}(:,3) = 1;
	for islice = 1:nslice{icyl}+1
		vec = Rotx((islice-1)/nslice{icyl}*2*pi*0.15)*normal{icyl}(islice,:)';
		normal{icyl}(islice,:) = vec';
	end

end


fid = fopen(fname,'w');
disp(['Opening ',fname]);
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