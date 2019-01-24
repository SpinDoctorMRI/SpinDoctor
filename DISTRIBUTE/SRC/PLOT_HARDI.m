function PLOT_HARDI(points,ADC_allcmpts_alldir)

ngdir = size(points,1);

[sph_pts,C_sph] = spheresurface_regularpoints(1,900);   
YY = spherical_harmonics(points(:,1),points(:,2),points(:,3),ones(size(points,1),1));
sh_coeff = YY\ADC_allcmpts_alldir;
    
YY_sph = spherical_harmonics(sph_pts(:,1),sph_pts(:,2),sph_pts(:,3),ones(size(sph_pts,1),1));
ADC_interp = YY_sph*sh_coeff;
ADC_alldir = ADC_interp.*sph_pts;
        
figure; hold on;
ADC_Ellipsoid = trisurf(C_sph,ADC_alldir(:,1),ADC_alldir(:,2),ADC_alldir(:,3),ADC_interp); view(3); colorbar; axis equal;
ADC_Ellipsoid.EdgeColor = 'none';
xlabel('x'); ylabel('y'); zlabel('z');
title(['ADC in ',num2str(ngdir),' directions']);
plot3(ADC_allcmpts_alldir.*points(:,1),ADC_allcmpts_alldir.*points(:,2),ADC_allcmpts_alldir.*points(:,3),'k.',...
'markersize',30);