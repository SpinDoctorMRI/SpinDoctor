function [I_total,I_cmpts] = get_inertia_tensor(femesh)
    
    ncompartment = femesh.ncompartment;

    I_cmpts = zeros(ncompartment,3,3);
    for icmpt = 1:ncompartment
        p = femesh.points{icmpt};
        e = femesh.elements{icmpt};
        [~, volumes, centers] = get_volume_mesh(p, e);
        r2 = vecnorm(centers,2,1).^2;
        
        Ixx = (r2 -centers(1,:).^2)*volumes';
        Iyy = (r2-centers(2,:).^2)*volumes';
        Izz = (r2-centers(3,:).^2)*volumes';
        Ixy = centers(1,:).*centers(2,:)*volumes';
        Ixz = centers(1,:).*centers(3,:)*volumes';
        Iyz = centers(2,:).*centers(3,:)*volumes';
        I_cmpts(icmpt,:,:) = [Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];

    end
    
    p = [femesh.points{:}];
    e = [femesh.elements{:}];
    [~, volumes, centers] = get_volume_mesh(p, e);
    r2 = vecnorm(centers,2,1).^2;
    Ixx = (r2 -centers(1,:).^2)*volumes';
    Iyy = (r2-centers(2,:).^2)*volumes';
    Izz = (r2-centers(3,:).^2)*volumes';
    Ixy = (centers(1,:).*centers(2,:))*volumes';
    Ixz = (centers(1,:).*centers(3,:))*volumes';
    Iyz = (centers(2,:).*centers(3,:))*volumes';
    I_total = [Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];

