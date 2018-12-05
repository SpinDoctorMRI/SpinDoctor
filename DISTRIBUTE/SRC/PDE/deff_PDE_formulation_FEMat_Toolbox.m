function [model_FEM_matrices] ...
    = deff_PDE_formulation_FEMat(Ncmpt,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts)

global UG


for icmpt = 1:Ncmpt
    
    model_orig{icmpt} = createpde();
    
    geometryFromMesh(model_orig{icmpt},Pts_cmpt_reorder{icmpt},Ele_cmpt_reorder{icmpt});
    
    specifyCoefficients(model_orig{icmpt},'m',0,'d',1,'c',DIFF_cmpts(icmpt),'a', 0,'f',0);
    
    nface = model_orig{icmpt}.Geometry.NumFaces; %number of face
    %     normal dot (c grad u) + q*u = g
    for iface = 1:nface
        applyBoundaryCondition(model_orig{icmpt},'neumann','face',iface,'g',@BTDeffNeumann_notime,'q',0,'Vectorized','on'); % 3-D geometry
    end
    
    model_FEM_matrices{icmpt} = assembleFEMatrices(model_orig{icmpt});
    % The 6 FE matrices can be obtained in the following way.
end