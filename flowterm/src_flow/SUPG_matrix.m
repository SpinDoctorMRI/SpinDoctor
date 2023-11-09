function [C_stable,M_stable,A_stable] = SUPG_matrix(vv_cmpts,...
                                                    mv_cmpts,...
                                                    Jxx_cmpts,...
                                                    Jy_cmpts,...
                                                    Jz_cmpts,...
                                                    hk_cmpts,...
                                                    norm2_v_cmpts,...
                                                    ug,...
                                                    thetadt,...
                                                    diffusivity,...
                                                    elements_cmpts,...
                                                    tautype)
                                                            
%SUPG_MATRIX Summary of this function goes here
%   Detailed explanation goes here
    ncompartment = size(vv_cmpts,2);
    C_stable_cmpts = {cell(1, ncompartment)};
    M_stable_cmpts = {cell(1, ncompartment)};
    Jx_stable_cmpts = repmat({cell(1, ncompartment)}, 1, 3);
    for icmpt = 1:ncompartment
        elements = elements_cmpts{icmpt};
        elements = elements.';
        
        NE=size(elements,1);
        NLB=4; %number of local basic functions, it must be known!
        
        vv = vv_cmpts{icmpt};

        mv = mv_cmpts{icmpt};

        Jx = Jxx_cmpts{icmpt};
        Jy = Jy_cmpts{icmpt};
        Jz = Jz_cmpts{icmpt};
        
        hk = hk_cmpts{icmpt};
        norm2_v = norm2_v_cmpts{icmpt};
        
        tau_k = zeros(NE,1);
        
        dif = diffusivity(1, 1, icmpt);
        for indele = 1:NE
            if tautype==1
                % type 1
                tau_k(indele,1) = hk(indele)*hk(indele)/(4*dif*thetadt+ 2*thetadt*norm2_v(indele)*hk(indele) );
            elseif tautype==2
                % type 2
                xitmp2 = hk(indele)*norm2_v(indele)/(6*dif);
                if xitmp2>1
                    xitmp2=1;
                end
                tau_k(indele,1) = hk(indele)*xitmp2/( 2*thetadt*norm2_v(indele) );
            elseif tautype==3
                % type 3
                tau_k(indele,1) = min([hk(indele)/(2*thetadt*norm2_v(indele)), hk(indele)*hk(indele)/(thetadt*dif)  ]);
            else
                error("only three types of stabilization paramters supported.")
            end
        end

        Z = astam(tau_k,vv);
        massZ = astam(tau_k,mv);
        JxZ = astam(tau_k,Jx);
        JyZ = astam(tau_k,Jy);
        JzZ = astam(tau_k,Jz);

        Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);
        X=permute(Y,[2 1 3]);
        
        C_stable_cmpts{icmpt}=sparse(X(:),Y(:),Z(:)); 
        M_stable_cmpts{icmpt}=sparse(X(:),Y(:),massZ(:)); 
        Jx_stable_cmpts{1}{icmpt}=sparse(X(:),Y(:),JxZ(:)); 
        Jx_stable_cmpts{2}{icmpt}=sparse(X(:),Y(:),JyZ(:)); 
        Jx_stable_cmpts{3}{icmpt}=sparse(X(:),Y(:),JzZ(:)); 
    end

    C_stable = blkdiag(C_stable_cmpts{:});
    M_stable = blkdiag(M_stable_cmpts{:});
    Jx_stable = cellfun(@(J) blkdiag(J{:}), Jx_stable_cmpts, "UniformOutput", false);
    A_stable = ug(1) * Jx_stable{1} + ug(2) * Jx_stable{2} + ug(3) * Jx_stable{3};
end

