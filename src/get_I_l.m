function I_l = get_I_l(D,R0,l,setup)
seq = setup.gradient.sequences{1};
t = seq.diffusion_time;

beta = R0./sqrt(D*t);






I_l = A(l,beta).*exp(-beta.^2/4)./(4*pi*D*t).^(3/2) + B(l,beta).*erf(beta/2)./(4*pi/R0^3);

end

function a = A(l,beta)
    if l==0
        a = 1;
    elseif l==2
        a= -(1+6*beta.^(-2));
    elseif l==4
        a = 1 + 20*beta.^(-2) + 210*beta.^(-4);
    elseif l==6
        a=-(1+42*beta.^(-2) + 1574/2 + beta.^(-4) + 10395*beta.^(-6));
    elseif l==8
        a = 1 + 72*beta.^(-2) + 10395/4 * beta.^(-4) + 45045*beta.^(-6) +675675*beta.^(-8);
    else
        error('A_%d not implemented yet',l);
    end
end
function b = B(l,beta)
    if l==0
        b = 0;
    elseif l==2
        b= 3;
    elseif l==4
        b = 15/2 * (1 - 14*beta.^(-2));
    elseif l==6
        b= 105/8 * (1-36*beta.^(-2) + 396* beta.^(-4));
    elseif l==8
        b = 315/16 * (1 - 66*beta.^(-2) + 1716*beta.^(-4) - 17160*beta.^(-6));
    else
        error('B_%d not implemented yet',l);
    end
end