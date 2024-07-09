function Y = harmonic(u,l,m)
    i = sqrt(-1);
    x = u(:,1); y = u(:,2); z = u(:,3);
    em = zeros(1,l+1); em(abs(m)+1) = 1;
    leg = @(t) squeeze(pagemtimes(em,legendre(l,t,'norm')))';
    if m< 0
        Y =(-1)^m*((x+1i*y)./sqrt(1-z.^2)).^abs(m) .* leg(z)/sqrt(2*pi);

    else
        Y = ((x+1i*y)./sqrt(1-z.^2)).^m .* leg(z)/sqrt(2*pi);

    end

end