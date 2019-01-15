function [YY] = spherical_harmonics(x,y,z,r)

% this generates spherical harmonics (basis functions on the sphere)
if length(x) < 4
    sizeYY = 1;
elseif (length(x) >= 4 && length(x) < 9)
    sizeYY = 4;
elseif length(x) >= 9
    sizeYY = 9;
end

YY = zeros(length(x),sizeYY);

% Y00
YY(:,1) = 1/2*sqrt(1/pi);

if (length(x) >= 4)
    % Y1,{-1,0,1}
    YY(:,2) = sqrt(3/(4*pi))*y./r;
    YY(:,3) = sqrt(3/(4*pi))*z./r;
    YY(:,4) = sqrt(3/(4*pi))*x./r;
    if (length(x) >= 9)
        % Y2,{-2,-1,0,1,2}
        YY(:,5) = 1/2*sqrt(15/pi)*x.*y./(r.^2);
        YY(:,6) = 1/2*sqrt(15/pi)*y.*z./(r.^2);
        YY(:,7) = 1/2*sqrt(15/pi)*x.*z./(r.^2);
        YY(:,8) = 1/4*sqrt(5/pi)*(-x.^2-y.^2+2*z.^2)./(r.^2);
        YY(:,9) = 1/4*sqrt(15/pi)*(x.^2-y.^2)./(r.^2);
    end
end


