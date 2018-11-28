function [y] = mynum2str_fname(x,n) 
  
  nx = length(x);
  y = '';
  
  for ix = 1:nx
    if (abs(x(ix))<1e-14)
      p = '0';
    else
      m = floor(log10(x(ix)));
      if (abs(m) >= 3)
        %p = [num2str(round((x(ix)).*10.^(-m)*10^n)./10^n),'x10^{',int2str(m),'}'];
        p = [num2str(round((x(ix)).*10.^(-m)*10^n)./10^n),'e',int2str(m)];
      else
        p = [num2str(round((x(ix)).*10^n)/10^n)];  
      end
    end
    if (ix > 1)
      y = cat(2,y,',',p);
    else
      y = p;
    end
  end
  