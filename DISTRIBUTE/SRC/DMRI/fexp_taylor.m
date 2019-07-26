function y = fexp_taylor(const,time)


  y1 = quad(@(sime) (1+const*(sime-time)+(const*(sime-time)).^2/2+(const*(sime-time)).^3/6).*seqprofile(sime),0,time,1e-2);

  y = const*seqintprofile(time).*y1;