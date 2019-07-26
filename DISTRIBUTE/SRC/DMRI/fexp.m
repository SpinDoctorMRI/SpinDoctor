function y = fexp(const,time,TE)

  y1 = quad(@(sime) exp(const*(sime-time)).*seqprofile(sime),0,time,1e-2);

  y = const*seqintprofile(time).*y1;