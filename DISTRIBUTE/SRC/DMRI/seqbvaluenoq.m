function bvaluenoq = seqbvaluenoq
 
	SEQ_DEFINITIONS
  
  global BDELTA SDELTA SEQ OGSEPER 
  global PGSE OGSEsin OGSEcos dPGSE
  
  
  omega = 1./OGSEPER*2*pi;
    
  if (SEQ == PGSE)
    bvaluenoq = SDELTA^2*(BDELTA-SDELTA/3);
  elseif (SEQ == OGSEcos)
    bvaluenoq = -(cos(omega * SDELTA) ^ 2 * BDELTA * omega - omega * BDELTA - cos(omega * SDELTA) * sin(omega * SDELTA) - omega * SDELTA + 0.2e1 * sin(omega * SDELTA)) / omega ^ 3;
  elseif (SEQ == OGSEsin)
    bvaluenoq = (cos(omega * SDELTA) ^ 2 * BDELTA * omega - 0.2e1 * cos(omega * SDELTA) * BDELTA * omega + 0.2e1 * cos(omega * SDELTA) * omega * SDELTA + omega * BDELTA - cos(omega * SDELTA) * sin(omega ...
                                                      * SDELTA) + omega * SDELTA - 0.2e1 * sin(omega * SDELTA)) / omega ^ 3;
  elseif (SEQ == dPGSE)
    bvaluenoq = 2*SDELTA^2*(BDELTA-SDELTA/3);
  else
%     global sym_s sym_ft sym_SDELTA sym_BDELTA
%     syms sym_s sym_SDELTA sym_BDELTA t
%     ft=subs(sym_ft, {sym_SDELTA, sym_BDELTA}, {SDELTA, BDELTA});
%     F2 = int(ft, sym_s, 0, t)^2;     
%     bvaluenoq = double(int(F2,t,0,2*(SDELTA+BDELTA)));
    F2 = @(t) seqintprofile(t).^2;
    bvaluenoq = integral(F2,0,SDELTA+BDELTA);   
  end