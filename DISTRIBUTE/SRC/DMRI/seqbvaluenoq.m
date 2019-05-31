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
    F2 = @(t) seqintprofile(t).^2;
    bvaluenoq = integral(F2,0,SDELTA+BDELTA);   
  end