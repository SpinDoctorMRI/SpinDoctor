function difftime = seqdifftime
 
  
  global BDELTA SDELTA SEQ OGSEPER 
  global PGSE OGSEsin OGSEcos dPGSE
  
  
  omega = 1./OGSEPER*2*pi;
    
  if (SEQ == PGSE)
    difftime = (BDELTA-SDELTA/3);
  elseif (SEQ == OGSEcos)
    difftime = 1/8*OGSEPER;
  elseif (SEQ == OGSEsin)
    difftime = 3/8*OGSEPER;
  elseif (SEQ == dPGSE)
    difftime = 2*(BDELTA-SDELTA/3);
  else
    disp('error in seqdifftime');
    stop
  end