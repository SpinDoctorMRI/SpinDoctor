function difftime = seqdifftime
 
  
  global BDELTA SDELTA SEQ OGSEPER 
  global PGSE OGSEsin OGSEcos
  
  
  omega = 1./OGSEPER*2*pi;
    
  if (SEQ == PGSE)
    difftime = (BDELTA-SDELTA/3);
  elseif (SEQ == OGSEcos)
    difftime = 1/8*OGSEPER;
  elseif (SEQ == OGSEsin)
    difftime = 3/8*OGSEPER;
  else
    disp('error in seqdifftime');
    stop
  end