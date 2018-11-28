function ft = seqprofile(time);
  
  global BDELTA SDELTA SEQ OGSEPER
  global PGSE OGSEsin OGSEcos
  PGSE = 1;
  OGSEsin = 2;
  OGSEcos = 3;
    
  ft = zeros(size(time));

  
  if (SEQ == PGSE)  

    ft(time >= 0 & time <= SDELTA) = 1;
    ft(time >= BDELTA & time <= BDELTA+SDELTA) = -1;
  elseif (SEQ == OGSEsin)
    [ii] = find(time >= 0 & time <= SDELTA);
    ft(ii) = sin(time(ii)/OGSEPER*2*pi);
    [ii] = find(time >= BDELTA & time <= BDELTA+SDELTA);
    ft(ii) = -sin((time(ii)-BDELTA)/OGSEPER*2*pi);  
  elseif (SEQ == OGSEcos)
    [ii] = find(time >= 0 & time <= SDELTA);
    ft(ii) = cos(time(ii)/OGSEPER*2*pi);
    [ii] = find(time >= BDELTA & time <= BDELTA+SDELTA);
    ft(ii) = -cos((time(ii)-BDELTA)/OGSEPER*2*pi);  
  else
      
    disp('error in seqprofile');
    stop    
  end

  %disp([time,ft]);