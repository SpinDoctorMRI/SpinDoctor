function Ft = seqintprofile(time);
  % will be fixed late
  global BDELTA SDELTA SEQ OGSEPER
  global PGSE OGSEsin OGSEcos
  
  Ft = zeros(size(time));

  if (SEQ == PGSE)  
    Ft(time >= 0 & time <= SDELTA) = time(time >= 0 & time <= SDELTA);
    Ft(time < BDELTA & time > SDELTA ) = SDELTA ;
    Ft(time >= BDELTA & time <= BDELTA+SDELTA) = SDELTA-(time(time >= BDELTA & time <= BDELTA+SDELTA)-BDELTA);
  
  elseif (SEQ == OGSEsin)
    [ii] = find(time >= 0 & time <= SDELTA);
    %ft(ii) = sin(time(ii)/OGSEPER*2*pi);
    
    Ft(ii) = (1-cos(time(ii)/OGSEPER*2*pi))*OGSEPER/(2*pi);
    tmp = (1-cos(SDELTA/OGSEPER*2*pi))*OGSEPER/(2*pi);
    
    [ii] = find(time >= BDELTA & time <= BDELTA+SDELTA);
    %ft(ii) = -sin((time(ii)-BDELTA)/OGSEPER*2*pi);      
    Ft(ii) = (-1+cos((time(ii)-BDELTA)/OGSEPER*2*pi))*OGSEPER/(2*pi)+tmp;
    
  elseif (SEQ == OGSEcos)
    [ii] = find(time >= 0 & time <= SDELTA);
    %ft(ii) = cos(time(ii)/OGSEPER*2*pi);
    Ft(ii) = sin(time(ii)/OGSEPER*2*pi)*OGSEPER/(2*pi); 
    
    tmp = sin(SDELTA/OGSEPER*2*pi)*OGSEPER/(2*pi);
    
    [ii] = find(time >= BDELTA & time <= BDELTA+SDELTA);
    %ft(ii) = -cos((time(ii)-BDELTA)/OGSEPER*2*pi);  
    Ft(ii) = -sin((time(ii)-BDELTA)/OGSEPER*2*pi)*OGSEPER/(2*pi)+tmp;  
  else
    disp('error in seqprofile');
    stop    
  end

  %disp([time,ft]);