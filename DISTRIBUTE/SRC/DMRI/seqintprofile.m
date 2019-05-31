function Ft = seqintprofile(time);
  % will be fixed late
  global BDELTA SDELTA SEQ OGSEPER
  global PGSE OGSEsin OGSEcos dPGSE
  
  % dPGSE = 4;
  
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
    
  elseif (SEQ == dPGSE)  
    Ft(time >= 0 & time <= SDELTA) = time(time >= 0 & time <= SDELTA);
    Ft(time <= BDELTA & time >= SDELTA ) = SDELTA ;
    Ft(time > BDELTA & time <= BDELTA+SDELTA) = SDELTA-(time(time >= BDELTA & time <= BDELTA+SDELTA)-BDELTA);    
    
    displacement = BDELTA+SDELTA;
    Ft(time > 0 + displacement & time <= SDELTA + displacement) = time(time >= 0 + displacement & time <= SDELTA + displacement) - displacement;
    Ft(time <= BDELTA + displacement & time > SDELTA + displacement ) = SDELTA;
    Ft(time > BDELTA + displacement & time <= BDELTA+SDELTA + displacement) = SDELTA-(time(time >= BDELTA + displacement & time <= BDELTA+SDELTA + displacement)-(BDELTA + displacement));        
  else  
    ft = @(t) seqprofile(t);
    Ft = integral(ft,0,time,'AbsTol',1e-6,'RelTol',1e-3);   
  end