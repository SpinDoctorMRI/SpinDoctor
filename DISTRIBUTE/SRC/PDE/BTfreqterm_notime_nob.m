function coeff = BTfreqterm_notime_nob(region,state)

  global DIFFUSIVITY BDELTA SDELTA UG QVAL

  
  coeff = sqrt(-1)*(region.x*UG(1) + region.y*UG(2)+region.z*UG(3))+1e-16*sqrt(-1);


  