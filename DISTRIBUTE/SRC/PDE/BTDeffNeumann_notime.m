function coeff = BTfreqterm_notime_nob(region,state)

  global UG

  
  coeff = (region.nx*UG(1) + region.ny*UG(2)+region.nz*UG(3));


  