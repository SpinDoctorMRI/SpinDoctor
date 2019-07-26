function coeff = BTDeffNeumann_notime(region,state)

  global UG

  
  coeff = (region.nx*UG(1) + region.ny*UG(2)+region.nz*UG(3));


  