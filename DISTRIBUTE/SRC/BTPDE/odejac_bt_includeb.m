function Yout= odejac_bt_includeb(t,Y)

  global FEM_K FEM_A FEM_Q 
  global QVAL
  
  %Yout = -FEM_M_U\FEM_M_L\(FEM_K+FEM_A*seqprofile(t)+FEM_Q);
  Yout = -(FEM_K+FEM_A*QVAL*seqprofile(t)+FEM_Q);
  
  %Yout = FEM_JAC;
  %%Yout = pcg(FEM_M,-(FEM_K+FEM_A+FEM_Q),1e-6,150,[],[]);
  