function Yout= odefun_bt_includeb(t,Y)


  global FEM_K FEM_A FEM_Q FEM_G
  global QVAL
    
  
  Yout = -(FEM_K*Y+FEM_A*Y*QVAL*seqprofile(t)+FEM_Q*Y)+FEM_G*seqintprofile(t);    
  

  %[Yout,flag,relres,iter,resvec] = pcg(FEM_M,Yout,1e-6,150,FEM_M_Prec,FEM_M_Prec');
 
  
  
  %Yout = -(FEM_K*Y+FEM_A*Y*seqprofile(t)+FEM_Q*Y);
  %Yout = FEM_M_U\(FEM_M_L\(Yout));
  %Yout = Yout(:);a