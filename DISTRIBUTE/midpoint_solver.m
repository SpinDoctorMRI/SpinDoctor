function sol=midpoint_solver(TLIST,ICC,options)
global FEM_M FEM_K FEM_A FEM_Q QVAL

tp = TLIST(1); % previous time
dt = options.dt;
t = tp; % current time
u = ICC; 
tol = options.tol; maxit = 10000;
u_tarray = [u];
tarray = [tp];
KQ = -FEM_K - FEM_Q;
while t < TLIST(2) + dt
    Af_t  = (KQ - QVAL*seqprofile(t )*FEM_A );
    Af_tp = (KQ - QVAL*seqprofile(tp)*FEM_A );
    Al = 1/dt*FEM_M - 0.5*Af_t;
    Ar = 1/dt*FEM_M + 0.5*Af_tp;
    [u,flag,relres,iter,resvec] = bicgstab(Al,Ar*u, tol, maxit, [], [], u);
    u_tarray = [u_tarray, u];
    tarray = [tarray, t];
    tp = t;
    t = tp + dt;
end;
sol.y = u_tarray;
sol.x = tarray;
