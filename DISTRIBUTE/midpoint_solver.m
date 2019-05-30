function sol=midpoint_solver(TLIST,ICC,options)
global FEM_M FEM_K FEM_A FEM_Q QVAL

tp = TLIST(1); % previous time
dt = options.dt;
t = tp;
u = ICC; 
tol = options.tol; maxit = 10000;
u_tarray = [u];
tarray = [t];
while t < TLIST(2) + dt
    Al = FEM_M - dt/2*(-FEM_K - QVAL*seqprofile(t) *FEM_A - FEM_Q);
    Ar = FEM_M + dt/2*(-FEM_K - QVAL*seqprofile(tp)*FEM_A - FEM_Q);
    [u,flag,relres,iter,resvec] = bicgstab(Al,Ar*u, tol, maxit, [], [], u);
    u_tarray = [u_tarray, u];
    tp = t;
    t = t + dt;
    tarray = [tarray, t];
end;
sol.y = u_tarray;
sol.x = tarray;

