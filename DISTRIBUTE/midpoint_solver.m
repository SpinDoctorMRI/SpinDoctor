function sol=midpoint_solver(TLIST,ICC,options)
global FEM_M FEM_K FEM_A FEM_Q QVAL

tp = TLIST(1); % previous time
dt = options.dt;
t = tp + dt; % current time
u = ICC; 
tol = options.tol; maxit = 10000;
u_tarray = [u];
tarray = [tp];
while t < TLIST(2) + dt
    tm = tp + 0.5*dt;
    Af = 0.5*dt*(-FEM_K - QVAL*seqprofile(tm)*FEM_A - FEM_Q);
    Al = FEM_M - Af;
    Ar = FEM_M + Af;
    [u,flag,relres,iter,resvec] = bicgstab(Al,Ar*u, tol, maxit, [], [], u);
    u_tarray = [u_tarray, u];
    tarray = [tarray, t];
    tp = t;
    t = tp + dt;
end;
sol.y = u_tarray;
sol.x = tarray;

