function A=build_A(meshes, Paras, flag)
% STOKESMATRIX constructs a staggered difference matrix for Stokes
%   A=StokesMatrix(nu,m,n); computes a finite difference operator on a
%   staggered grid for the Stokes operator 
%
%     [-nup*Delta           dx]  
%   A=[          -nup*Delta dy]
%     [    dx        dy       ]
%
%   with homogeneous boundary conditions on a rectangular domain of
%   size (0,1)*(0,n*h) where the mesh size h=1/(m+1).  The velocity
%   unknowns (u,v) and pressure unknowns p are stored in lexicographic
%   ordering starting in the y direction, and the staggered grid uses
%   locations u(i*h,(j+1/2)*h), v((i+1/2)*h,j*h) and
%   p((i+1/2)*h,(j+1/2)*h).

%%
% paras;
nu1 = Paras.nu1; % the dynamic viscosity. for the PDEs.
nu2 = Paras.nu2;
% kappa = Paras.kappa; % drag coefficients. for interface conditions.
sigma = Paras.sigma; % the reaction term. for the PDEs. \sigma in paper.


if flag==1
    nu=nu1;
else
    nu=nu2;
end

dx=meshes.dx;
dy=meshes.dy;
m=meshes.m;
n=meshes.n;


e=ones(m-1,1);
Duxx=(1/(dx*dx))*spdiags([e, (-2*e), e],[-1 0 1],m-1,m-1);
e=ones(n,1);
Duyy=(1/(dy*dy))*spdiags([e, (-2*e), e], [-1, 0, 1], n, n);
eye_ux=speye(size(Duxx));
eye_uy=speye(size(Duyy));
Lu=kron(eye_uy,Duxx)+kron(Duyy,eye_ux);
Lu=-nu*Lu+sigma*speye((m-1)*n);

e=ones(m,1);
Dvxx=(1/(dx*dx))*spdiags([e, (-2*e), e],[-1 0 1],m,m);
e=ones(n-1,1);
Dvyy=(1/(dy*dy))*spdiags([e, (-2*e), e], [-1, 0, 1], n-1, n-1);
eye_vx=speye(size(Dvxx));
eye_vy=speye(size(Dvyy));
Lv=kron(eye_vy,Dvxx)+kron(Dvyy,eye_vx);
Lv=-nu*Lv+sigma*speye((n-1)*m);

e=ones(m,1);
Dpx=spdiags([-e, e], [0, 1], m-1,m)./(dx);
eye_py=speye(n);
Lpx=kron(eye_py, Dpx);            

e=ones(n, 1);
Dpy=spdiags([-e, e], [0, 1], n-1, n)./(dy);
eye_px=speye(m);
Lpy=kron(Dpy, eye_px);

e=ones(m, 1);
Dux=spdiags([-e, e], [-1, 0], m, m-1)./dx;
eye_vy=speye(n);
Lux=kron(eye_vy, Dux);

e=ones(n,1);
Dvy=spdiags([-e,e], [-1,0], n, n-1)./dy;
eye_ux=speye(m);
Lvy=kron(Dvy, eye_ux);

A=[      Lu    ,                   sparse((m-1)*n, m*(n-1)),     Lpx;
   sparse((n-1)*m, n*(m-1))   ,    Lv      ,                          Lpy;
           Lux     ,                        Lvy        ,                      sparse(n*m, m*n) ];

end

