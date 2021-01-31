function [A,f]=modify_interface(A, f, u_star, meshes, Paras, flag)
%     for interface condition;

% paras();
nu1 = Paras.nu1; % the dynamic viscosity. for the PDEs.
nu2 = Paras.nu2;
kappa = Paras.kappa; % drag coefficients. for interface conditions.

if flag==1
    nu=nu1;
else
    nu=nu2;
end

dy=meshes.dy;
m=meshes.m;
n=meshes.n;


for i=1:m-1
    temp=(m-1)*(n-1)+i;
    if flag==1
        A(i,i)=A(i,i)-(nu*(nu-0.5*dy*kappa))/(dy*dy*(nu+0.5*dy*kappa));
        f(i)=f(i)+( 2*kappa*dy*nu * ((3/2)*u_star(temp)-0.5*u_star(temp-m+1)) )/(dy*dy*(2*nu+kappa*dy));

    elseif flag==2
        A(temp, temp)=A(temp, temp)-(nu*(nu-0.5*dy*kappa))/(dy*dy*(nu+0.5*dy*kappa));
        f(temp)=f(temp)+(2*kappa*dy*nu* ((3/2)*u_star(i) - 0.5*u_star(i+m-1)) )/(dy*dy*(2*nu+kappa*dy));

    end
end

