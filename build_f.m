function [f] = build_f(domain,mesh, Paras, flag)
%BUILD_F Summary of this function goes here
%   Detailed explanation goes here

% paras;
% nu1 = Paras.nu1; % the dynamic viscosity. for the PDEs.
% nu2 = Paras.nu2;
% kappa = Paras.kappa; % drag coefficients. for interface conditions.
% tau = Paras.tau; % the reaction term. for the PDEs. \sigma in paper.


dx=mesh.dx;
dy=mesh.dy;
m=mesh.m;
n=mesh.n;

exact_solution;

if flag==1
    x = [domain.left:dx/2:domain.right];
    y = [domain.bottom:dy/2:domain.top];
    f1=solution.f11;
    f2=solution.f12;
elseif flag==2
    x = [domain.left:dx/2:domain.right];
    y = [domain.bottom:dy/2:domain.top];
    f1=solution.f21;
    f2=solution.f22;
else
    disp('The flag must be 1 or 2. Please check again.')
end

f=sparse(3*m*n-m-n,1);
counter = 1;
for j=1:1:n
    for i=1:1:m-1
        x0=x(i*2);
        x1=x(i*2+2);
        y0=y(j*2-1);
        y1=y(j*2+1);
        f(counter) = integral2(f1, x0, x1, y0, y1)/(dx*dy);
        counter=counter+1;
    end
end
for j=1:1:n-1
    for i=1:1:m
        x0=x(i*2-1);
        x1=x(i*2+1);
        y0=y(j*2);
        y1=y(j*2+2);
        f(counter)=integral2(f2, x0, x1, y0, y1)/(dx*dy);
        counter=counter+1;
    end
end
end

