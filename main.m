clear all;
close all;
% clc;

m=4;
n=4;

paras;

% domain size cannot change in this project.
domain1.left = 0;
domain1.right=1;
domain1.bottom=0;
domain1.top=1;

domain2.left = domain1.left;
domain2.right=domain1.right;
domain2.top=domain1.bottom;
domain2.bottom=-domain1.top;

dx=(domain1.right-domain1.left)/m;
dy=(domain1.top-domain1.bottom)/n;

mesh.DOF_all=3*m*n-m-n;
mesh.DOF_u=(m-1)*n;
mesh.DOF_v=(n-1)*m;
mesh.DOF_p=m*n;

exact_solution;

A1=build_A(nu1,tau,dx,dy,m,n);
f1=build_f(dx,dy,m,n,tau,domain1,1);
[A1,f1]=modify_boundary_dirichlet(A1,f1,dx,dy, m,n,nu1,mesh,domain1, 1);
A1=[A1; zeros(1,2*m*n-m-n), ones(1,m*n)]; % pressure constrain, mean=0.
f1=[f1; 0];
f1_origin=f1;
A1_origin=A1;
% initial interface condition;
x=linspace(domain1.left, domain1.right, m-1); 
u2_interface=sin(pi*x/(domain1.right-domain1.left));  %% set the frequency of the initial value 
u2=sparse(mesh.DOF_u,1);
u2( ((m-1)*(n-2)+1):(m-1)*(n-1), 1 ) = u2_interface;
u2(((m-1)*(n-1)+1):end, 1) = u2_interface;
% u2=solution.u21([dx:dx:1-dx], 0); % %% ----test----

% % % % 
A2=build_A(nu2,tau, dx,dy, m,n);
f2=build_f(dx,dy,m,n,tau,domain2,2);
[A2,f2]=modify_boundary_dirichlet(A2,f2,dx,dy, m,n,nu2,mesh,domain2, 2);
A2=[A2; zeros(1, 2*m*n-m-n), ones(1, m*n)];  % pressure constrain, mean=0.
f2=[f2; 0];
f2_origin=f2;
A2_origin=A2;
x=linspace(domain1.left, domain1.right, m-1);
u1=sin(pi*x/(domain1.right-domain1.left)); %% Set the frequency of the initial value 
u1=[u1, u1];
% u1=rand(size(f2));  % initial interface condition;
% u1=solution.u11([dx:dx:1-dx], 0); % test


% % save 
% filename=['./u1_store_m', num2str(m), '.mat'];
% save(filename, 'u1_store');
% filename=['./u2_store_m', num2str(m), '.mat'];
% save(filename, 'u2_store');

% 
epsilon=1e-6;
u1_store=sparse(mesh.DOF_all,1);
u2_store=sparse(mesh.DOF_all,1);
for i=0:1:100
    
    [A1,f1]=modify_interface(A1_origin,f1_origin, u2, nu1, kappa, dx, dy, m, n, 1);
    [A2,f2]=modify_interface(A2_origin,f2_origin, u1, nu2, kappa, dx,dy, m,n, 2);

    u2=A2\f2;
    u1=A1\f1;
    u1_store(:, i+1)=u1;
    u2_store(:, i+1)=u2;
    
    if i~=0
        i
        err_u1=norm((u1_store(:, i)-u1))
        err_u2=norm((u2_store(:, i)-u2))
%         if max(abs(u1_store(:, i) - u1))<epsilon && max(abs(u2_store(:, i) - u2))<epsilon
        if norm((u1_store(:, i) - u1))<epsilon && norm(u2_store(:, i) - u2) <epsilon
            disp('=========the iteration: ')
            i
%             pause
            break
        end
    end
end
% % save 
% filename=['./u1_store_m', num2str(m), '.mat'];
% save(filename, 'u1_store');
% filename=['./u2_store_m', num2str(m), '.mat'];
% save(filename, 'u2_store');


plot_solution





