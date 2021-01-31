% clear all;
% close all;
% clc;

paras;
nu1 = Paras.nu1; % the dynamic viscosity. for the PDEs.
nu2 = Paras.nu2;
kappa = Paras.kappa; % drag coefficients. for interface conditions.
% sigma = Paras.tau; % the reaction term. for the PDEs. \sigma in paper.

% domain size cannot change in this project.
Domain1.left = 0;
Domain1.right= 1;
Domain1.bottom= 0;
Domain1.top= 1;

Domain2.left = Domain1.left;
Domain2.right=Domain1.right;
Domain2.top=Domain1.bottom;
Domain2.bottom=-Domain1.top;

Meshes.m=64;
m=Meshes.m;

Meshes.dx=(Domain1.right-Domain1.left)/m;
dx=Meshes.dx;

Meshes.dy=Meshes.dx;
dy=Meshes.dy;

Meshes.n = (Domain1.top-Domain1.bottom)/dy;
n=Meshes.n;


Meshes.DOF_all=3*m*n-m-n;
Meshes.DOF_u=(m-1)*n;
Meshes.DOF_v=(n-1)*m;
Meshes.DOF_p=m*n;

exact_solution;

%%
A1=build_A(Meshes, Paras, 1);
f1=build_f(Domain1, Meshes, Paras, 1);
[A1,f1]=modify_boundary_dirichlet(A1, f1 ,Domain1, Meshes, Paras, 1);
A1=[A1; sparse(1,Meshes.DOF_u+Meshes.DOF_v), ones(1,Meshes.DOF_p)]; % pressure constrain, mean==0.
f1=[f1; 0];
f1_origin=f1;
A1_origin=A1;

%% ------- initial interface condition;
x=linspace(Domain1.left, Domain1.right, m-1);
rng(1);  % Control random numbers generator
u2_interface=2*rand(size(x))-1; % random number between (-1, 1)
u2=sparse(Meshes.DOF_u,1);
u2( ((m-1)*(n-2)+1):(m-1)*(n-1), 1 ) = u2_interface;
u2(((m-1)*(n-1)+1):end, 1) = u2_interface;

% % ------- initial interface condition;
% x=linspace(domain1.left, domain1.right, m-1);
% % u2_interface=rand(size(x));
% k=1;
% u2_interface=sin(k*pi*x/(domain1.right-domain1.left));  %% set the frequency of the initial value
% u2=sparse(meshes.DOF_u,1);
% u2( ((m-1)*(n-2)+1):(m-1)*(n-1), 1 ) = u2_interface;
% u2(((m-1)*(n-1)+1):end, 1) = u2_interface;
% % u2=solution.u21([dx:dx:1-dx], 0); % %% ----test----

%% % % %
A2=build_A(Meshes, Paras, 2);
f2=build_f(Domain2, Meshes, Paras, 2);
[A2,f2]=modify_boundary_dirichlet(A2,f2,Domain2, Meshes, Paras, 2);
A2=[A2; sparse(1, 2*m*n-m-n), ones(1, m*n)];  % pressure constrain, mean=0.
f2=[f2; 0];
f2_origin=f2;
A2_origin=A2;

%------- initial interface condition;

x=linspace(Domain2.left, Domain2.right, m-1);
rng(2);  % Control random numbers generator
u1=2*rand(size(x)) -1 ; % random number between (-1, 1)
rng(3);  % Control random numbers generator
u1_2=2*rand(size(x)) -1 ; % for second line. random number between (-1, 1)
u1=[u1, u1_2];

% % ------- initial interface condition;
% x=linspace(domain2.left, domain2.right, m-1);
% % u1=rand(size(x));
% u1=sin(k*pi*x/(domain2.right-domain2.left)); %% Set the frequency of the initial value
% u1=[u1, u1];
% % u1=rand(size(f2));  % initial interface condition;
% % u1=solution.u11([dx:dx:1-dx], 0); % test


%%
epsilon=1e-8;
u1_store=sparse(Meshes.DOF_all, 2);
u2_store=sparse(Meshes.DOF_all, 2);
for i=0:1:500

    [A1, f1] = modify_interface(A1_origin,f1_origin, u2, Meshes, Paras, 1);
    [A2, f2] = modify_interface(A2_origin,f2_origin, u1, Meshes, Paras, 2);

% -----Old version-----
%     [A1,f1]=modify_interface_relaxed_case2(A1_origin,f1_origin, u2, meshes, 1);
%     [A2,f2]=modify_interface_relaxed_case2(A2_origin,f2_origin, u1, meshes, 2);
    u1=A1\f1;
    u2=A2\f2;

    u1_store(:, i+1)=u1;
    u2_store(:, i+1)=u2;

    if i~=0
%         i
%         err_u1=norm((u1_store(:, i)-u1))
%         err_u2=norm((u2_store(:, i)-u2))
        if max(abs(u1_store(:, i) - u1))<epsilon && max(abs(u2_store(:, i) - u2))<epsilon
%         if max(abs(u1))<epsilon && max(abs(u2))<epsilon
%         if norm((u1_store(:, i) - u1))<epsilon && norm(u2_store(:, i) - u2) <epsilon
            disp('=========the iteration: '); i=i+1
            break
        end
    end
end

% test_plot()
% plot_numerical_and_exact()
plot_parameters()
% plot_hmax
% plot_solution
% plot_sigma()

%% save
% filename=['./u1_store_m', num2str(m), '.mat'];
% save(filename, 'u1_store');
% filename=['./u2_store_m', num2str(m), '.mat'];
% save(filename, 'u2_store');






