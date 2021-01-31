
paras;
nu1 = Paras.nu1; % the dynamic viscosity. for the PDEs.
nu2 = Paras.nu2;
kappa = Paras.kappa; % drag coefficients. for interface conditions.
sigma = Paras.sigma; % the reaction term. for the PDEs. \sigma in paper.

%% The constructed solutions
% a=1;
% solution.u11=@(x,y) -a*x.^2.*(x-1).^2.*(y-1);
% solution.u12=@(x,y)  a*x.*y.*(6*x+y-3*x.*y+2*x.^2.*y-4*x.^2-2);
% solution.p1=@(x,y) cos(pi*x).*sin(pi*y);
% solution.f11=@(x,y) -sigma*a*x.^2.*(x-1).^2.*(y-1) - nu1 * ( -2*a*(x-1).^2.*(y-1) - 8*a*x.*(x-1).*(y-1) - ...
%                                 2*a*x.^2.*(y-1) ) - pi*sin(pi*x).*sin(pi*y);
% solution.f12=@(x,y) sigma*a*x.*y.*(2*x.^2.*y-4*x.^2-3.*x.*y+6*x+y-2)-nu1*(6*a*y.*(2*x-1).*(y-2)+ ...
%                                 2*a*x.*(2*x.^2-3*x+1))+pi*cos(pi*x).*cos(pi*y);                             
% solution.p2=@(x,y) cos(pi*x).*sin(pi*y);
% solution.u21=@(x,y) -nu1 / nu2 .* a .* x .^ 2 .* (x - 1) .^ 2 .* y + (nu1 ./ kappa + 1) .* a .* x .^ 2 .* (x - 1) .^ 2;
% solution.u22=@(x,y) (2 .* x - 1) .* nu1 ./ nu2 .* a .* x .* (x - 1) .* y .^ 2 - (2 .* x - 1) .* (2 .* nu1 ./ kappa + 2) .* a .* x .* (x - 1) .* y;
% solution.f21=@(x,y) (sigma .* (-nu1 ./ nu2 .* a .* x .^ 2 .* (x - 1) .^ 2 .* y + (nu1 ./ kappa + 1) .* a .* x .^ 2 .* (x - 1) .^ 2)) - (nu2 .* (-2 .* nu1 / nu2 .* a .* (x - 1) .^ 2 .* y - 8 .* nu1 ./ nu2 .* a .* x .* (x - 1) .* y - 2 .* nu1 ./ nu2 .* a .* x .^ 2 .* y + 2 .* (nu1 ./ kappa + 1) .* a .* (x - 1) .^ 2 + 8 .* (nu1 ./ kappa + 1) .* a .* x .* (x - 1) + 2 .* (nu1 ./ kappa + 1) .* a .* x .^ 2)) - pi .* sin(pi .* x) .* sin(pi .* y);
% solution.f22=@(x,y) (sigma .* ((2 .* x - 1) .* nu1 ./ nu2 .* a .* x .* (x - 1) .* y .^ 2 - (2 .* x - 1) .* (2 .* nu1 ./ kappa + 2) .* a .* x .* (x - 1) .* y)) - (nu2 .* (4 .* nu1 ./ nu2 .* a .* (x - 1) .* y .^ 2 + 4 .* nu1 ./ nu2 .* a .* x .* y .^ 2 + 2 .* (2 .* x - 1) .* nu1 ./ nu2 .* a .* y .^ 2 - 4 .* (2 .* nu1 ./ kappa + 2) .* a .* (x - 1) .* y - 4 .* (2 .* nu1 ./ kappa + 2) .* a .* x .* y - 2 .* (2 .* x - 1) .* (2 .* nu1 ./ kappa + 2) .* a .* y + 2 .* (2 .* x - 1) .* nu1 ./ nu2 .* a .* x .* (x - 1))) + cos(pi .* x) .* pi .* cos(pi .* y);

%%  ----TEST-----
solution.u11=@(x,y) zeros(size(x)).*zeros(size(y));
solution.u12=@(x,y) zeros(size(x)).*zeros(size(y));
solution.u21=@(x,y) zeros(size(x)).*zeros(size(y));
solution.u22=@(x,y) zeros(size(x)).*zeros(size(y));
solution.p1=@(x,y) zeros(size(x)).*zeros(size(y));
solution.f11=@(x,y) zeros(size(x)).*zeros(size(y));
solution.f12=@(x,y) zeros(size(x)).*zeros(size(y));
solution.p2=@(x,y) zeros(size(x)).*zeros(size(y));
solution.f21=@(x,y) zeros(size(x)).*zeros(size(y));
solution.f22=@(x,y) zeros(size(x)).*zeros(size(y));

%% -------plot---------
% m=20;
% n=20;
% x_test=linspace(0,10,m);
% y_test=linspace(-1,1,n);
% [X_test,Y_test] = meshgrid(x_test,y_test);
% 
% figure('name', 'exact_solution')
% subplot(2,3,1)
% surf(X_test,Y_test,solution.u11(X_test,Y_test));
% subplot(2,3,2)
% surf(X_test, Y_test, solution.u12(X_test,Y_test));
% subplot(2,3,3)
% surf(X_test, Y_test, solution.p1(X_test, Y_test));
% % 
% y_test=-y_test;
% [X_test,Y_test] = meshgrid(x_test,y_test);
% 
% subplot(2,3,4)
% surf(X_test, Y_test, solution.u21(X_test, Y_test));
% subplot(2,3,5)
% surf(X_test, Y_test, solution.u22(X_test, Y_test));
% subplot(2,3,6)
% surf(X_test, Y_test, solution.p2(X_test, Y_test))
% 
% figure()
% surf(X_test,Y_test,solution.f11(X_test,Y_test));
% figure()
% surf(X_test,Y_test,solution.f12(X_test,Y_test));

% figure()
% surf(X,Y,f21(X,Y));
% figure()
% surf(X,Y,solution.f22(X,Y));

% shg

