%% plot_solution
x=[domain1.left: dx/2: domain1.right];
y1=[domain1.left: dy/2: domain1.right];
[meshes.X1_all, meshes.Y1_all]=meshgrid(x, y1);

y2=[domain2.bottom: dy/2: domain2.top];
[meshes.X2_all,meshes.Y2_all]=meshgrid(x,y2);

[meshes.X1_U, meshes.Y1_U]=meshgrid(x(1:2:length(x)), y1(2:2:length(y1)));
[meshes.X1_V, meshes.Y1_V]=meshgrid(x(2:2:length(x)), y1(1:2:length(y1)));
[meshes.X1_P, meshes.Y1_P]=meshgrid(x(2:2:length(x)), y1(2:2:length(y1)));
[meshes.X2_U, meshes.Y2_U]=meshgrid(x(1:2:length(x)), y2(2:2:length(y1)));
[meshes.X2_V, meshes.Y2_V]=meshgrid(x(2:2:length(x)), y2(1:2:length(y1)));
[meshes.X2_P, meshes.Y2_P]=meshgrid(x(2:2:length(x)), y2(2:2:length(y1)));

pde.U1_in=reshape(u1(1: meshes.DOF_u), [m-1, n]);
pde.U1_in=pde.U1_in.';
pde.V1_in=reshape(u1((meshes.DOF_u+1):(meshes.DOF_u+meshes.DOF_v) ), [m, n-1]);
pde.V1_in=pde.V1_in.';
pde.P1=reshape(u1((meshes.DOF_u+meshes.DOF_v+1):end), [m, n]);
pde.P1=pde.P1.';

pde.U2_in=reshape(u2(1: m*n-n), [m-1, n]);
pde.U2_in=pde.U2_in.';
pde.V2_in=reshape(u2((m*n-n+1):(2*n*m-n-m) ), [m, n-1]);
pde.V2_in=pde.V2_in.';
pde.P2=reshape(u2((2*n*m-n-m+1):end), [m, n]);
pde.P2=pde.P2.';

pde.U1=solution.u11(meshes.X1_U, meshes.Y1_U);
pde.U1(1:n, 2:m)=pde.U1_in;
pde.V1=solution.u12(meshes.X1_V, meshes.Y1_V);
pde.V1(2:n, 1:m)=pde.V1_in;

pde.U2=solution.u21(meshes.X2_U, meshes.Y2_U);
pde.U2(1:n, 2:m)=pde.U2_in;
pde.V2=solution.u22(meshes.X2_V, meshes.Y2_V);
pde.V2(2:n, 1:m)=pde.V2_in;

% %% plot_meshes;
% figure('name', 'meshes')
% hold on
% plot(meshes.X1_P,meshes.Y1_P, 'go');
% plot(meshes.X1_U,meshes.Y1_U, 'k+')
% plot(meshes.X1_V,meshes.Y1_V, 'rx')
% hold off

%% 
figure('name', 'numerical solution')
% Velocity u
subplot(2,3,1)
surf(meshes.X1_U, meshes.Y1_U, pde.U1);
title('velocity u (in x - direction)')

subplot(2,3,4)
surf(meshes.X2_U, meshes.Y2_U, pde.U2)

% Velocity v
subplot(2,3,2)
surf(meshes.X1_V, meshes.Y1_V, pde.V1);
title('velocity v (in y - direction)')

subplot(2,3,5)
surf(meshes.X2_V, meshes.Y2_V, pde.V2);

% Pressure p
subplot(2,3,3)
surf(meshes.X1_P,meshes.Y1_P,pde.P1);
title('pressure p')
subplot(2,3,6)
surf(meshes.X2_P, meshes.Y2_P, pde.P2);

% figure(2)
% subplot(2,3,5)
% surf(meshes.X1_all,meshes.Y1_all, sqrt(pde.U1_all.^2+pde.V1_all.^2));
% hold on
% surf(meshes.X2_all,meshes.Y2_all, sqrt(pde.U2_all.^2+pde.V2_all.^2));
% hold off
% title('velocity magnitude |v|')


%% -------- Flow with arrows
    % interpolating with splines to achieve values on the gridpoints of 
    % [meshes.calculated.XView,meshes.calculated.YView] 
figure('name', 'velocity flow')
pde.U1_all = interp2(meshes.X1_U, meshes.Y1_U, full(pde.U1), meshes.X1_all, meshes.Y1_all, 'spline');
pde.V1_all = interp2(meshes.X1_V, meshes.Y1_V, full(pde.V1), meshes.X1_all, meshes.Y1_all, 'spline');
quiver(meshes.X1_all, meshes.Y1_all, pde.U1_all, pde.V1_all);
hold on
pde.U2_all = interp2(meshes.X2_U, meshes.Y2_U, full(pde.U2), meshes.X2_all, meshes.Y2_all);
pde.V2_all = interp2(meshes.X2_V, meshes.Y2_V, full(pde.V2), meshes.X2_all, meshes.Y2_all);
quiver(meshes.X2_all, meshes.Y2_all, pde.U2_all, pde.V2_all);
hold off
title('flow')
%% ---------plot contour
figure('name', 'velocity contour')
subplot(1,2,1)
[C, h] = contour(meshes.X1_all, meshes.Y1_all, (pde.U1_all.^2+pde.V1_all.^2).^(0.5));
hold on
clabel(C, h)
[C, h] = contour(meshes.X2_all, meshes.Y2_all, (pde.U2_all.^2+pde.V2_all.^2).^(0.5));
hold off
clabel(C, h)

subplot(1,2,2)
contour(meshes.X1_P,meshes.Y1_P,pde.P1);
hold on
contour(meshes.X2_P,meshes.Y2_P,pde.P2);
hold off

%% plot error of solution
[~, iter] = size(u1_store);
solution.U11=solution.u11(meshes.X1_U,meshes.Y1_U);
solution.U12=solution.u12(meshes.X1_V,meshes.Y1_V);
solution.U21=solution.u21(meshes.X2_U, meshes.Y2_U);
solution.U22=solution.u22(meshes.X2_V, meshes.Y2_V);
solution.P1=solution.p1(meshes.X1_P, meshes.Y1_P);
solution.P2=solution.p2(meshes.X2_P, meshes.Y2_P);

figure('name', 'absolute error')
subplot(2,3,1)
surf(meshes.X1_U, meshes.Y1_U, solution.U11-pde.U1);
title('error of U1 (in horizontal)')
subplot(2,3,2)
surf(meshes.X1_V, meshes.Y1_V, solution.U12-pde.V1);
title('error of V1 (in vertical)')
subplot(2,3,3)
surf(meshes.X1_P, meshes.Y1_P, solution.P1-pde.P1);
title('error of P1')

subplot(2,3,4)
surf(meshes.X2_U, meshes.Y2_U, solution.U21-pde.U2);
title('error of U2')
subplot(2,3,5)
surf(meshes.X2_V, meshes.Y2_V, solution.U22-pde.V2);
title('error of V2')
subplot(2,3,6)
surf(meshes.X2_P, meshes.Y2_P, solution.P2-pde.P2);
title('error of P2')

%% plot L2 error for 
x=[(domain1.left+dx):dx:(domain1.right-dx)];
y=[(domain1.bottom+dy/2):dy:(domain1.top-dy/2)];
[X,Y]=meshgrid(x,y);
U11_exact = solution.u11(X,Y).';
u11_exact = reshape(U11_exact, [meshes.DOF_u, 1]);

x=[(domain1.left+dx/2):dx:(domain1.right-dx/2)];
y=[(domain1.bottom+dy):dy:(domain1.top-dy)];
[X,Y] = meshgrid(x,y);
U12_exact = solution.u12(X,Y).';
u12_exact = reshape(U12_exact, [meshes.DOF_v, 1]);

P1_exact=solution.P1.';
p1_exact = reshape(P1_exact, [meshes.DOF_p, 1]);

x=[(domain2.left+dx):dx:(domain2.right-dx)];
y=[(domain2.bottom+dy/2):dy:(domain2.top-dy/2)];
[X,Y]=meshgrid(x,y);
U21_exact = solution.u21(X,Y).';
u21_exact = reshape(U21_exact, [meshes.DOF_u, 1]);

x=[(domain2.left+dx/2):dx:(domain2.right-dx/2)];
y=[(domain2.bottom+dy):dy:(domain2.top-dy)];
[X,Y] = meshgrid(x,y);
U22_exact = solution.u22(X,Y).';
u22_exact = reshape(U22_exact, [meshes.DOF_v, 1]);

P2_exact=solution.P2.';
p2_exact = reshape(P2_exact, [meshes.DOF_p, 1]);

for i=1:iter

    err1_u = u11_exact - u1_store(1:meshes.DOF_u, i);
    err1_v = u12_exact - u1_store(meshes.DOF_u+1:meshes.DOF_u+meshes.DOF_v, i);
    err1_p = p1_exact - u1_store(meshes.DOF_u+meshes.DOF_v+1:end, i);
    
    err1_U(i) = sqrt((sum(err1_u.^2)+sum(err1_v.^2))*dx*dy*0.5);
    err1_P(i) = sqrt(sum(err1_p.^2)*dx*dy);
    
%     err2_u = err2_all(1:meshes.DOF_u);
%     err2_v = err2_all(meshes.DOF_u+1:meshes.DOF_u+meshes.DOF_v);
%     err2_p = err2_all(meshes.DOF_u+meshes.DOF_v+1:end);
    err2_u = u21_exact - u2_store(1:meshes.DOF_u, i);
    err2_v = u22_exact - u2_store(meshes.DOF_u+1:meshes.DOF_u+meshes.DOF_v, i);
    err2_p = p2_exact - u2_store(meshes.DOF_u+meshes.DOF_v+1:end, i);
  
    err2_U(i) = sqrt((sum(err2_u.^2)+sum(err2_v.^2))*dx*dy*0.5);
    err2_P(i) = sqrt(sum(err2_p.^2)*dx*dy);
    
end

figure()
semilogy(err1_U, 'k-o')
hold on
semilogy(err1_P, 'b-*')
semilogy(err2_U, 'k-.')
semilogy(err2_P, 'b-.')
k=1:100;
plot(k, (0.96).^(k), 'r-')
hold off
shg
disp('the last error')
err1_U(iter)
err1_P(iter)
err2_U(iter)
err2_P(iter)
% computing the convergence rate
rate_1U = err1_U(iter)^(1/iter);
rate_1P = err1_P(iter)^(1/iter);
rate_2U = err2_U(iter)^(1/iter);
rate_2P = err2_P(iter)^(1/iter);

%%  -------The error between the current and the previous solutions
for i=1:iter
    
err1_all=u1_store(:, end)-u1_store(:, i);
err2_all=u2_store(:, end)-u2_store(:, i);
err1_u = err1_all(1:meshes.DOF_u, 1);
err1_v = err1_all(meshes.DOF_u+1:meshes.DOF_u+meshes.DOF_v, 1);
err1_p = err1_all(meshes.DOF_u+meshes.DOF_v+1:end, 1);

err1_U(i) = sqrt((sum(err1_u.^2)+sum(err1_v.^2))*dx*dy*0.5);
err1_P(i) = sqrt(sum(err1_p.^2)*dx*dy);

err2_u = err2_all(1:meshes.DOF_u);
err2_v = err2_all(meshes.DOF_u+1:meshes.DOF_u+meshes.DOF_v);
err2_p = err2_all(meshes.DOF_u+meshes.DOF_v+1:end);

err2_U(i) = sqrt((sum(err2_u.^2)+sum(err2_v.^2))*dx*dy*0.5);
err2_P(i) = sqrt(sum(err2_p.^2)*dx*dy);

end
figure('name', 'L2 error test version')
semilogy(err1_U, 'k-o', 'displayname', 'd1-velocity')
hold on
semilogy(err1_P, 'b-*', 'displayname', 'd1-pressure')
semilogy(err2_U, 'k-.', 'displayname', 'd2-velocity')
semilogy(err2_P, 'b-.', 'displayname', 'd2-pressure')
k=1:10;
plot(k, (10^(-1)).^(k), 'r-', 'displayname', '0.1^k')
hold off
legend(gca, 'show')
shg

% ------test
% figure()
% plot(u1_store(:, 1));
% hold on
% for i=2:5
%     plot(u1_store(:, i));
% end
% hold off

    
%% plot exact solutions
x_test=linspace(domain1.left,domain1.right,m);
y_test=linspace(domain1.bottom,domain1.top,n);
[X_test,Y_test] = meshgrid(x_test,y_test);
figure('name', 'exact_solution')
subplot(2,3,1)
surf(X_test,Y_test,solution.u11(X_test,Y_test));
subplot(2,3,2)
surf(X_test, Y_test, solution.u12(X_test,Y_test));
subplot(2,3,3)
surf(X_test, Y_test, solution.p1(X_test, Y_test))

x_test=linspace(domain2.left,domain2.right,m);
y_test=linspace(domain2.bottom,domain2.top,n);
[X_test,Y_test] = meshgrid(x_test,y_test);

subplot(2,3,4)
surf(X_test, Y_test, solution.u21(X_test, Y_test));
subplot(2,3,5)
surf(X_test, Y_test, solution.u22(X_test, Y_test));
subplot(2,3,6)
surf(X_test, Y_test, solution.p2(X_test, Y_test))


shg


