%% Plot numerical solution and exact solution
% close all
% load numerical_exact/numerical_exact_32.mat
%% Reshape the numerical data
x=[Domain1.left: dx/2: Domain1.right];
y1=[Domain1.bottom: dy/2: Domain1.top];
[Meshes.X1_all, Meshes.Y1_all]=meshgrid(x, y1);

y2=[Domain2.bottom: dy/2: Domain2.top];
[Meshes.X2_all,Meshes.Y2_all]=meshgrid(x,y2);

[Meshes.X1_U, Meshes.Y1_U]=meshgrid(x(1:2:length(x)), y1(2:2:length(y1)));
[Meshes.X1_V, Meshes.Y1_V]=meshgrid(x(2:2:length(x)), y1(1:2:length(y1)));
[Meshes.X1_P, Meshes.Y1_P]=meshgrid(x(2:2:length(x)), y1(2:2:length(y1)));
[Meshes.X2_U, Meshes.Y2_U]=meshgrid(x(1:2:length(x)), y2(2:2:length(y1)));
[Meshes.X2_V, Meshes.Y2_V]=meshgrid(x(2:2:length(x)), y2(1:2:length(y1)));
[Meshes.X2_P, Meshes.Y2_P]=meshgrid(x(2:2:length(x)), y2(2:2:length(y1)));

pde.U1_in=reshape(u1(1: Meshes.DOF_u), [m-1, n]);
pde.U1_in=pde.U1_in.';
pde.V1_in=reshape(u1((Meshes.DOF_u+1):(Meshes.DOF_u+Meshes.DOF_v) ), [m, n-1]);
pde.V1_in=pde.V1_in.';
pde.P1=reshape(u1((Meshes.DOF_u+Meshes.DOF_v+1):end), [m, n]);
pde.P1=pde.P1.';

pde.U2_in=reshape(u2(1: m*n-n), [m-1, n]);
pde.U2_in=pde.U2_in.';
pde.V2_in=reshape(u2((m*n-n+1):(2*n*m-n-m) ), [m, n-1]);
pde.V2_in=pde.V2_in.';
pde.P2=reshape(u2((2*n*m-n-m+1):end), [m, n]);
pde.P2=pde.P2.';

pde.U1=solution.u11(Meshes.X1_U, Meshes.Y1_U);
pde.U1(1:n, 2:m)=pde.U1_in;
pde.V1=solution.u12(Meshes.X1_V, Meshes.Y1_V);
pde.V1(2:n, 1:m)=pde.V1_in;

pde.U2=solution.u21(Meshes.X2_U, Meshes.Y2_U);
pde.U2(1:n, 2:m)=pde.U2_in;
pde.V2=solution.u22(Meshes.X2_V, Meshes.Y2_V);
pde.V2(2:n, 1:m)=pde.V2_in;

pde.U1_all = interp2(Meshes.X1_U, Meshes.Y1_U, full(pde.U1), Meshes.X1_all, Meshes.Y1_all, 'spline');
pde.V1_all = interp2(Meshes.X1_V, Meshes.Y1_V, full(pde.V1), Meshes.X1_all, Meshes.Y1_all, 'spline');
pde.U2_all = interp2(Meshes.X2_U, Meshes.Y2_U, full(pde.U2), Meshes.X2_all, Meshes.Y2_all, 'spline');
pde.V2_all = interp2(Meshes.X2_V, Meshes.Y2_V, full(pde.V2), Meshes.X2_all, Meshes.Y2_all, 'spline');

%% Exact solutions
solution.U1=solution.u11(Meshes.X1_U,Meshes.Y1_U);
solution.V1=solution.u12(Meshes.X1_V,Meshes.Y1_V);
solution.U2=solution.u21(Meshes.X2_U, Meshes.Y2_U);
solution.V2=solution.u22(Meshes.X2_V, Meshes.Y2_V);
solution.P1=solution.p1(Meshes.X1_P, Meshes.Y1_P);
solution.P2=solution.p2(Meshes.X2_P, Meshes.Y2_P);
solution.U1_all = interp2(Meshes.X1_U, Meshes.Y1_U, full(solution.U1), Meshes.X1_all, Meshes.Y1_all, 'spline');
solution.V1_all = interp2(Meshes.X1_V, Meshes.Y1_V, full(solution.V1), Meshes.X1_all, Meshes.Y1_all, 'spline');

solution.U2_all = interp2(Meshes.X2_U, Meshes.Y2_U, full(solution.U2), Meshes.X2_all, Meshes.Y2_all, 'spline');
solution.V2_all = interp2(Meshes.X2_V, Meshes.Y2_V, full(solution.V2), Meshes.X2_all, Meshes.Y2_all, 'spline');

% %% plot_meshes;
% figure('name', 'meshes')
% hold on
% plot(meshes.X1_P,meshes.Y1_P, 'go');
% plot(meshes.X1_U,meshes.Y1_U, 'k+')
% plot(meshes.X1_V,meshes.Y1_V, 'rx')
% hold off
%% TEST---plot-numerical solution
% figure('name', 'numerical solution')
% % Velocity u
% subplot(2,3,1)
% surf(Meshes.X1_U, Meshes.Y1_U, pde.U1);
% title('velocity u (in x - direction)')
% hold on 
% % subplot(2,3,4)
% surf(Meshes.X2_U, Meshes.Y2_U, pde.U2)
% hold off
% % Velocity v
% subplot(2,3,2)
% surf(Meshes.X1_V, Meshes.Y1_V, pde.V1);
% title('velocity v (in y - direction)')
% hold on
% % subplot(2,3,5)
% surf(Meshes.X2_V, Meshes.Y2_V, pde.V2);
% hold off
% % Pressure p
% subplot(2,3,3)
% surf(Meshes.X1_P,Meshes.Y1_P,pde.P1);
% title('pressure p')
% hold on
% % subplot(2,3,6)
% surf(Meshes.X2_P, Meshes.Y2_P, pde.P2);
% hold off
% 
% figure()
% surf(x, y1, sqrt(solution.U1_all.^2+solution.V1_all.^2) ) 
% hold on
% surf(x, y2, sqrt(solution.U2_all.^2+solution.V2_all.^2) ) 
% hold off
% shg

%% --middle values

close all
nu1=0.1; % the dynamic viscosity. for the PDEs.
nu2=0.2;
kappa=0.1; % drag coefficients. for interface conditions.
tau=0.1; % the reaction term. for the PDEs. \sigma in paper.

y1_middle = Meshes.Y1_P(:, 1);
y2_middle = Meshes.Y2_P(:, 1);
p1_middle = pde.P1(:, 16);
p2_middle = pde.P2(:, 16);
y_p_middle = cat(1, y2_middle, y1_middle);
p_middle = cat(1, p2_middle, p1_middle);

y1_u_middle = Meshes.Y1_all(:, 1);
y2_u_middle = Meshes.Y2_all(:, 1);
u1_middle = pde.U1_all(:, 33);
v1_middle = pde.V1_all(:, 33);
u2_middle = pde.U2_all(:, 33);
v2_middle = pde.V2_all(:, 33);
u1_middle = sqrt(u1_middle.^2 + v1_middle.^2);
u2_middle = sqrt(u2_middle.^2 + v2_middle.^2);

y_u_middle = cat(1, y2_u_middle, y1_u_middle);
u_middle = cat(1, u2_middle, u1_middle);


figure()
plot(y_p_middle, p_middle, '-k', 'LineWidth', 1.5)
xlabel('$y$', 'Interpreter', 'latex', 'FontSize',17)
ylabel('$p$', 'Interpreter', 'latex', 'FontSize',17)
title('Pressure')
set(gca,'FontSize',17)

figure()
plot(y_u_middle, u_middle, '-k',  'LineWidth', 1.5)
xlabel('$y$', 'FontSize',17, 'Interpreter', 'latex')
ylabel('$\mathbf{u}$', 'FontSize',17, 'Interpreter', 'latex')
title('Velocity')
set(gca,'FontSize',17)
shg
%% -------- Velocity --- heatmap + arrows

figure('name', 'velocity flow')
% exact solutions
subplot(1,3,1)
imagesc(x, y1, sqrt(solution.U1_all.^2+solution.V1_all.^2) ) 
hold on
imagesc(x, y2, sqrt(solution.U2_all.^2+solution.V2_all.^2) ) 
% 
quiver(Meshes.X1_all, Meshes.Y1_all, solution.U1_all, solution.V1_all);
quiver(Meshes.X2_all, Meshes.Y2_all, solution.U2_all, solution.V2_all);
hSlices1_stream = streamslice(Meshes.X1_all, Meshes.Y1_all, solution.U1_all, solution.V1_all);
hSlices2_stream = streamslice(Meshes.X2_all, Meshes.Y2_all, solution.U2_all, solution.V2_all);
set(hSlices1_stream, 'LineWidth', 1, 'color', 'k')
set(hSlices2_stream, 'LineWidth', 1, 'color', 'k')
% 
colormap(jet)
% colorbar()
colorbar( 'Ticks', [0 0.1 0.2 0.3 0.4])

axis equal
axis([0,1,-1,1])
axis xy
hold off
title('Exact')
set(gca, 'FontSize', 14)

hold off

% numerical
subplot(1,3,2)
imagesc(x, y1, sqrt(pde.U1_all.^2+pde.V1_all.^2) ) 
hold on
imagesc(x, y2, sqrt(pde.U2_all.^2+pde.V2_all.^2) ) 
% 
hSlices1_stream = streamslice(Meshes.X1_all, Meshes.Y1_all, pde.U1_all, pde.V1_all);
hSlices2_stream = streamslice(Meshes.X2_all, Meshes.Y2_all, pde.U2_all, pde.V2_all);

quiver(Meshes.X1_all, Meshes.Y1_all, pde.U1_all, pde.V1_all);
quiver(Meshes.X2_all, Meshes.Y2_all, pde.U2_all, pde.V2_all);

set(hSlices1_stream, 'LineWidth', 1, 'color', 'k')
set(hSlices2_stream, 'LineWidth', 1, 'color', 'k')
% 
colormap(jet)
colorbar( 'Ticks', [0 0.1 0.2 0.3 0.4])
axis equal
axis([0,1,-1,1])
axis xy
title('Numerical')
set(gca, 'FontSize', 14)

hold off
% 
subplot(1,3,3)
imagesc(x, y1, abs( sqrt(pde.U1_all.^2+pde.V1_all.^2) -sqrt(solution.U1_all.^2+solution.V1_all.^2))) 
hold on
imagesc(x, y2, abs( sqrt(pde.U2_all.^2+pde.V2_all.^2) -sqrt(solution.U2_all.^2+solution.V2_all.^2) )) 

% 
colormap(jet)
colorbar()
axis equal
axis([0,1,-1,1])
axis xy


set(gca, 'FontSize', 14)

hold off
title('Error of velocity')
%%  heatmap pressure
figure()
subplot(1,3,1)
imagesc(x(2:2:length(x)), y1(2:2:length(y1)), solution.P1 ) 
hold on
imagesc(x(2:2:length(x)), y2(2:2:length(y1)), solution.P2 ) 
colormap(jet)
colorbar()
colorbar( 'Ticks', [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9])
axis equal
axis([0,1,-1,1])
axis xy
hold off
set(gca, 'FontSize', 14)

title('Exact')

subplot(1,3,2)
imagesc(x(2:2:length(x)), y1(2:2:length(y1)), pde.P1 ) 
hold on
imagesc(x(2:2:length(x)), y2(2:2:length(y1)), pde.P2 ) 
colormap(jet)
colorbar()
colorbar( 'Ticks', [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9])
axis equal
axis([0,1,-1,1])
axis xy
hold off

title('Numerical')
set(gca, 'FontSize', 14)

subplot(1,3,3)
imagesc(x(2:2:length(x)), y1(2:2:length(y1)), abs(pde.P1 -solution.P1) )
hold on
imagesc(x(2:2:length(x)), y2(2:2:length(y1)), abs(pde.P2 -solution.P2)) 
colormap(jet)
colorbar()
% colorbar( 'Ticks', [0, 0.005, 0.01, 0.015, 0.02, 0.025])
axis equal
axis([0,1,-1,1])
axis xy
hold off
title('Error of pressure')
set(gca, 'FontSize', 14)

%%  --velocity+quiver-- Plot contour 
% figure()
% % exact solutions
% subplot(1,2,1)
% [C, h] = contour(meshes.X1_all, meshes.Y1_all, sqrt(solution.U1_all.^2+solution.V1_all.^2));
% hold on
% clabel(C, h)
% [C, h] = contour(meshes.X2_all, meshes.Y2_all, sqrt(solution.U2_all.^2+solution.V2_all.^2));
% clabel(C, h)
% %
% quiver(meshes.X1_all, meshes.Y1_all, solution.U1_all, solution.V1_all);
% quiver(meshes.X2_all, meshes.Y2_all, solution.U2_all, solution.V2_all);
% hSlices1_stream = streamslice(meshes.X1_all, meshes.Y1_all, solution.U1_all, solution.V1_all);
% hSlices2_stream = streamslice(meshes.X2_all, meshes.Y2_all, solution.U2_all, solution.V2_all);
% set(hSlices1_stream, 'LineWidth', 1, 'color', 'k')
% set(hSlices2_stream, 'LineWidth', 1, 'color', 'k')
% % 
% axis equal
% axis([0,1,-1,1])
% axis xy
% hold off
% title('Velocity')
% % numerical solution
% subplot(1,2,2)
% [C, h] = contour(meshes.X1_all, meshes.Y1_all, sqrt(pde.U1_all.^2+pde.V1_all.^2));
% hold on
% clabel(C, h)
% [C, h] = contour(meshes.X2_all, meshes.Y2_all, sqrt(pde.U2_all.^2+pde.V2_all.^2));
% clabel(C, h)
% %
% hSlices1_stream = streamslice(meshes.X1_all, meshes.Y1_all, pde.U1_all, pde.V1_all);
% hSlices1_quiver = quiver(meshes.X1_all, meshes.Y1_all, pde.U1_all, pde.V1_all);
% hSlices2_stream = streamslice(meshes.X2_all, meshes.Y2_all, pde.U2_all, pde.V2_all);
% hSlices2_quiver = quiver(meshes.X2_all, meshes.Y2_all, pde.U2_all, pde.V2_all);
% set(hSlices1_stream, 'LineWidth', 1, 'color', 'k');
% set(hSlices1_stream, 'LineWidth',1, 'Color', 'k');
% 
% axis equal
% axis([0,1,-1,1])
% axis xy
% hold off
% title('Velocity')

%% ---------plot contour pressure
% figure('name', 'Pressure')
% subplot(1,2,1)
% [C, h] = contour(meshes.X1_P,meshes.Y1_P,solution.P1);
% clabel(C, h)
% hold on
% [C, h] = contour(meshes.X2_P,meshes.Y2_P,solution.P2);
% clabel(C, h)
% hold off
% axis equal
% title('Pressure')
% 
% subplot(1,2,2)
% [C, h] = contour(meshes.X1_P,meshes.Y1_P,pde.P1);
% clabel(C, h)
% hold on
% [C, h] = contour(meshes.X2_P,meshes.Y2_P,pde.P2);
% clabel(C, h)
% hold off
% title('Pressure')


