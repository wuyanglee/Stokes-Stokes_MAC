% plot domain I 
pde.u=u1(1:mesh.DOF_u);
pde.v=u1((mesh.DOF_u+1):(mesh.DOF_v+mesh.DOF_u));
pde.p=u1((mesh.DOF_u+mesh.DOF_v+1):end);

pde.U = reshape(pde.u, [m-1,n]).';
pde.V = reshape(pde.v, [m, n-1]).';
pde.P = reshape(pde.p, [m, n]).';

mesh.x_u=[dx:dx:(1-dx)];
mesh.y_u=[dy/2:dy:(1-dy/2)]; % different for domain I or domain II
[mesh.X_u, mesh.Y_u]=meshgrid(mesh.x_u, mesh.y_u);
solution.U = solution.u11(mesh.X_u, mesh.Y_u);

mesh.x_v=[dx/2:dx:(1-dx/2)];
mesh.y_v=[dy:dy:(1-dy)]; % different for domain I or domain II
[mesh.X_v, mesh.Y_v]=meshgrid(mesh.x_v, mesh.y_v);
solution.V = solution.u12(mesh.X_v, mesh.Y_v);

mesh.x_p=[dx/2:dx:(1-dx/2)];
mesh.y_p=[dy/2:dy:(1-dy/2)]; % different for domain I or domain II
[mesh.X_p, mesh.Y_p]=meshgrid(mesh.x_p, mesh.y_p);
solution.P = solution.p1(mesh.X_p, mesh.Y_p);

err_u = sqrt(sum(sum((pde.U-solution.U).^2))*dx*dy)
err_v = sqrt(sum(sum((pde.V-solution.V).^2))*dx*dy)
err_p = sqrt(sum(sum((pde.P-solution.P).^2))*dx*dy)

figure('name', 'error')
subplot(1,3,1)
surf(mesh.X_u, mesh.Y_u, solution.U-pde.U)
subplot(1,3,2)
surf(mesh.X_v, mesh.Y_v, solution.V-pde.V)
subplot(1,3,3)
surf(mesh.X_p, mesh.Y_p, solution.P-pde.P)
shg
figure('name', 'numerical')
subplot(1,3,1)
surf(mesh.X_u, mesh.Y_u, pde.U)
subplot(1,3,2)
surf(mesh.X_v, mesh.Y_v, pde.V)
subplot(1,3,3)
surf(mesh.X_p, mesh.Y_p, pde.P)
shg
figure('name', 'exact')
subplot(1,3,1)
surf(mesh.X_u, mesh.Y_u, solution.U)
subplot(1,3,2)
surf(mesh.X_v, mesh.Y_v, solution.V)
subplot(1,3,3)
surf(mesh.X_p, mesh.Y_p, solution.P)

% -----convert domain I to domain II
u1=u2;
solution.u11=solution.u21;
solution.u12=solution.u22;
solution.p1=solution.p2;

% -----plot_solution_stokes
pde.u=u1(1:mesh.DOF_u);
pde.v=u1((mesh.DOF_u+1):(mesh.DOF_v+mesh.DOF_u));
pde.p=u1((mesh.DOF_u+mesh.DOF_v+1):end);

pde.U = reshape(pde.u, [m-1,n]).';
pde.V = reshape(pde.v, [m, n-1]).';
pde.P = reshape(pde.p, [m, n]).';

mesh.x_u=[dx:dx:(1-dx)];
mesh.y_u=[(-1+dy/2):dy:-dy/2]; % different for domain I or domain II
[mesh.X_u, mesh.Y_u]=meshgrid(mesh.x_u, mesh.y_u);
solution.U = solution.u11(mesh.X_u, mesh.Y_u);

mesh.x_v=[dx/2:dx:(1-dx/2)];
mesh.y_v=[(-1+dy):dy:-dy]; % different for domain I or domain II
[mesh.X_v, mesh.Y_v]=meshgrid(mesh.x_v, mesh.y_v);
solution.V = solution.u12(mesh.X_v, mesh.Y_v);

mesh.x_p=[dx/2:dx:(1-dx/2)];
mesh.y_p=[(-1+dy/2):dy:-dy/2]; % different for domain I or domain II
[mesh.X_p, mesh.Y_p]=meshgrid(mesh.x_p, mesh.y_p);
solution.P = solution.p1(mesh.X_p, mesh.Y_p);

err_u = sqrt(sum(sum((pde.U-solution.U).^2))*dx*dy)
err_v = sqrt(sum(sum((pde.V-solution.V).^2))*dx*dy)
err_p = sqrt(sum(sum((pde.P-solution.P).^2))*dx*dy)

figure('name', 'error')
subplot(1,3,1)
surf(mesh.X_u, mesh.Y_u, solution.U-pde.U)
subplot(1,3,2)
surf(mesh.X_v, mesh.Y_v, solution.V-pde.V)
subplot(1,3,3)
surf(mesh.X_p, mesh.Y_p, solution.P-pde.P)
shg
figure('name', 'numerical')
subplot(1,3,1)
surf(mesh.X_u, mesh.Y_u, pde.U)
subplot(1,3,2)
surf(mesh.X_v, mesh.Y_v, pde.V)
subplot(1,3,3)
surf(mesh.X_p, mesh.Y_p, pde.P)
shg
figure('name', 'exact')
subplot(1,3,1)
surf(mesh.X_u, mesh.Y_u, solution.U)
subplot(1,3,2)
surf(mesh.X_v, mesh.Y_v, solution.V)
subplot(1,3,3)
surf(mesh.X_p, mesh.Y_p, solution.P)
shg