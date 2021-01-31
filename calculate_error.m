%% plot L2 error for 
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

x=[(domain1.left+dx):dx:(domain1.right-dx)];
y=[(domain1.bottom+dy/2):dy:(domain1.top-dy/2)];
[X,Y]=meshgrid(x,y);
U11_exact = solution.u11(X,Y).';
u11_exact = reshape(U11_exact, [meshes.DOF_u, 1]); %

x=[(domain1.left+dx/2):dx:(domain1.right-dx/2)];
y=[(domain1.bottom+dy):dy:(domain1.top-dy)];
[X,Y] = meshgrid(x,y);
U12_exact = solution.u12(X,Y).';
u12_exact = reshape(U12_exact, [meshes.DOF_v, 1]); %

solution.P1=solution.p1(meshes.X1_P, meshes.Y1_P);
solution.P2=solution.p2(meshes.X2_P, meshes.Y2_P);

P1_exact=solution.P1.';
p1_exact = reshape(P1_exact, [meshes.DOF_p, 1]); %

x=[(domain2.left+dx):dx:(domain2.right-dx)];
y=[(domain2.bottom+dy/2):dy:(domain2.top-dy/2)];
[X,Y]=meshgrid(x,y);
U21_exact = solution.u21(X,Y).';
u21_exact = reshape(U21_exact, [meshes.DOF_u, 1]); %

x=[(domain2.left+dx/2):dx:(domain2.right-dx/2)];
y=[(domain2.bottom+dy):dy:(domain2.top-dy)];
[X,Y] = meshgrid(x,y);
U22_exact = solution.u22(X,Y).';
u22_exact = reshape(U22_exact, [meshes.DOF_v, 1]); %
 
P2_exact=solution.P2.';
p2_exact = reshape(P2_exact, [meshes.DOF_p, 1]); %

[~, iter] = size(u1_store);

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
k=1:iter;
% plot(k, (0.01).^(k), 'r-')
hold off
title('L2 error wrt iteration')

shg
disp('the last error')
last_err1_U = err1_U(end)
err1_P(iter);
err2_U(iter);
err2_P(iter);
% computing the convergence rate
rate_1U = err1_U(iter)^(1/iter);
rate_1P = err1_P(iter)^(1/iter);
rate_2U = err2_U(iter)^(1/iter);
rate_2P = err2_P(iter)^(1/iter);
