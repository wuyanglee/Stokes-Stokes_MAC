% plot_sigma(tau)

%% reshape the data
x=[(Domain1.left+dx):dx:(Domain1.right-dx)];
y=[(Domain1.bottom+dy/2):dy:(Domain1.top-dy/2)];
[X,Y]=meshgrid(x,y);
U11_exact = solution.u11(X,Y).';
u11_exact = reshape(U11_exact, [Meshes.DOF_u, 1]);

x=[(Domain1.left+dx/2):dx:(Domain1.right-dx/2)];
y=[(Domain1.bottom+dy):dy:(Domain1.top-dy)];
[X,Y] = meshgrid(x,y);
U12_exact = solution.u12(X,Y).';
u12_exact = reshape(U12_exact, [Meshes.DOF_v, 1]);

x=[(Domain1.left+dx/2):dx:(Domain1.right-dx/2)];
y=[(Domain1.bottom+dy/2):dy:(Domain1.top-dy/2)];
[X,Y] = meshgrid(x,y);
solution.P1=solution.p1(X, Y);
P1_exact=solution.P1.';
p1_exact = reshape(P1_exact, [Meshes.DOF_p, 1]);

x=[(Domain2.left+dx):dx:(Domain2.right-dx)];
y=[(Domain2.bottom+dy/2):dy:(Domain2.top-dy/2)];
[X,Y]=meshgrid(x,y);
U21_exact = solution.u21(X,Y).';
u21_exact = reshape(U21_exact, [Meshes.DOF_u, 1]);

x=[(Domain2.left+dx/2):dx:(Domain2.right-dx/2)];
y=[(Domain2.bottom+dy):dy:(Domain2.top-dy)];
[X,Y] = meshgrid(x,y);
U22_exact = solution.u22(X,Y).';
u22_exact = reshape(U22_exact, [Meshes.DOF_v, 1]);

x=[(Domain2.left+dx/2):dx:(Domain2.right-dx/2)];
y=[(Domain2.bottom+dy/2):dy:(Domain2.top-dy/2)];
[X,Y] = meshgrid(x,y);
solution.P2=solution.p2(X, Y);
P2_exact=solution.P2.';
p2_exact = reshape(P2_exact, [Meshes.DOF_p, 1]);

%% calculate the error
area1 = (Domain1.right - Domain1.left)*(Domain1.top-Domain1.bottom);
area2 = (Domain2.right - Domain2.left)*(Domain2.top-Domain2.bottom);

[~, iter] = size(u1_store);
err1_U = zeros(1, iter);
err2_U = zeros(1, iter);

for i=1:iter
    err1_u = u11_exact - u1_store(1:Meshes.DOF_u, i);
    err1_v = u12_exact - u1_store(Meshes.DOF_u+1:Meshes.DOF_u+Meshes.DOF_v, i);
    err1_p = p1_exact - u1_store(Meshes.DOF_u+Meshes.DOF_v+1:end, i);

    err1_U(i) = sqrt( (sum(err1_u.^2)+sum(err1_v.^2))*dx*dy + sum(err1_p.^2)*dx*dy);
%     err1_U(i) = err1_U(i) / sqrt( (sum(u11_exact.^2)+sum(u12_exact.^2))*dx*dy*(1/4) + sum(p1_exact.^2)*dx*dy);
%     err1_P(i) = sqrt( sum(err1_p.^2)*dx*dy );

    err2_u = u21_exact - u2_store(1:Meshes.DOF_u, i);
    err2_v = u22_exact - u2_store(Meshes.DOF_u+1:Meshes.DOF_u+Meshes.DOF_v, i);
    err2_p = p2_exact - u2_store(Meshes.DOF_u+Meshes.DOF_v+1:end, i);
  
    err2_U(i) = sqrt( (sum(err2_u.^2)+sum(err2_v.^2) )*dx*dy + sum(err2_p.^2)*dx*dy);
%     err2_U(i) = err2_U(i) / sqrt( (sum(u21_exact.^2)+sum(u22_exact.^2) )*dx*dy*0.25 + sum(p2_exact.^2)*dx*dy);

%     err2_P(i) = sqrt(sum(err2_p.^2)*dx*dy);
end
%% Save data
folder='sigma/';
baseMatFileName= sprintf('sigma%.2f_err1_U.csv', sigma);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
csvwrite(fullMatFileName,err1_U);

folder='sigma/';
baseMatFileName= sprintf('sigma%.2f_err2_U.csv', sigma);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
csvwrite(fullMatFileName, err2_U);

%%  ------- Load saved data 
% load sigma/sigma100.00_err1_U.csv
% load sigma/sigma10.00_err1_U.csv
% load sigma/sigma1.00_err1_U.csv
% load sigma/sigma0.10_err1_U.csv
% nu1=0.1;
% nu2=0.2;
% kappa=1;
% 
% %  Plot
% kmin=pi;
% rho_analytic = sqrt( (kappa^2)/((kmin*nu2+kappa)*(kmin*nu1+kappa)) );
% % k=1:iter;  % for current data
% k=1:length(sigma1_00_err1_U);  % for saved data
% 
% figure()
% semilogy(k, (rho_analytic).^(k), '-k', 'LineWidth', 2)
% hold on
% 
% semilogy(sigma100_00_err1_U, '-o', 'LineWidth', 2)
% semilogy(sigma10_00_err1_U, '--o', 'LineWidth', 2)
% semilogy(sigma1_00_err1_U, '-.o', 'LineWidth', 2)
% semilogy(sigma0_10_err1_U, ':o', 'LineWidth', 2)
% 
% xlim([1, length(sigma1_00_err1_U)])
% % hold off
% leg1 = legend( '$\rho^*$', '$\sigma=100$', '$\sigma=10$', '$\sigma=1$', '$\sigma=0.1$');
% set(leg1, 'Interpreter', 'latex')
% set(leg1, 'FontSize', 17)
% % xlabel('Iteration', 'FontSize', 17)
% % ylabel('$\left\Vert\mathbf{e}_{1,h}^m\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)
% % title()
% set(gca,'FontSize',17)
% shg

%% Save data for kappa
% folder='kappa/';
% baseMatFileName= sprintf('kappa%.2f_err1_U.csv', kappa);
% fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
% csvwrite(fullMatFileName,err1_U);
% 
% folder='kappa/';
% baseMatFileName= sprintf('kappa%.2f_err2_U.csv', kappa);
% fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
% csvwrite(fullMatFileName, err2_U);
%%  Plot for kappa
% load ./kappa/kappa0.10_err1_U.csv
% load ./kappa/kappa1.00_err1_U.csv
% load ./kappa/kappa10.00_err1_U.csv
% load ./kappa/kappa100.00_err1_U.csv
% nu1=0.1;
% nu2=0.05
% % kappa =1;
% 
% % kmin=pi;
% % rho_analytic = sqrt( (kappa^2)/((kmin*nu2+kappa)*(kmin*nu1+kappa)) );
% % k=1:iter;  % for current data
% % k=1:length(kappa0_10_err1_U);  % for saved data
% 
% figure()
% % semilogy(k, (rho_analytic).^(k), '-k', 'LineWidth', 2)
% % hold on
% semilogy(kappa100_00_err1_U, '-.o', 'LineWidth', 2)
% hold on
% semilogy(kappa10_00_err1_U, '--o', 'LineWidth', 2)
% 
% semilogy(kappa1_00_err1_U, '-o', 'LineWidth', 2)
% semilogy(kappa0_10_err1_U, ':o', 'LineWidth', 2)
% 
% xlim([1, 100])
% xlabel('Iteration')
% ylabel('Error')
% % hold off
% leg1 = legend( '$\kappa=100$', '$\kappa=10$', '$\kappa=1$', '$\kappa=0.1$');
% set(leg1, 'Interpreter', 'latex')
% set(leg1, 'FontSize', 17)
% % xlabel('Iteration', 'FontSize', 17)
% % ylabel('$\left\Vert\mathbf{e}_{1,h}^m\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)
% % title()
% set(gca,'FontSize',17)
% shg




