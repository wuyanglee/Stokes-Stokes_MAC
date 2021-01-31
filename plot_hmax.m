% Plot_hmax
%% Reshape data
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

%% calculate errors
[~, iter] = size(u1_store);
for i=1:iter
    err1_u = u11_exact - u1_store(1:Meshes.DOF_u, i);
    err1_v = u12_exact - u1_store(Meshes.DOF_u+1:Meshes.DOF_u+Meshes.DOF_v, i);
    err1_p = p1_exact - u1_store(Meshes.DOF_u+Meshes.DOF_v+1:end, i);

    err1_U(i) = sqrt( (sum(err1_u.^2)+sum(err1_v.^2))*dx*dy + sum(err1_p.^2)*dx*dy);
%     err1_P(i) = sqrt( sum(err1_p.^2)*dx*dy );

    err2_u = u21_exact - u2_store(1:Meshes.DOF_u, i);
    err2_v = u22_exact - u2_store(Meshes.DOF_u+1:Meshes.DOF_u+Meshes.DOF_v, i);
    err2_p = p2_exact - u2_store(Meshes.DOF_u+Meshes.DOF_v+1:end, i);
  
    err2_U(i) = sqrt( (sum(err2_u.^2)+sum(err2_v.^2) )*dx*dy + sum(err2_p.^2)*dx*dy);
%     err2_P(i) = sqrt(sum(err2_p.^2)*dx*dy);
end

err1_U(iter)


%% computing the convergence rate
% disp('the numerical convegence rate, rho=')
if mod(iter, 2)==1  % If iter is Odd, 
    rate_1U = (err1_U(iter)/err1_U(1))^(1/(iter-1));
    rate_2U =  (err2_U(iter)/err2_U(1))^(1/(iter-1));
else   % IF iter is Even, the first term in vectors err do NOT be used.
    rate_1U = (err1_U(iter-1)/err1_U(1))^(1/(iter-2));
    rate_2U =  (err2_U(iter-1)/err2_U(1))^(1/(iter-2));
end

kmin=pi;

rho_analytic = sqrt ( kappa^2 / ( (kappa+ sigma/( sqrt( kmin^2+sigma/nu1 ) -kmin ) )...
                     *(kappa+sigma/ (sqrt(kmin^2 +sigma/nu2)-kmin ))));

% lambda1 = sqrt(kmin^2 + sigma/nu1);
% lambda2 = sqrt(kmin^2 + sigma/nu2);
% rho_2 = (sqrt((1+sigma/kappa/(lambda1-kmin)) * (1+sigma/kappa/(lambda2-kmin))))^(-1);

            
ratio = rate_1U/rho_analytic;

relative_ratio = abs(rate_1U-rho_analytic)/rho_analytic;


%% Save Data
folder='hmax/';
baseMatFileName= sprintf('hmax%d_err1_U.csv', 1/dx);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
csvwrite(fullMatFileName,err1_U);

baseMatFileName= sprintf('hmax%d_err2_U.csv', 1/dx);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
csvwrite(fullMatFileName,err2_U);

%% Plot hmax
load hmax/hmax128_err1_U.csv
load hmax/hmax64_err1_U.csv
load hmax/hmax32_err1_U.csv
load hmax/hmax16_err1_U.csv
load hmax/hmax8_err1_U.csv
nu1=0.1;
nu2=0.01;
kappa=0.1;
sigma=1;

close all
figure()
kmin=pi;

rho_analytic = sqrt ( kappa^2 / ( (kappa+ sigma/( sqrt( kmin^2+sigma/nu1 ) -kmin ) )...
                     *(kappa+sigma/ (sqrt(kmin^2 +sigma/nu2)-kmin ))));
                 

k=1:length(hmax128_err1_U);
semilogy(k, (rho_analytic).^(k), 'k-', 'LineWidth', 2)
hold on
semilogy(hmax8_err1_U, '--s', 'LineWidth', 2)
semilogy(hmax16_err1_U, '--*', 'LineWidth', 2)
semilogy(hmax32_err1_U, '--<', 'LineWidth', 2)
semilogy(hmax64_err1_U, '--x', 'LineWidth', 2)

semilogy(hmax128_err1_U, '--o', 'LineWidth', 2)
% xlim([1,length(hmax128_err1_U)])

hold off
leg1 = legend('$\rho^*$',  '$h=1/8$','$h=1/16$', '$h=1/32$',  '$h=1/64$', '$h=1/128$');
set(leg1, 'Interpreter', 'latex')
set(leg1, 'FontSize', 17)
set(gca, 'FontSize', 17)


xlabel('Number of iteration', 'FontSize', 17)
ylabel('$\left\Vert\mathbf{E}_{1}\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)

shg

%%  Plot e1, e2, epsilon1, epsilon2
load hmax/hmax64_err1_U.csv
load hmax/hmax64_err2_U.csv
% nu1=0.1;
% nu2=0.01;
% kappa=0.1;
% tau=1;

% close all
figure()
% kmin=pi;
% rho_analytic =sqrt( (kappa^2)/((kmin*nu2+kappa)*(kmin*nu1+kappa)) );
k=1:length(hmax64_err1_U);
semilogy(hmax64_err1_U, 'r--o', 'LineWidth', 2)
hold on
semilogy(hmax64_err2_U, 'b-.+', 'LineWidth', 2) 
hold off
xlim([1, length(hmax64_err1_U)])
leg1 = legend('$\left\Vert \mathbf{E}_1 \right\Vert_{2}$', '$\left\Vert\mathbf{E}_2\right\Vert_{2}$');
set(leg1, 'Interpreter', 'latex')
set(leg1, 'FontSize', 17)
set(gca, 'FontSize', 17)
xlabel('Number of iteration', 'FontSize', 17)
% ylabel('$\left\Vert\mathbf{E}_{1}\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)
shg



