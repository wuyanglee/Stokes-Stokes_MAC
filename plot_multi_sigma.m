%%  ------- Load saved data 
load sigma_domain_1x1/sigma100.00_err1_U.csv
load sigma_domain_1x1/sigma10.00_err1_U.csv
load sigma_domain_1x1/sigma1.00_err1_U.csv
load sigma_domain_1x1/sigma0.10_err1_U.csv
nu1=0.1;
nu2=0.2;
kappa=1;

%  Plot
kmin=pi;
rho_analytic = sqrt( (kappa^2)/((kmin*nu2+kappa)*(kmin*nu1+kappa)) );
% k=1:iter;  % for current data
k=1:length(sigma1_00_err1_U);  % for saved data

figure()
% semilogy(k, (rho_analytic).^(k), '-k', 'LineWidth', 2)

semilogy(sigma100_00_err1_U, '-o', 'LineWidth', 2)
hold on

semilogy(sigma10_00_err1_U, '--o', 'LineWidth', 2)
semilogy(sigma1_00_err1_U, '-.o', 'LineWidth', 2)
semilogy(sigma0_10_err1_U, ':o', 'LineWidth', 2)

xlimit = length(sigma1_00_err1_U);
ylimit = 1e-12;
xlim([1, xlimit])
ylim([ylimit, 1])
% hold off
leg1 = legend( '$\sigma=100$', '$\sigma=10$', '$\sigma=1$', '$\sigma=0.1$');
set(leg1, 'Interpreter', 'latex')
set(leg1, 'FontSize', 17)
% xlabel('Iteration', 'FontSize', 17)
% ylabel('$\left\Vert\mathbf{e}_{1,h}^m\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)
% title()
set(gca,'FontSize',17)
shg
%% 
load sigma_domain_1x10/sigma100.00_err1_U.csv
load sigma_domain_1x10/sigma10.00_err1_U.csv
load sigma_domain_1x10/sigma1.00_err1_U.csv
load sigma_domain_1x10/sigma0.10_err1_U.csv
nu1=0.1;
nu2=0.2;
kappa=1;

%  Plot
kmin=pi;
rho_analytic = sqrt( (kappa^2)/((kmin*nu2+kappa)*(kmin*nu1+kappa)) );
% k=1:iter;  % for current data
k=1:length(sigma1_00_err1_U);  % for saved data

figure()
% semilogy(k, (rho_analytic).^(k), '-k', 'LineWidth', 2)
% hold on

semilogy(sigma100_00_err1_U, '-o', 'LineWidth', 2)
hold on
semilogy(sigma10_00_err1_U, '--o', 'LineWidth', 2)
semilogy(sigma1_00_err1_U, '-.o', 'LineWidth', 2)
semilogy(sigma0_10_err1_U, ':o', 'LineWidth', 2)

% xlim([1, length(sigma1_00_err1_U)])

xlim([1, xlimit])
ylim([ylimit, 1])

% hold off
leg1 = legend('$\sigma=100$', '$\sigma=10$', '$\sigma=1$', '$\sigma=0.1$');
set(leg1, 'Interpreter', 'latex')
set(leg1, 'FontSize', 17)
% xlabel('Iteration', 'FontSize', 17)
% ylabel('$\left\Vert\mathbf{e}_{1,h}^m\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)
% title()
set(gca,'FontSize',17)
shg
%% 
load sigma_domain_1x40/sigma100.00_err1_U.csv
load sigma_domain_1x40/sigma10.00_err1_U.csv
load sigma_domain_1x40/sigma1.00_err1_U.csv
load sigma_domain_1x40/sigma0.10_err1_U.csv
nu1=0.1;
nu2=0.2;
kappa=1;

%  Plot
kmin=pi;
rho_analytic = sqrt( (kappa^2)/((kmin*nu2+kappa)*(kmin*nu1+kappa)) );
% k=1:iter;  % for current data
k=1:length(sigma1_00_err1_U);  % for saved data

figure()
% semilogy(k, (rho_analytic).^(k), '-k', 'LineWidth', 2)
% hold on

semilogy(sigma100_00_err1_U, '-o', 'LineWidth', 2)
hold on
semilogy(sigma10_00_err1_U, '--o', 'LineWidth', 2)
semilogy(sigma1_00_err1_U, '-.o', 'LineWidth', 2)
semilogy(sigma0_10_err1_U, ':o', 'LineWidth', 2)

% xlim([1, length(sigma1_00_err1_U)])
xlim([1, xlimit])
ylim([ylimit, 1])
% hold off
% leg1 = legend( '$\rho^*$', '$\sigma=100$', '$\sigma=10$', '$\sigma=1$', '$\sigma=0.1$');
leg1 = legend( '$\sigma=100$', '$\sigma=10$', '$\sigma=1$', '$\sigma=0.1$');

set(leg1, 'Interpreter', 'latex')
set(leg1, 'FontSize', 17)
% xlabel('Iteration', 'FontSize', 17)
% ylabel('$\left\Vert\mathbf{e}_{1,h}^m\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)
% title()
set(gca,'FontSize',17)
shg
%%
load sigma_domain_1x100/sigma100.00_err1_U.csv
load sigma_domain_1x100/sigma10.00_err1_U.csv
load sigma_domain_1x100/sigma1.00_err1_U.csv
load sigma_domain_1x100/sigma0.10_err1_U.csv
nu1=0.1;
nu2=0.2;
kappa=1;

%  Plot
kmin=pi;
rho_analytic = sqrt( (kappa^2)/((kmin*nu2+kappa)*(kmin*nu1+kappa)) );
% k=1:iter;  % for current data
k=1:length(sigma1_00_err1_U);  % for saved data

figure()
% semilogy(k, (rho_analytic).^(k), '-k', 'LineWidth', 2)
% hold on

semilogy(sigma100_00_err1_U, '-o', 'LineWidth', 2)
hold on
semilogy(sigma10_00_err1_U, '--o', 'LineWidth', 2)
semilogy(sigma1_00_err1_U, '-.o', 'LineWidth', 2)
semilogy(sigma0_10_err1_U, ':o', 'LineWidth', 2)

% xlim([1, length(sigma1_00_err1_U)])
xlim([1, xlimit])
ylim([1e-12, 1])
% hold off
% leg1 = legend( '$\rho^*$', '$\sigma=100$', '$\sigma=10$', '$\sigma=1$', '$\sigma=0.1$');
leg1 = legend( '$\sigma=100$', '$\sigma=10$', '$\sigma=1$', '$\sigma=0.1$');

set(leg1, 'Interpreter', 'latex')
set(leg1, 'FontSize', 17)
% xlabel('Iteration', 'FontSize', 17)
% ylabel('$\left\Vert\mathbf{e}_{1,h}^m\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)
% title()
set(gca,'FontSize',17)
shg