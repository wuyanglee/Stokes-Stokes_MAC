%% drive sigma (tau)

clear all

sigmaS=[0, 0.1, 1, 10, 100];

paras;
nu1 = Paras.nu1;
nu2 = Paras.nu2;
kappa = Paras.kappa;
kmin = pi;



index=0;
rho_analyticS = zeros(1, length(sigmaS));

for sigma=sigmaS
    Paras.sigma = sigma;
    rho_analytic = sqrt ( kappa^2 / ( (kappa+ sigma/( sqrt( kmin^2+sigma/nu1 ) -kmin ) )...
        *(kappa+sigma/ (sqrt(kmin^2 +sigma/nu2)-kmin ))));
    index=index+1;
    rho_analyticS(index) = rho_analytic;
    main();
    
end
%%
folder='sigma/';
baseMatFileName= sprintf('rho_analyticS_%d.csv', index);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
csvwrite(fullMatFileName, rho_analyticS);

%%  ------- Load saved data

close all

load sigma/rho_analyticS_5.csv
rho_analyticS = rho_analyticS_5;

load sigma/sigma100.00_err1_U.csv
load sigma/sigma10.00_err1_U.csv
load sigma/sigma1.00_err1_U.csv
load sigma/sigma0.10_err1_U.csv
load sigma/sigma0.00_err1_U.csv

nu1=0.1;
nu2=0.2;
kappa=1;

%  Plot
kmin=pi;

rho_limit = sqrt(kappa^2 / (kappa+2*nu1*kmin) / (kappa+2*nu2*kmin));

% rho_analytic = sqrt ( kappa^2 / ( (kappa+ sigma/( sqrt( kmin^2+sigma/nu1 ) -kmin ) )...
%     *(kappa+sigma/ (sqrt(kmin^2 +sigma/nu2)-kmin ))));

% k=1:iter;  % for current data

figure()
k=1:length(sigma0_10_err1_U);  % for saved data
fig1 = semilogy(k, 0.05*( rho_limit ).^(k), 'd-.r', 'LineWidth', 2, 'MarkerSize',8);    
hold on
fig2 = semilogy(sigma0_00_err1_U, 'r-.o', 'LineWidth', 2, 'MarkerSize',8);

k=1:length(sigma0_10_err1_U);  % for saved data
fig1 = semilogy(k, 0.05*( rho_analyticS(2) ).^(k), 'g:', 'LineWidth', 2);
hold on
fig2 = semilogy(sigma0_10_err1_U, 'g:^', 'LineWidth', 2, 'MarkerSize',8);


k=1:length(sigma1_00_err1_U);  % for saved data
fig1 = semilogy(k, 0.05*( rho_analyticS(3) ).^(k), '-m', 'LineWidth', 2, 'MarkerSize',8);
hold on
fig2 = semilogy(sigma1_00_err1_U, 'm-*', 'LineWidth', 2, 'MarkerSize',8);

k=1:length(sigma10_00_err1_U);  % for saved data
semilogy(k, 0.05*( rho_analyticS(4) ).^(k), '-.b', 'LineWidth', 2)
semilogy(sigma10_00_err1_U, 'b-.s', 'LineWidth', 2, 'MarkerSize',8)


k=1:length(sigma100_00_err1_U);  % for saved data
semilogy(k, 0.02* ( rho_analyticS(5) ).^(k), '--k', 'LineWidth', 2)
semilogy(sigma100_00_err1_U, 'k--x', 'LineWidth', 2, 'MarkerSize',8)

% semilogy(k, 0.05*rho_limit.^k )



% semilogy(sigma0_10_err1_U, ':o', 'LineWidth', 2)

% xlim([1, length(sigma1_00_err1_U)])
% hold off

leg1 = legend( '$~$', '~~~~0', '$~$', '~~~~0.1', '$~$', '~~~~1','$~$', '~~~~10','$~$', '~~~~100');
% leg1.NumColumns = 2;
% [leg1, ~,~,~ ]= columnlegend(2, { '$\rho^*$', '$\sigma=1$','$\rho^*$', '$\sigma=10$','$\rho^*$', '$\sigma=100$'}); 
% leg1 = gridLegend(fig1, 1,'$\rho^*$', 'Orientation','Horizontal');
title_leg= strcat('$~~\rho^*~~~Error~~~ \sigma~~~~$');
title(leg1, title_leg, 'FontSize', 17)

set(leg1, 'Interpreter', 'latex')
set(leg1, 'FontSize', 17)
set(leg1, 'Orientation','horizontal')
set(leg1, 'NumColumns', 2)

xlabel('Number of iteration', 'FontSize', 17)
ylabel('$\left\Vert\mathbf{E}_{1}\right\Vert_{2}$', 'Interpreter', 'latex', 'FontSize', 17)
set(gca,'FontSize',17)
hold off
shg

%%
figure
semilogy(sigma10_00_err1_U(1:2:end))
shg
hold on
semilogy(sigma10_00_err1_U(2:2:end))
shg

%% Save Data
% folder='sigma/';
% baseMatFileName= sprintf('sigma%d_high%d_err1_U.csv', sigma, high);
% fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
% csvwrite(fullMatFileName, err1_U);
%
% baseMatFileName= sprintf('sigma%d_err1_P.csv', 1/dx);
% fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
% csvwrite(fullMatFileName,err1_P);
%
% baseMatFileName= sprintf('sigma%d_err2_U.csv', 1/dx);
% fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
% csvwrite(fullMatFileName,err2_U);
%
% baseMatFileName= sprintf('sigma%d_err2_P.csv', 1/dx);
% fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
% csvwrite(fullMatFileName,err2_P);