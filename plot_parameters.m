% Plot_parameters_table
close all
%% === generate exact solution data

%% ==compute errors
[~, iter] = size(u1_store);
err1_U = zeros(iter, 1);
err1_P = zeros(iter, 1);
err2_U = zeros(iter, 1);
err2_P = zeros(iter, 1);
for i=1:iter
    err1_u = u1_store(1:Meshes.DOF_u, i);
    err1_v = u1_store(Meshes.DOF_u+1:Meshes.DOF_u+Meshes.DOF_v, i);
    err1_p = u1_store(Meshes.DOF_u+Meshes.DOF_v+1:end, i);

% numerical test (using the Last numerical solution)
%     err1_u = u1_store(1:meshes.DOF_u, end) - u1_store(1:meshes.DOF_u, i);
%     err1_v = u1_store(meshes.DOF_u+1:meshes.DOF_u+meshes.DOF_v, end)  - u1_store(meshes.DOF_u+1:meshes.DOF_u+meshes.DOF_v, i);
%     err1_p = u1_store(meshes.DOF_u+meshes.DOF_v+1:end, end)  - u1_store(meshes.DOF_u+meshes.DOF_v+1:end, i);
    
    err1_U(i) = sqrt( (sum(err1_u.^2)+sum(err1_v.^2))*dx*dy + sum(err1_p.^2)*dx*dy);

    err2_u = u2_store(1:Meshes.DOF_u, i);
    err2_v = u2_store(Meshes.DOF_u+1:Meshes.DOF_u+Meshes.DOF_v, i);
    err2_p = u2_store(Meshes.DOF_u+Meshes.DOF_v+1:end, i);
  
    err2_U(i) = sqrt( (sum(err2_u.^2)+sum(err2_v.^2) )*dx*dy + sum(err2_p.^2)*dx*dy);
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
folder='parameters/';
baseMatFileName= sprintf('nu1%.2f_nu2%.2f_kappa%.1f_sigma%d_err1_U.csv', nu1, nu2,  kappa, sigma);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
csvwrite(fullMatFileName,full(err1_U));

% baseMatFileName= sprintf('nu1%.2f_err2_U.csv', nu1);
% fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
% csvwrite(fullMatFileName,full(err2_U));


%%  Plot
% figure()
% semilogy(err1_U, 'r-o')
% hold on
% semilogy(err1_P, 'b-o')
% semilogy(err2_U, 'r-*')
% semilogy(err2_P, 'b-*')
% k=1:iter;
% plot(k, (rho_analytic).^(k), 'k-')
% hold off
% shg





