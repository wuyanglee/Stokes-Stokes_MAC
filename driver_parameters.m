% Drive the parameters' main
% You need commend first three lines scripts in the "para" file.
% Change the exact solution.
%%
clear all
close all


nu1S = [0.01, 0.1, 1];
nu2S = [0.01, 0.1, 1];
kappaS = [0.1, 1,  10];
sigmaS = [1, 10, 100];
% kappas = [0.1];
% nu1s=[0.1];
% nu2s=[1];
% paras;
% Paras.sigma = 1;
% sigma = Paras.sigma;
% kmin = pi;

ratioS = [];
relative_ratioS=[];
rate_1US=[];
rho_analyticS=[];
iterS = [];

for kappa=kappaS
    Paras.kappa=kappa;
    
    for nu1=nu1S
        Paras.nu1=nu1;
        for sigma=sigmaS
            Paras.sigma=sigma;
            for nu2=nu2S
                Paras.nu2=nu2;
                if nu2==nu1
                    continue   % jump out this iteration
                end
                kappa
                nu1
                sigma
                nu2
                main;  % Important

                ratioS=[ratioS, ratio];
                %                 relative_ratioS=[relative_ratioS, relative_ratio];
                rho_analyticS=[rho_analyticS, rho_analytic];
                rate_1US=[rate_1US, rate_1U];
                iterS = [iterS, iter];
            end
        end
    end
end


%% Save Data
folder='parameters/';
baseMatFileName= sprintf('parameters_sigma%d_ratios.csv', sigma);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
ratioS = full(ratioS);
csvwrite(fullMatFileName, ratioS);

% folder='parameters/';
% baseMatFileName= sprintf('parameters_sigma%d_relative_ratios.csv', sigma);
% fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
% relative_ratioS = full(relative_ratioS);
% csvwrite(fullMatFileName, relative_ratioS);

baseMatFileName= sprintf('parameters_sigma%d_rhoS.csv', sigma);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
rho_analyticS = full(rho_analyticS);
csvwrite(fullMatFileName, rho_analyticS);

baseMatFileName= sprintf('parameters_sigma%d_rateS.csv', sigma);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
rate_1US = full(rate_1US);
csvwrite(fullMatFileName, rate_1US);

folder='parameters/';
baseMatFileName= sprintf('parameters_sigma%d_iterS.csv', sigma);
fullMatFileName  = fullfile(folder, baseMatFileName); % folder = pwd or wherever...
iterS = full(iterS);
csvwrite(fullMatFileName, iterS);
