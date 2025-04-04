clearvars
clc
%% Question 1
rng(42); % set seed to reproduce the same charts    
addpath("utils") % add utils functions

N = 100; % sample size
beta = [1, 2, 1,-1]; % beta coefficients
sigma_u = [1, 2, 10]; % Variance values for U
conflev = [1e-2, 5e-2, 1e-1]; % confidence leves
R = [0, 0, 0, 1]';

% initialize estimation
n_iter = 5; % number of trials
beta_estimates = zeros(length(sigma_u), n_iter, length(beta));
errs_estimates = zeros(length(sigma_u), n_iter);
bounds         = zeros(length(sigma_u), n_iter, length(conflev));

for j=1:length(sigma_u)
    for i=1:n_iter
        % define input matrix 
        X = generate_data(N, [0, 0, 0], [1, 2, 1]);
        U = sample_normal(N, 0, sigma_u(j));
        
        % calculate target variable
        Y = X*beta' + U;
        
        % compute parameter estimates
        beta_hat = X \ Y;
        beta_estimates(j, i, :) = beta_hat;
    
        % compute estimated error (residuals)
        residuals = Y - X*beta_hat;
        s2 = residuals'*residuals./(N-3);
        
        invx = R'*inv(X'*X)*R;
        
        bound = s2*invx;
        %bound = s2 * R'/(X' * X)*R;

        for z=1:length(conflev)
            b = sqrt(bound)*tinv(conflev(z), N-3); %check
            bounds(j, i, z) = b;
        end
        
        errs_estimates(j, i, :) = sqrt(s2);
    end
end

% calculating confidence intervals
ic_up  = beta_estimates(:, :, 4) - bounds;
ic_low = beta_estimates(:, :, 4) + bounds;


c1 = sum(ic_low <= beta(4), 2);

