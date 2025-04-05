clearvars
clc
%% Question 3
rng(42); % set seed to reproduce the same charts    
addpath("utils") % add utils functions

N       = 100; % sample size
n_iter  = 500;
beta    = [1, 2, 1, -1]; % beta coefficients
K       = length(beta);
sigma_e = 1; % Variance values for e
alpha_p22 = 5e-2 ; % significance level for t-test

rejections_homo = zeros(n_iter, 1); 
rejections_white = zeros(n_iter, 1); 

for i = 1:n_iter
    
    X = generate_data(N, [0, 0, 0], [1, 2, 1]); 
    epsilon = sample_normal(N, 0, sqrt(sigma_e)); 
    U = X(:,3) .* epsilon; 
   
    Y = X * beta' + U;

    beta_hat = X \ Y;
    residuals = Y - X * beta_hat;

    sigma2_hat = (residuals' * residuals) / (N - K);
    V_homo = sigma2_hat * inv(X' * X);
    se_homo = sqrt(V_homo(3,3));
    t_homo = (beta_hat(3) - 1) / se_homo;

    S_white = zeros(K);
    for j = 1:N
        x_j = X(j,:)';
        S_white = S_white + (residuals(j)^2) * (x_j * x_j');
    end
    V_white = inv(X' * X) * S_white * inv(X' * X);
    se_white = sqrt(V_white(3,3));
    t_white = (beta_hat(3) - 1) / se_white;

    t_crit = tinv(1 - alpha_p22/2, N - K);

    if abs(t_homo) > t_crit
        rejections_homo(i) = 1;
    end

    if abs(t_white) > t_crit
        rejections_white(i) = 1;
    end
end

fprintf('Rechazos con matriz homocedástica: %d de %d (%.2f%%)\n', ...
    sum(rejections_homo), n_iter, 100*mean(rejections_homo));

fprintf('Rechazos con matriz robusta White: %d de %d (%.2f%%)\n', ...
    sum(rejections_white), n_iter, 100*mean(rejections_white));