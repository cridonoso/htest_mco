function var_hat = estimated_variance(residuals, X, R)
    %ESTIMATED_VARIANCE Summary of this function goes here
    %   Calculates estimated variance
    K  = length(R); % number of parameters
    N = size(X, 1);
    s2 = residuals'*residuals./(N-K); % residual variance
    var_hat = s2*R'*((X'*X)\R); % estimated variance
end

