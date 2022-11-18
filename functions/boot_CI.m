function CI = boot_CI(X)
% bootstrap CIs based on matrix of values. 
% CI = boot_CI(X)
% rows: cases
% columns: variables 

NBoot = 1000; % number of bootstrap samples 
[M, N] = size(X); 
for n = 1:NBoot
    ind = randi(M, M,1); 
    bootX(n,:) = mean(X(ind, :)); 
end
CI(1,:) = quantile(bootX, .025);
CI(2,:) = quantile(bootX, .975);