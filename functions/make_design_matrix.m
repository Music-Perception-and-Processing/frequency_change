function [des_mat] = make_design_matrix(des_struc);
% takes structure with fields = factors of the design matrix, makes design
% matrix 

% prelims
nam = fieldnames(des_struc);
L = length(nam);

for n = 1:L % get no of levels for each factor 
    no_fact_levels(n) = length(des_struc.(nam{n})); 
end

% generate full matrix 
matL = prod(no_fact_levels); 

des_mat = ones(matL, L);

for n = 1:L % get lengths of each factor
    reps = matL/prod(no_fact_levels(1:n)); % number of repetitions
    for k = 1:no_fact_levels(n)
        cont = des_struc.(nam{n})(k); % factor content
        xx((k-1)*reps+1:k*reps,1) = repmat(cont, reps, 1);
    end
    des_mat(:, n) = repmat(xx, matL/size(xx,1),1);
    clear xx
end