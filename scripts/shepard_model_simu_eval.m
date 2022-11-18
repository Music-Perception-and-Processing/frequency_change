% evaluate model simulation 

% go to top level folder 
%cd('/Users/kaisiedenburg/Dropbox (Personal)/2dShepard/ExperimentScriptsData/xScriptsData') % please adjust accordingly! 
cd(fullfile('..', 'data'))

clear all
load simu_resp_mat_v8c.mat
load simu_indmat_v8c.mat
load res_mat_v8c.mat

% empirical measurements from participants 
nHP = 1; 
d.e1 = resp_mat.e1; 
d.e2 = resp_mat.e2; 
d.e3a = resp_mat.e3a; 
d.e3b = resp_mat.e3b;

% simulated data
d.e1s = resp_mat.pc(nHP).e1_simu;
d.e2s = resp_mat.pc(nHP).e2_simu;
d.e3as = resp_mat.pc(nHP).e3a_simu;
d.e3bs = resp_mat.pc(nHP).e3b_simu;


%% generate triplets of full parameter space 
N = 50; 
tripl = []; 
for n1 = 0:N
    for n2 = 0:N-n1
        tripl = [tripl; [n1, n2, N-n1-n2]];
    end
end
tripl = tripl/N; 

% get different weights 
alpha_vec = tripl;

% get projection of triplet into 2d space 
tripl_transl = tripl - [1/3 1/3 1/3]; % center this to unity
[COEFF, SCORE] = pca(tripl_transl); % use pca as transform
figure; plot3(SCORE(:,1), SCORE(:,2), SCORE(:,3), '.')
grid on

alpha_trans = alpha_vec*COEFF; % alpha_vec contains all possible combinations on the grid. 
alpha_min = min(alpha_trans); 
% and translation/scaling: 
alpha_trans = alpha_trans - alpha_min; 
alpha_max = max(alpha_trans); 
alpha_trans = alpha_trans./repmat(alpha_max, length(alpha_trans),1); 


%% do the simulation across all weightings

% SIMULATION [no bootstrapping]
for nL = 1
    nL
    for nAlpha = 1:size(alpha_vec,1)
        nAlpha/size(alpha_vec,1)
        a1 = alpha_vec(nAlpha,1); % flip to have FS left and EN right
        a2 = alpha_vec(nAlpha,2);
        a3 = alpha_vec(nAlpha,3);
        
        % exp 1
        NP = 12; 
        % continuous weighting
        diff.ac = d.e1s(:, 11);
        diff.std.ac = var(diff.ac); 
        diff.cc_low = d.e1s(:, 14);
        diff.std.cc_low = var(diff.cc_low); 
        diff.cc_high = d.e1s(:, 15);
        diff.std.cc_high = var(diff.cc_high); 
        
        for nPart = 1:NP
             for nSh = 1:11 % shift
                 for nAt = 1:2 % attr
                     for nHa = 1:2 % harm
                        ii = ind.e1.part(:,nPart) & ind.e1.shift(:,nSh) & ind.e1.attr(:,nAt) & ind.e1.harm(:,nHa);
                        zz_cont = (a1*diff.ac(ii)/diff.std.ac + a2*diff.cc_low(ii)/diff.std.cc_low + a3*diff.cc_high(ii)/diff.std.cc_high > 0); % continuous weighting, then decision 
                        res.prop_e1.sim(nAlpha, nL).sum_cont(nPart, nSh,nHa,nAt) = nanmean(zz_cont > 0); % implement lapse rate here
                     end
                 end
             end
        end

        % exp 2
        % continuous weighting
        diff.ac = d.e2s(:, 14);
        diff.std.ac = std(diff.ac); 
        diff.cc_low = d.e2s(:, 17);
        diff.std.cc_low = std(diff.cc_low); 
        diff.cc_high = d.e2s(:, 18);
        diff.std.cc_high = std(diff.cc_high); 
        
        %get data
        % clear resp X
        for nPart = 1:NP
            for nHarm = 1:2 % harmonicity
                for nSFS = 1:3 % sfs shift 
                    for nENV = 1:3 % env shift 
                        ii = ind.e2.part(:,nPart) & ind.e2.harm(:, nHarm) & ind.e2.sfs_shift(:,nSFS) & ind.e2.env_shift(:,nENV);
                        zz_cont = (a1*diff.ac(ii)/diff.std.ac + a2*diff.cc_low(ii)/diff.std.cc_low + a3*diff.cc_high(ii)/diff.std.cc_high > 0); % continuous weighting, then decision 
                        res.prop_e2.sim(nAlpha, nL).sum_cont(nPart, nHarm,nSFS,nENV) = nanmean(zz_cont>0); % implement lapse rate here
                    end
                end
            end
        end

        % exp 3 
        % continuous weighting
        diff.ac = d.e3as(:, 11);
        diff.std.ac = std(diff.ac); 
        diff.cc_low = d.e3as(:, 14);
        diff.std.cc_low = std(diff.cc_low); 
        diff.cc_high = d.e3as(:, 15);
        diff.std.cc_high = std(diff.cc_high); 
        
        %Exp 3A
        for nPart = 1:NP 
            for nSFS = 1:4
                for nENV = 1:4
                    ii = ind.e3a.part(:,nPart) & ind.e3a.sfs(:, nSFS) & ind.e3a.env(:, nENV);
                    zz_cont = (a1*diff.ac(ii)/diff.std.ac + a2*diff.cc_low(ii)/diff.std.cc_low + a3*diff.cc_high(ii)/diff.std.cc_high > 0); % continuous weighting, then decision 
                    res.prop_e3a.sim(nAlpha, nL).sum_cont(nPart, nSFS, nENV) = nanmean(zz_cont>0); % implement lapse rate here
                end
            end
        end

        %Exp 3B
        % continuous weighting
        diff.ac = d.e3bs(:, 11);
        diff.std.ac = std(diff.ac); 
        diff.cc_low = d.e3bs(:, 14);
        diff.std.cc_low = std(diff.cc_low); 
        diff.cc_high = d.e3bs(:, 15);
        diff.std.cc_high = std(diff.cc_high); 
        
        for nPart = 1:NP
            for nSFS = 1:3
                for nENV = 1:3
                    ii = ind.e3b.part(:,nPart) & ind.e3b.sfs(:, nSFS) & ind.e3b.env(:, nENV);
                    zz_cont = (a1*diff.ac(ii)/diff.std.ac + a2*diff.cc_low(ii)/diff.std.cc_low + a3*diff.cc_high(ii)/diff.std.cc_high > 0); % continuous weighting, then decision 
                    res.prop_e3b.sim(nAlpha,  nL).sum_cont(nPart, nSFS, nENV) = nanmean(zz_cont>0); % implement lapse rate here
                end
            end
        end
    end
end

%% Bootstrap EMPIRICAL DATA [no alpha]

NBoot = 2000;
% first bootstrap sample is actually full data set!
for nB = 1:NBoot
    nB
    % exp 1
    NP = 12; 

    for nPart = 1:NP
         for nSh = 1:11 % shift
             for nAt = 1:2 % attr
                 for nHa = 1:2 % harm
                    ii = ind.e1.part(:,nPart) & ind.e1.shift(:,nSh) & ind.e1.attr(:,nAt) & ind.e1.harm(:,nHa);
                    if nB == 1 % first bootstrap sample is actually full data set!
                        ii_boot = ii; 
                    else
                        ii_num = find(ii); 
                        ii_boot = randsample(ii_num, length(ii_num), 'true');
                    end
                    empresp = 1-d.e1s(ii_boot,7); 
                    empresp(empresp > 1 | empresp < 0) = nan; 
                    res.prop_e1.emp(nPart, nSh,nHa,nAt,nB) = nanmean(empresp); % empirical data
                 end
             end
         end
    end
    % exp 2
    for nPart = 1:NP
        for nHarm = 1:2 % harmonicity
            for nSFS = 1:3 % sfs shift 
                for nENV = 1:3 % env shift 
                    ii = ind.e2.part(:,nPart) & ind.e2.harm(:, nHarm) & ind.e2.sfs_shift(:,nSFS) & ind.e2.env_shift(:,nENV);
                    if nB == 1
                        ii_boot = ii; 
                    else
                        ii_num = find(ii); 
                        ii_boot = randsample(ii_num, length(ii_num), 'true');
                    end
                    ii_num = find(ii); 
                    ii_boot = randsample(ii_num, length(ii_num), 'true');
                    res.prop_e2.emp(nPart, nHarm, nSFS, nENV,nB) = nanmean(d.e2s(ii_boot,8)>2);
                end
            end
        end
    end

    % exp 3 
    %Exp 3A
    for nPart = 1:NP 
        for nSFS = 1:4
            for nENV = 1:4
                ii = ind.e3a.part(:,nPart) & ind.e3a.sfs(:, nSFS) & ind.e3a.env(:, nENV);
                if nB == 1
                    ii_boot = ii; 
                else
                    ii_num = find(ii); 
                    ii_boot = randsample(ii_num, length(ii_num), 'true');
                end
                ii_num = find(ii); 
                ii_boot = randsample(ii_num, length(ii_num), 'true');
                res.prop_e3a.emp(nPart, nSFS, nENV,nB) = nanmean(d.e3as(ii_boot, 6)-1);
            end
        end
    end
    %Exp 3B
    for nPart = 1:NP
        for nSFS = 1:3
            for nENV = 1:3
                ii = ind.e3b.part(:,nPart) & ind.e3b.sfs(:, nSFS) & ind.e3b.env(:, nENV);
                if nB == 1 % original data for first bootstrap sample! 
                    ii_boot = ii; 
                else
                    ii_num = find(ii); 
                    ii_boot = randsample(ii_num, length(ii_num), 'true');
                end
                ii_num = find(ii); 
                ii_boot = randsample(ii_num, length(ii_num), 'true');
                res.prop_e3b.emp(nPart, nSFS, nENV, nB) = nanmean(d.e3bs(ii_boot, 6) > 2);
            end
        end
    end
end

%% get data in one format !! 

% simulations are averaged across individual trials of simulated subjects,
% i.e. raw data for simulations is the same 
NBoot = 2000;
res.all = []; 
nL = 1; 
nB = 1;
for nB = 1:NBoot
    nB
    for nPart = 1:12
        for nAl = 1:size(res.prop_e1.sim,1)
            for nHarm = 1:2
                for nAtt = 1:2
                    x = squish(mean(res.prop_e1.sim(nAl, nL).sum_cont(:, :, nHarm, nAtt)));
                    res.all.sim(nAl).dat(nHarm,nAtt).x(nPart, :) = x; 
                    y = squeeze((res.prop_e1.emp(nPart, :, nHarm, nAtt, nB)));
                    res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :) = y; 
                end
            end
            % exp 2 harm
            nAtt = 3; 
            for nHarm = 1:2
                x = squish(mean(res.prop_e2.sim(nAl, nL).sum_cont(:, nHarm, :, :)));
                res.all.sim(nAl).dat(nHarm,nAtt).x(nPart, :) = x; 
                y = squish((res.prop_e2.emp(nPart, nHarm, :, :, nB)));
                res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :) = y;  
            end

            % exp 3a1
            nHarm = 3; 
            nAtt = 1;
            x = squish(mean(res.prop_e3a.sim(nAl, nL).sum_cont(:, 2:4, 1)));
            y = squish((res.prop_e3a.emp(nPart, 2:4, 1, nB)));
            res.all.sim(nAl).dat(nHarm,nAtt).x(nPart, :) = x; 
            res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :) = y;  

            nAtt = 2;
            x = squish(mean(res.prop_e3a.sim(nAl, nL).sum_cont(:, 1,2:4)));
            y = squish((res.prop_e3a.emp(nPart, 1, 2:4, nB)));
            res.all.sim(nAl).dat(nHarm,nAtt).x(nPart, :) = x; 
            res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :) = y;  

            % exp 3b 
            nAtt = 3; 
            x = squish(mean(res.prop_e3b.sim(nAl, nL).sum_cont(:, :, :)));
            y = squish((res.prop_e3b.emp(nPart, :, :, nB)));
            res.all.sim(nAl).dat(nHarm,nAtt).x(nPart, :) = x; 
            res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :) = y;  
        end
    end
end

%% compute model performance!!!
% takes around 7h on my macbook
 
disp('compute model performance')
nL = 1; 
for nB = 1:NBoot
    nB
    for nAl = 1:size(alpha_vec,1)
        for nPart = 1:12
            for nHarm = 1:3
                for nAtt = 1:3
                    x = squish(res.all.sim(nAl).dat(nHarm,nAtt).x(nPart,:)); 
                    y = squish(res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :)); 
                    res.eval(nAl, nHarm,nAtt, nB, nPart) = corr(x,y); 
                end
            end
        end
    end
end

%% evaluate these correlations 
% get maxima
for nB = 1:size(res.eval,4)-2; 
    nB
    for nHarm = 1:3
        for nAtt = 1:3
            for nPart = 1:12
                [m, mi] = max((res.eval(:,nHarm,nAtt,nB,nPart)));% get index for best alpha across mean 
                res.max.ind(nPart,nHarm,nAtt, nB) = mi; 
                res.max.val(nPart,nHarm,nAtt, nB) = m; 
            end
        end
    end
end

% get mean maxima across bootstrap
for nHarm = 1:3
    for nAtt = 1:3
        for nPart = 1:12
            res.max.q1(nPart,nHarm,nAtt, :) = quantile(alpha_trans(res.max.ind(nPart,nHarm,nAtt, :),:), .1); % extract  
            res.max.q2(nPart,nHarm,nAtt, :) = quantile(alpha_trans(res.max.ind(nPart,nHarm,nAtt, :),:), .9);
            res.max.m(nPart,nHarm,nAtt, :) = quantile(alpha_trans(res.max.ind(nPart,nHarm,nAtt, :),:), .5); % median 
             
        end
    end
end
% save('res_mat_v8c.mat', 'res',  '-v7.3')

%% evaluate parsimony of model fit

% generate alpha indices 
alpha_vals.mod_ind{1} = find(sum(abs(alpha_vec - [1 0 0]),2) == 0);
alpha_vals.mod_ind{2} = find(sum(abs(alpha_vec - [0 1 0]),2) == 0);
alpha_vals.mod_ind{3} = find(sum(abs(alpha_vec - [0 0 1]),2) == 0);

alpha_vals.mod_ind{4} = find(alpha_vec(:,3) == 0);
alpha_vals.mod_ind{5} = find(alpha_vec(:,2) == 0);
alpha_vals.mod_ind{6} = find(alpha_vec(:,1) == 0);

alpha_vals.mod_ind{7} = 1:length(alpha_vec);  

% step 1: 1d models
for nAl = 1:7
    for nDim = 1:3
        for nHarm = 1:3
            for nPart = 1:12
                res.r2_check(nDim, nHarm, nPart, nAl) = max(res.eval(alpha_vals.mod_ind{nAl}, nHarm,nDim, 1, nPart));
            end
        end
    end
end


%% get mean maxima across bootstrap
for nHarm = 1:3
    for nAtt = 1:3
        for nPart = 1:12
            res.max.qq1(nPart,nHarm,nAtt).x = quantile(alpha_vec(res.max.ind(nPart,nHarm,nAtt, :),:), .05); % extract  
            res.max.qq2(nPart,nHarm,nAtt).x = quantile(alpha_vec(res.max.ind(nPart,nHarm,nAtt, :),:), .95);
            res.max.mm(nPart,nHarm,nAtt).x = quantile(alpha_vec(res.max.ind(nPart,nHarm,nAtt, :),:), .5); % median       
            Q = [res.max.qq1(nPart,nHarm,nAtt).x; res.max.qq2(nPart,nHarm,nAtt).x];
            L1 = [Q([1 2],1), Q([2 1], 2), Q([2 1],3)]; % because data on simplex is dependent, pair smallest with largest value and vice versa 
            L2 = [Q([2 1],1), Q([1 2], 2), Q([2 1],3)];
            L3 = [Q([2 1],1), Q([2 1], 2), Q([1 2],3)];
            
            res.max.trans(nPart, nHarm, nAtt).tl1 = (L1*COEFF(:,1:2) - alpha_min(1:2))./alpha_max(1:2); % lines in 3d transformed into 2d
            res.max.trans(nPart, nHarm, nAtt).tl2 = (L2*COEFF(:,1:2) - alpha_min(1:2))./alpha_max(1:2);
            res.max.trans(nPart, nHarm, nAtt).tl3 = (L3*COEFF(:,1:2) - alpha_min(1:2))./alpha_max(1:2);
            res.max.trans(nPart, nHarm, nAtt).m = (res.max.mm(nPart,nHarm,nAtt).x*COEFF(:,1:2) - alpha_min(1:2))./alpha_max(1:2); % transform in 2d plotting space
        end
    end
end



%% INDEX PERMUTATION TEST 

NPERM = 1000;
res.permcorr = [];
res.perm = [];
nL = 1; 
clear r
nB = 1; 
tic
% compute correlation with the permuted data as basis for fitting 
for nPerm = 1:NPERM
    nPerm
    % instantiate permutation, one for each stimulus condition 
    for k = 1:2
        for kk = 1:2
            r(k,kk).perm =  randperm(11); 
        end
    end
    for k = 1:3
        r(k,3).perm = randperm(9);
    end
    for k = 1:2
        r(3,k).perm = randperm(3);
    end
   
    for nHarm = 1:3
        for nAtt = 1:3
            clear x0; 
            for nPart = 1:12 
                for nAlpha = 1:size(alpha_vec,1)
                    x = squish(res.all.sim(nAlpha).dat(nHarm,nAtt).x(nPart,:)); 
                    y = squish(res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :));
                    % correlation for model with permuted shift indices 
                    res.perm(nAlpha, nHarm,nAtt, nL, nPart) = corr(x(r(nHarm,nAtt).perm),y); 
                end
                % get maxima for permuted model 
                [m, mi] = max((res.eval(:,nHarm,nAtt,nB,nPart)));% get max across alphas 
                res.max.indp(nPart,nHarm,nAtt,nB) = mi; 
                res.max.valp(nPart,nHarm,nAtt,nB) = m; 
                x0(nHarm,nAtt).x(nPart,:) = squish(res.all.sim(mi).dat(nHarm,nAtt).x(nPart,:))';% assemble best individual models for permuted data
            end
            xx = mean(x0(nHarm,nAtt).x,1)'; % model average (across indiv fitted versions)
            yy = squish(mean(res.all.emp(nB).dat(nHarm,nAtt).y(:, :))); % mean response pattern
            res.permcorr(nHarm,nAtt,nPerm) = corr(xx(r(nHarm, nAtt).perm),yy); % correlation of models optimized to permuted index
        end
    end
end
toc

% extract comparison matrix, .95 percentiles of results with 
% randomly shuffled indices. If method makes it easy to overfit, then there
% should be a lot of instances where performance is above model performance
% (e.g. > .95) 

%% fit model on restricted SUBSETS OF SHIFTS 
% leave experiments as they are, just leave out conditions for trainings

nL = 1; 
nB = 1; 

% generate index sets : n choose k 

n = 9; 
for k = 1:n
    ind.n9(k).nck = nchoosek(1:n,k);
end
n = 11; 
for k = 1:n
    ind.n11(k).nck = nchoosek(1:n,k);
end
n = 3; 
for k = 1:n
    ind.n3(k).nck = nchoosek(1:n,k);
end
    
%% evaluate correlations for every possible index set     

for nAlpha = 1:size(alpha_vec,1)
    nAlpha/size(alpha_vec,1)
    for nPart = 1:12
        for nHarm = 1:3 
            for nAtt = 1:3 
                x = res.all.sim(nAlpha).dat(nHarm,nAtt).x(nPart, :)';
                y = res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :)';
                Lx = length(x);
                for k = 1:Lx % different lengths of the data depending on condition
                    if Lx == 11
                        for nInd = 1:size(ind.n11(k).nck,1) % check all different combinations of index sets given k
                            iii = ind.n11(k).nck(nInd,:);
                            res.pred(nAlpha, nHarm,nAtt, k, nInd, nPart) = corr(x(iii),y(iii));
                        end
                    elseif Lx == 9
                        for nInd = 1:size(ind.n9(k).nck,1) % check all different combinations of index sets given k
                            iii = ind.n9(k).nck(nInd,:);
                            res.pred(nAlpha, nHarm,nAtt, k, nInd, nPart) = corr(x(iii),y(iii));
                        end
                    elseif Lx == 3
                        for nInd = 1:size(ind.n3(k).nck,1) % check all different combinations of index sets given k
                            iii = ind.n3(k).nck(nInd,:);
                            res.pred(nAlpha, nHarm,nAtt, k, nInd, nPart) = corr(x(iii),y(iii));
                        end
                    end
                end
            end
        end
    end
end

%% get maxima
for k = 1:11 %1:size(res.eval,4)-2; 
     for nInd = 1:size(ind.n11(k).nck,1) % check all different combinations of index sets given k
        for nHarm = 1:3
            for nAtt = 1:3
                for nPart = 1:12
                    [m, mi] = max((res.pred(:,nHarm,nAtt,k, nInd, nPart)));% get index for best alpha across mean 
                    res.max_pred.ind(nPart,nHarm,nAtt,k,nInd) = mi; 
                    res.max_pred.val(nPart,nHarm,nAtt,k,nInd) = m; 
                end
            end
        end
     end
end

% save('res_mat_v8c.mat', 'res',  '-v7.3')

