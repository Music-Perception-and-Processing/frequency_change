function [shift_stats, mod_par] = shepard_model_2d(sig1, sig2, fs, mod_par)
% model spectral shift perception where sig1 and sig2 are the two signals to be compared, 
% fs is the audio sampling frequency, and mod_par are the model parameters

% gammatone fb decomposition 
GT1 = mod_par.gtf(sig1);
GT2 = mod_par.gtf(sig2);

for nSplit = 1:2
    n_split = max(find(mod_par.f_erb < mod_par.f_cutoff(nSplit))); % split parameter
    erbs1 = rms(GT1); % erb spectrum based on gammatone fb object 
    erbs2 = rms(GT2); 

    % do the banana split 
    if nSplit == 1
        erbs1(n_split:end) = 0;
        erbs2(n_split:end) = 0;
    elseif nSplit == 2
        erbs1(1:n_split) = 0;
        erbs2(1:n_split) = 0;
    end

    erbs1_log = 20*log10(erbs1/max(erbs1)); % convert to decibel scale
    erbs2_log = 20*log10(erbs2/max(erbs2));

    erbs1_log(erbs1_log < mod_par.erb_thresh) = mod_par.erb_thresh; % thresholding
    erbs2_log(erbs2_log < mod_par.erb_thresh) = mod_par.erb_thresh; 

    % cross-correlation
    [CC, lags] = xcorr(erbs1_log, erbs2_log, 'normalized'); % CROSS CORRELATE
    CC = CC + mod_par.noise_lev.cc*randn(size(CC));% add noise 
    [CC_max, CC_ind] = sort(CC, 'descend'); % find largest non-zero peak
    lag_best = lags(CC_ind(1)); 
    if lags(CC_ind(1)) == 0 % algo must "decide" whether up down, 0 is not an option
        lag_best = lags(CC_ind(2));
    end
    shift_stats.cc.best_lag = lag_best; % the lag with highest cross-correlation
    ccind = -mod_par.search_cc_lag < lags & lags < mod_par.search_cc_lag; % SEARCH RANGE FOR CC LAGS!!!
    CC = CC'; lags = lags'; 

    if nSplit == 1
        shift_stats.cc_low.centroid = sum(lags(ccind).*CC(ccind))/sum(CC(ccind)); % CC centroid 
    elseif nSplit == 2
        shift_stats.cc_high.centroid = sum(lags(ccind).*CC(ccind))/sum(CC(ccind)); % CC centroid 
    end
    
    % spectral centroid 
    %sc1 = sum(mod_par.f_erb.*(erbs1_log-mod_par.erb_thresh))/sum(erbs1_log-mod_par.erb_thresh); % spectral centroid 
    %sc2 = sum(mod_par.f_erb.*(erbs2_log-mod_par.erb_thresh))/sum(erbs2_log-mod_par.erb_thresh); % spectral centroid 
    %shift_stats.sc.c1 = sc1;
    %shift_stats.sc.c2 = sc2;
    %shift_stats.sc.c_diff = sc1 - sc2; % T1 higher
    shift_stats.sc.c_diff = -100; % T1 higher
end

% autocorrelation computation
for nFB = 1:mod_par.numb_fb % ac in individual filterbanks 
    [xc_mat1(nFB,:), plag1] = xcorr(GT1(:,nFB)); % autocorrelation
    [xc_mat2(nFB,:), plag2] = xcorr(GT2(:,nFB)); % autocorrelation
end

sacf1 = sum(xc_mat1)/max(sum(xc_mat1)); % sum and normalize
sacf2 = sum(xc_mat2)/max(sum(xc_mat2));

auc1 = sacf1 + mod_par.noise_lev.ac*randn(size(sacf1)); % add noise to summed autocorrelation functions
auc2 = sacf2 + mod_par.noise_lev.ac*randn(size(sacf2));
plag1_ms = 1000*plag1/fs; % time lag in ms 
plag2_ms = 1000*plag2/fs; % time lag in ms
search_range1 = 1000/mod_par.ac.fmax  <= plag1_ms & plag1_ms < 1000/mod_par.ac.fmin; % search for pitch between 64 and 128 Hz
search_range2 = 1000/mod_par.ac.fmax <= plag2_ms & plag2_ms < 1000/mod_par.ac.fmin; % search for pitch between 64 and 128 Hz
[~, plag1_ind] = max(auc1.*search_range1); % find aucor max in search range 
[~, plag2_ind] = max(auc2.*search_range2); % find aucor max in search range 
shift_stats.ac.p1 = 1000/plag1_ms(plag1_ind); % estimated pitch
shift_stats.ac.p2 = 1000/plag2_ms(plag2_ind); % estimated pitch
shift_stats.ac.p_diff = shift_stats.ac.p1 - shift_stats.ac.p2;

% get peak of sacf in search range as a feature of pitch salience 
shift_stats.ac.ac_peak1 = max(sum(xc_mat1).*search_range1);
shift_stats.ac.ac_peak2 = max(sum(xc_mat2).*search_range2);

end