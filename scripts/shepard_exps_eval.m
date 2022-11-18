% shepard empirical data analysis 
clear all

% go to top level folder first 
% cd('/Users/kaisiedenburg/Dropbox (Personal)/2dShepard/ExperimentScriptsData/xScriptsData') % please adjust accordingly! 
cd('scripts')
run shepard_load_data.m % gather data
cd(fullfile('..', 'functions'))

%% settings 
set(0,'DefaultAxesFontSize',22)
set(0, 'DefaultAxesLineWidth', 2)
set(0, 'DefaultLineLineWidth', 2);
set(0,'DefaultAxesTitleFontWeight','normal');

%% EXP 1
%% indices 

clear ind

NP = 13; 
for k = 1:NP % participant
    ind.e1.part(:,k) = resp_mat.e1(:, 1) == k; 
end
for k = 1:2 % session
    ind.e1.sess(:,k) = resp_mat.e1(:, 2) == k; 
end

for k = 1:2 % harmonicity
    ind.e1.harm(:,k) = resp_mat.e1(:, 3) == k; 
end

for k = 1:2 % attribute
    ind.e1.attr(:,k) = resp_mat.e1(:, 4) == k; 
end

for k = 1:12 % start
    ind.e1.start(:,k) = resp_mat.e1(:, 5) == k; 
end

for k = 1:11 % shift
    ind.e1.shift(:,k) = resp_mat.e1(:, 6) == k; 
end

corr_respo = zeros(length(ind.e1.shift),1); 
corr_respo(ind.e1.shift(:,6)) = nan; 
corr_respo(sum(ind.e1.shift(:,[1:5]),2)==1) = 1; 
corr_respo(sum(ind.e1.shift(:,[7:11]),2)==1) = 0; 

%% extract data
respo = resp_mat.e1(:,7);
respo(respo<0 | respo > 1) = nan; 
sum(isnan(respo))
restime = resp_mat.e1(:,8);
restime(respo<0 | respo > 1) = nan; 
restime(restime>3) = nan;
restime(restime<.25) = nan;
iscorr = respo == corr_respo; % quantify how well people are doing 
for nP = 1:NP
     for nSh = 1:11 % shift
         for nAt = 1:2 % attr
             for nHa = 1:2 % harm
                ii = ind.e1.part(:,nP) & ind.e1.shift(:,nSh) & ind.e1.attr(:,nAt) & ind.e1.harm(:,nHa);
                res.prop(nSh,nAt,nHa,nP) = nanmean(respo(ii));
                for nSta = 1:12
                    res.propx(nSh,nAt,nHa,nSta,nP) = nanmean(respo(ii & ind.e1.start(:,nSta)));
                end
                res.rt(nSh,nAt,nHa,nP) = nanmedian((restime(ii)));
             end
         end
     end
end
% quantify performance
for nP = 1:NP
    for nAt = 1:2 % attr
         for nHa = 1:2 % harm
             res.corr(nAt, nHa, nP) = mean(iscorr(ind.e1.attr(:,nAt) & ind.e1.harm(:, nHa) & ind.e1.part(:,nP))); 
         end
    end
end
parts = [1 2 4:13]; % participant 3 was at chance throughout!
res.prop_unbi = res.prop - repmat(mean(res.prop),11,1,1,1) + .5;



%% Exp 1 statistical model

X = resp_mat.e1; 

% response variable
response = 1-X(:,7);
response(response < 0) = nan; 
response(response > 1) = nan; 

% prepare data
ds = table(categorical(X(:,1)), categorical(X(:,2)), categorical(X(:,3)), categorical(X(:,4)), ...
    categorical(X(:,5)), (X(:,6)), response);
ds.Properties.VariableNames = {'part', 'session', 'harmonicity', 'attribute', 'start', 'shift', 'resp'};

% DEFINE VARIABLES 
tab = ds; 
DV = {'resp'}; 
IV = {'part', 'harmonicity', 'attribute', 'shift'};
NVars = length(IV);

% get factor levels: 
for nV = 1:NVars
    fac_lev.(IV{nV}) = unique(tab.(IV{nV})); 
    fac_num.(IV{nV}) = 1:length(unique(tab.(IV{nV})));
end

% get array dimensions
for nV = 1:NVars
    dims(nV) = length(fac_lev.(IV{nV}));
end

% get summary statistics for all participants 
[gtab] = grpstats(tab, IV, 'mean', 'DataVars', {'resp'});
emp_dat = permute(reshape(gtab.mean_resp, fliplr(dims)), fliplr([1:NVars])); 

%% EXP 1 GLME model
EXP1_M1 = fitglme(ds, 'resp ~ 1 +  harmonicity * attribute * shift + (1 | part) + (1 | start)', ...
    'Distribution', 'binomial', 'FitMethod', 'REMPL', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
EXP1_M1.Rsquared

% display
[beta, betname, s] = fixedEffects(EXP1_M1);
ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
    num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'}; 
ss

%% PLOT EXP 1
sh = [-.05 .05]; 
figure; 
subplot(1,2,1); hold on
mark_size = 10;
cm = colormap('lines'); 

% EXP 1 GLM prediction
harm_val = [1 -1]; % because using effects coding (where -1 is the last category)
attr_val = [1 -1]; % because numerical predictor is used 

% plot model response (fixed effects)
for nHarm = 1:2 
    for nAttr = 1:2
        subplot(1,2,nHarm)
        clear pred
        pred.harm = harm_val(nHarm);
        pred.attr = attr_val(nAttr); 
        pred.shift = linspace(min(ds.shift), max(ds.shift), 100)'; 
        Xf = make_design_matrix(pred); 
        Xf = [ones(length(Xf),1), Xf, Xf(:,1).*Xf(:,2), Xf(:,1).*Xf(:,3), Xf(:,2).*Xf(:,3), Xf(:,1).*Xf(:,2).*Xf(:,3)]; % add intercept and interactions
        beta_fix = fixedEffects(EXP1_M1);
        yhat = Xf*beta_fix(1:size(Xf,2));
        yhat_link = 1./(1 + exp(-yhat)); % use link function
        plot(pred.shift, yhat_link, '-', 'color', [cm(nAttr,:), .99], 'linewidth', 4); hold on
    end
    ylim([0 1])
    set(gca, 'box', 'off')
end

% add empirical data
co = colormap('lines');
tit = {'FS', 'EN', 'FS inharmonic', 'EN inharmonic'}; 
kk = 0;
sty_lin = {' ', ' '};
sty_mark = {'s', 'o'}; 
add = [-.15 .15];
nAt = 1;
nHa = 1; 
x = squeeze(1-res.prop(:,nAt,nHa, parts))';
for nHa = 1:2 % harmonicity
    kk = kk + 1;
    subplot(1,2,kk)
    for nAt = 1:2 % attribute    
        x = squeeze(1-res.prop(:,nAt,nHa, parts))';
        CI = boot_CI(x);
        mx = mean(x);
        medx = median(x);
        %errorbar([1:11]+add(nAt), mx, mx-CI(1,:), CI(2,:)-mx, sty_lin{nAt}, 'color', [0 0 0], 'linewidth', 1); hold on
        for nEbar = 1:11
            plot([nEbar+add(nAt), nEbar+add(nAt)], [CI(1,nEbar), CI(2,nEbar)], 'color', [0 0 0], 'linewidth', 3)
        end
        plot([1:11]+add(nAt), mx, sty_mark{nAt}, 'color', co(nAt, :), 'linewidth', 3, 'MarkerEdgeColor','k',...
        'MarkerFaceColor', co(nAt,:),'MarkerSize',12); hold on
        xlim([0 12])
        ylim([0 1])
        grid off
        title(tit{kk})
        set(gca,'XTick',[0:12])  % get x axis right 
        set(gca,'XTickLabel',{'' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' ''})
        set(gca, 'box', 'off')
    end
end
subplot(1,2,1)
title('Harmonic')
plot([0 12], [.5 .5], 'k:')
subplot(1,2,2)
plot([0 12], [.5 .5], 'k:')
subplot(1,2,1)
xlabel('Shift [st]')
ylabel('p(down)')
legend('SFS', 'SE')
subplot(1,2,2)
title('Inharmonic')

subplot(1,2,1)
text(-2, 1.1, 'A', 'fontsize', 24, 'fontweight', 'bold')
subplot(1,2,2)
text(-2, 1.1, 'B', 'fontsize', 24, 'fontweight', 'bold')


%% EXP 2

%% get indices

NP = 12;
for k = 1:NP
    ind.e2.part(:,k) = resp_mat.e2(:,1) == k; 
end

for k = 1:2
    ind.e2.harm(:,k) = resp_mat.e2(:,2) == k; 
end

shifts = [1 6 11];
for k = 1:3
    ind.e2.sfs_shift(:,k) = resp_mat.e2(:,4) == shifts(k); 
end

for k = 1:3
    ind.e2.env_shift(:,k) = resp_mat.e2(:,5) == shifts(k); 
end

%% get data
clear resp X
resp.all = resp_mat.e2(:, 8); 
resp.all(resp.all < 0) = nan; 
resp.all(resp.all > 4) = nan; 

for n = 1:NP
    for m = 1:2 % harmonicity 
        for k = 1:3 % FS 
            for l = 1:3  % EN  
                resp.mean(n, k, l, m) = nanmean(resp.all(ind.e2.part(:,n) & ind.e2.harm(:, m) & ind.e2.sfs_shift(:,k) & ind.e2.env_shift(:,l)));
                resp.mean_bin(n, k, l, m) = nanmean(resp.all(ind.e2.part(:,n) & ind.e2.harm(:, m) & ind.e2.sfs_shift(:,k) & ind.e2.env_shift(:,l)) <= 2);
                resp.mean_both(n, k, l, m) = nanmean(resp.all(ind.e2.part(:,n) & ind.e2.harm(:, m) & ind.e2.sfs_shift(:,k) & ind.e2.env_shift(:,l)) > 1 & ...
                    resp.all(ind.e2.part(:,n) & ind.e2.harm(:, m) & ind.e2.sfs_shift(:,k) & ind.e2.env_shift(:,l)) < 4);
                for nC = 1:4
                    mean_count.e2(n, k, l, m, nC) = nansum(resp.all(ind.e2.part(:,n) & ind.e2.harm(:, m) & ind.e2.sfs_shift(:,k) & ind.e2.env_shift(:,l)) == nC);
                end
            end
        end
    end
end

%% Exp 2 - statistical model 
X = resp_mat.e2; 

% response variable
response = X(:,8);
response(response < 0) = nan; 
response(response > 4) = nan; 
response_bin = 1-(response <= 2); 

% prepare data
ds = table(categorical(X(:,1)), categorical(X(:,2)), (X(:,4)), categorical(X(:,5)), ...
    X(:,6), (X(:,7)), response_bin);
ds.Properties.VariableNames = {'part', 'harmonicity', 'FS', 'EN', 'FSstart', 'ENstart', 'resp'};

% DEFINE VARIABLES 
tab = ds; 
DV = {'resp'}; 
IV = {'part', 'harmonicity', 'FS', 'EN'};
NVars = length(IV);

% get factor levels: 
for nV = 1:NVars
    fac_lev.(IV{nV}) = unique(tab.(IV{nV})); 
    fac_num.(IV{nV}) = 1:length(unique(tab.(IV{nV})));
end

% get array dimensions
for nV = 1:NVars
    dims(nV) = length(fac_lev.(IV{nV}));
end


%% EXP 2 GLME model
EXP2_M1 = fitglme(ds, 'resp ~ 1 +  harmonicity * FS * EN + (1 | part) + (1 | FSstart)', ...
    'Distribution', 'binomial', 'FitMethod', 'REMPL', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
EXP2_M1.Rsquared

% display
[beta, betname, s] = fixedEffects(EXP2_M1);
ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
    num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'}; 
ss

% PLOT EXP 2 stat model
sh = [-.05 .05]; 

figure; 
subplot(1,2,1); hold on
mark_size = 10;
cm = colormap('lines'); 
cm = cm(4:end,:); 

% EXP 2 GLM prediction
harm_val = [1 -1]; % because using effects coding (where -1 is the last category)
EN_1_val = [1 0 -1];
EN_6_val = [0 1 -1]; 
%en_val = [1 6 11]; % because numerical predictor is used 

%% plot Exp 2
for nHarm = 1:2 
    for nEN = 1:3
        subplot(1,2,nHarm)
        clear pred
        pred.harm = harm_val(nHarm);
        pred.fs = linspace(min(ds.FS), max(ds.FS), 100)'; 
        %pred.en = en_val(nEN); 
        pred.en1 = EN_1_val(nEN); 
        pred.en2 = EN_6_val(nEN); 
        Xf = make_design_matrix(pred); 
        Xf = [ones(length(Xf),1), Xf, ...
            Xf(:,1).*Xf(:,2), Xf(:,1).*Xf(:,3), Xf(:,1).*Xf(:,4), ...
            Xf(:,2).*Xf(:,3), Xf(:,2).*Xf(:,4), ...
            Xf(:,1).*Xf(:,2).*Xf(:,3), Xf(:,1).*Xf(:,2).*Xf(:,4)]; % add intercept and interactions
        beta_est = fixedEffects(EXP2_M1);
        yhat = Xf*beta_est(1:size(Xf,2)); 
        yhat_link = 1./(1 + exp(-yhat)); % use link function
        plot(pred.fs, yhat_link, '-', 'color', [cm(nEN,:), .99], 'linewidth', 4); hold on
    end
    ylim([0 1])
    set(gca, 'box', 'off')
end

% plot empirical data
sty_lin = {'-', '-', '-'};
sty_mark = {'s', 'o', '^'}; 
xshifts = [-.3 0 .3];
tit = {'Harmonic', 'Inharmonic'};

for nH = 1:2
    subplot(1,2,nH)
    for k = 1:3
        X = 1-squeeze((resp.mean_bin(:,:,k,nH)))
        CI = boot_CI(X);
        for nSH = 1:3
            plot([shifts(nSH)+xshifts(k), shifts(nSH)+xshifts(k)], [CI(1,nSH), CI(2,nSH)] , '-', 'color', [0 0 0], 'linewidth', 2); 
        end
        plot(shifts+xshifts(k), mean(X), sty_mark{k}, 'color', cm(k, :), 'linewidth', 3, 'MarkerEdgeColor','k',...
        'MarkerFaceColor', cm(k,:),'MarkerSize',12); hold on
    end
    set(gca, 'box', 'off')
    title(tit{nH})
    ylim([0 1])
    xlim([0.5 11.5])
    set(gca,'XTick',[1 6 11])  % This automatically sets
    %set(gca,'YTick',[1 2 3 4])  % This automatically sets
    %set(gca,'YTickLabel',{'Up' 'Both (up)' 'Both (down)' 'Down'})  % This automatically sets
    ylabel('p(down)')
    grid off
    plot([0 12], [.5 .5], 'k:')
    
end
subplot(1,2,1)
xlabel('SFS Shift [st]')
legend('SE 1', 'SE 6', 'SE 11')

subplot(1,2,1)
text(-1, 1.1, 'A', 'fontsize', 24, 'fontweight', 'bold')
subplot(1,2,2)
text(-1, 1.1, 'B', 'fontsize', 24, 'fontweight', 'bold')

%% EXP 3

%% indices Exp3A
resp.mean = [];
NP = 12; 
parts = [1:12];
for k = 1:length(parts)
    ind.e3a.part(:,k) = resp_mat.e3a(:,1) == parts(k); 
end

sfs_shifts = [0 1 2 3];
for k = 1:length(sfs_shifts)
    ind.e3a.sfs(:,k) = resp_mat.e3a(:,2) == sfs_shifts(k); 
end

env_shifts = [0 1 6 11];
for k = 1:length(env_shifts)
    ind.e3a.env(:,k) = resp_mat.e3a(:,3) == env_shifts(k); 
end

% get data
resp.all = resp_mat.e3a(:, 6) - 1; 
resp.all(resp.all < 0) = nan; 
resp.all(resp.all > 1) = nan; 

for n = 1:NP
    for m = 1:length(sfs_shifts)
        for k = 1:length(env_shifts)
            resp.mean(n, m, k) = nanmean(resp.all(ind.e3a.part(:,n) & ind.e3a.sfs(:, m) & ind.e3a.env(:, k)));
        end
    end
end


%% Exp 3A - statistical model 
X = resp_mat.e3a; 

% response variable
response = (X(:,6)>1)+0;

% preprocess data
X(X(:,3) == 6,3) = 2; 
X(X(:,3) == 11,3) = 3; 
ind_fs = X(:,2) > 0; 
shift_fac = [X(ind_fs,2); X(~ind_fs,3)];
X = [X, 2-ind_fs, shift_fac]; % attribute factor and shift factor 

% prepare data
ds = table(categorical(X(:,1)), (X(:,5)), categorical(X(:,7)), (X(:,8)), response);
ds.Properties.VariableNames = {'part', 'repetition', 'attr', 'shift', 'resp'};

% DEFINE VARIABLES 
tab = ds; 
DV = {'resp'}; 
IV = {'part', 'attr', 'shift'};
NVars = length(IV);

% get factor levels: 
for nV = 1:NVars
    fac_lev.(IV{nV}) = unique(tab.(IV{nV})); 
    fac_num.(IV{nV}) = 1:length(unique(tab.(IV{nV})));
end

% get array dimensions
for nV = 1:NVars
    dims(nV) = length(fac_lev.(IV{nV}));
end


%% EXP 3A GLME model
EXP3A_M1 = fitglme(ds, 'resp ~ 1 + attr * shift + (1 | part)', ...
    'Distribution', 'binomial', 'FitMethod', 'REMPL', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
EXP3A_M1.Rsquared

% display
[beta, betname, s] = fixedEffects(EXP3A_M1);
ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
    num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'}; 
ss

%% plot Exp 3A
attr_val = [1 -1]; % because using effects coding (where -1 is the last category)
shift_val = [1 2 3]; % because numerical predictor is used 

figure 
subplot(1,2,1); hold on
sh = [-.05 .05]; 
mark_size = 10;
cm = colormap('lines'); 
cm = cm(4:end,:); 

% plot Exp 3
xshift = [-.05, 0, .05];
pshift = [1,2,3]; 
sty_mark = {'s', 'o', '^'}; 
shifts = [1 2 3];
cm = colormap('lines'); 
cm = cm(4:end,:); 
co = colormap('lines'); 

% plot model
for nAtt = 1:2
    clear pred
    pred.attr = attr_val(nAtt); 
    pred.shift = linspace(min(ds.shift), max(ds.shift), 100)'; 
    Xf = make_design_matrix(pred); 
    Xf = [ones(length(Xf),1), Xf, Xf(:,1).*Xf(:,2)];
        % add intercept and interactions
    beta_est = fixedEffects(EXP3A_M1)
    yhat = Xf*beta_est(1:size(Xf,2));
    yhat_link = 1./(1 + exp(-yhat)); % use link function
    plot(pred.shift, yhat_link, '-', 'color', co(nAtt,:), 'linewidth', 4); hold on
end
ylim([0 1])
set(gca, 'box', 'off')


% plot empirical data 
pshift = [1,2,3]; 
xshift = [-.05, .05];
sty_mark = {'s', 'o', '^'}; 

% let's go
X_sfs = squeeze(resp.mean(:,2:end, 1)); 
X_env = squeeze(resp.mean(:,1,2:end)); 
XX = cat(3,X_sfs, X_env); 

CI_sfs = boot_CI(X_sfs);
CI_env = boot_CI(X_env);
CI = cat(3, CI_sfs, CI_env); 

for nAttr = 1:2 
    for nSH = 1:3
        plot([pshift(nSH)+xshift(nAttr), pshift(nSH)+xshift(nAttr)], [squeeze(CI(1,nSH, nAttr)), squeeze(CI(2,nSH, nAttr))], '-', 'color', [0 0 0], 'linewidth', 2); hold on
        plot(pshift+xshift(nAttr), mean(squeeze(XX(:,:,nAttr))), sty_mark{nAttr}, 'color', co(nAttr, :), 'linewidth', 3, 'MarkerEdgeColor','k',...
        'MarkerFaceColor', co(nAttr,:),'MarkerSize',12); hold on
    end
end
xlabel('SFS/SE Shift [st]')
ylabel('p(down)')
grid off
set(gca, 'box', 'off')
plot([0 4], [.5 .5], 'k:')
legend('SFS', 'SE')
set(gca, 'XTick', [1 2 3])
set(gca, 'XTickLabel', {'1/1', '2/6', '3/11'})
xlim([.5 3.5])
%ff = gca; % rotate xticks:
%ff.XTickLabelRotation = 30;


%% indices Exp 3B

NP = length(filens);
for k = 1:length(parts)
    ind.e3b.part(:,k) = resp_mat.e3b(:,1) == parts(k); 
end

env_shifts = [1 6 11];
for k = 1:length(env_shifts)
    ind.e3b.env(:,k) = resp_mat.e3b(:,3) == env_shifts(k); 
end

sfs_shifts = [1 2 3];
for k = 1:length(sfs_shifts)
    ind.e3b.sfs(:,k) = resp_mat.e3b(:,2) == sfs_shifts(k); 
end

% get data
clear resp
resp.all = resp_mat.e3b(:, 6); 
resp.all(resp.all < 0) = nan; 
resp.all(resp.all > 4) = nan;
ind.e3b.valid = ~isnan(resp.all);
mean_count.e3 = [];
for n = 1:NP
    for m = 1:length(env_shifts)
        for l = 1:length(sfs_shifts)
            resp.mean(n, l, m) = nanmean(resp.all(ind.e3b.valid & ind.e3b.part(:,n) & ind.e3b.env(:, m) & ind.e3b.sfs(:, l)));
            resp.mean_bin(n, l, m) = nanmean(resp.all(ind.e3b.valid & ind.e3b.part(:,n) & ind.e3b.env(:, m) & ind.e3b.sfs(:, l))>2);
            resp.mean_both(n, l, m) = nanmean(resp.all(ind.e3b.valid & ind.e3b.part(:,n) & ind.e3b.env(:, m) & ind.e3b.sfs(:, l)) > 1 & ...
                resp.all(ind.e3b.valid & ind.e3b.part(:,n) & ind.e3b.env(:, m) & ind.e3b.sfs(:, l)) < 4);
            for nC = 1:4
                mean_count.e3(n, l, m, nC) = nansum(resp.all(ind.e3b.valid & ind.e3b.part(:,n) & ind.e3b.env(:, m) & ind.e3b.sfs(:, l)) == nC);
            end
        end
    end
end



%% Exp 3B - statistical model 
X = resp_mat.e3b; 

% response variable
response = X(:,6);
response(response < 0) = nan; 
response(response > 4) = nan; 
response_bin = 1-(response <= 2); 

% prepare data
ds = table(categorical(X(:,1)), (X(:,2)), categorical(X(:,3)), (X(:,5)), response_bin);
ds.Properties.VariableNames = {'part', 'FS', 'EN', 'FSstart', 'resp'};

% DEFINE VARIABLES 
tab = ds; 
DV = {'resp'}; 
IV = {'part', 'FS', 'EN'};
NVars = length(IV);

% get factor levels: 
for nV = 1:NVars
    fac_lev.(IV{nV}) = unique(tab.(IV{nV})); 
    fac_num.(IV{nV}) = 1:length(unique(tab.(IV{nV})));
end

% get array dimensions
for nV = 1:NVars
    dims(nV) = length(fac_lev.(IV{nV}));
end


% EXP 3B GLME model
EXP3B_M1 = fitglme(ds, 'resp ~ 1 + FS * EN + (1 | part) + (1 | FSstart)', ...
    'Distribution', 'binomial', 'FitMethod', 'REMPL', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
EXP3B_M1.Rsquared

% display
[beta, betname, s] = fixedEffects(EXP3B_M1);
ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
    num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'}; 
ss

% PLOT EXP 2 stat model
sh = [-.05 .05]; 
subplot(1,2,1); hold on
mark_size = 10;
cm = colormap('lines'); 
cm = cm(4:end,:); 

% EXP 2 GLM prediction
harm_val = [1 -1]; % because using effects coding (where -1 is the last category)
EN_1_val = [1 0 -1];
EN_6_val = [0 1 -1]; 
en_val = [1 6 11]; % because numerical predictor is used 

%% plot Exp 3B
xshift = [-.05, 0, .05];
pshift = [1,2,3]; 
sty_mark = {'s', 'o', '^'}; 
shifts = [1 2 3];
cm = colormap('lines'); 
cm = cm(4:end,:); 

subplot(1,2,2)
for nEN = 1:3
    clear pred
    pred.fs = linspace(min(ds.FS), max(ds.FS), 100)'; 
    pred.en1 = EN_1_val(nEN); 
    pred.en2 = EN_6_val(nEN); 
    Xf = make_design_matrix(pred); 
    Xf = [ones(length(Xf),1), Xf, ...
        Xf(:,1).*Xf(:,2), Xf(:,1).*Xf(:,3)];
        % add intercept and interactions
    yhat = Xf*fixedEffects(EXP3B_M1);
    yhat_link = 1./(1 + exp(-yhat)); % use link function
    plot(pred.fs, yhat_link, '-', 'color', [cm(nEN,:), .99], 'linewidth', 4); hold on
end
ylim([0 1])
set(gca, 'box', 'off')

% plot empirical data
for k = 1:3
    X = squeeze((resp.mean_bin(:,:,k)))
    CI = boot_CI(X);
    for nSH = 1:length(shifts)
        plot([shifts(nSH)+xshift(k),shifts(nSH)+xshift(k)], [CI(1,nSH), CI(2,nSH)], '-', 'color', [0 0 0], 'linewidth', 2); hold on
    end
    %plot(shifts+xshift(k), mean(X), '-', 'color', cm(k,:), 'linewidth', 2); hold on
    plot(shifts+xshift(k), mean(X), sty_mark{k}, 'color', cm(k, :), 'linewidth', 3, 'MarkerEdgeColor','k',...
    'MarkerFaceColor', cm(k,:),'MarkerSize',12); hold on
end
plot([0 4], [.5 .5], 'k:')
set(gca, 'box', 'off')
xlabel('SFS Shift [st]')
legend('SE 1', 'SE 6', 'SE 11')
%title('Log-equidistant')
ylim([0 1])
set(gca,'XTick',[1 2 3])  % This automatically sets
ylabel('p(down)')
grid off 
xlim([.5 3.5])
title('')

subplot(1,2,1)
text(-.1, 1.05, 'A', 'fontsize', 24, 'fontweight', 'bold')
subplot(1,2,2)
text(-.1, 1.05, 'B', 'fontsize', 24, 'fontweight', 'bold')

%% look at histogram [Fig S1, supplementary materials]
figure;

subplot(1,3,1)
X = squeeze(mean(mean_count.e2(:,:,:,1,:)));
XX = squeeze([X(:,1,:); X(:,2,:); X(:,3,:)])
max(sum(XX(:,2:3)/36,2))
mean(sum(XX(:,2:3)/36,2))

mean(sum(XX([1 9],2:3)/36,2))
mean(sum(XX([3 7],2:3)/36,2))

bar(XX/36)
ylabel('Proportion responses')
set(gca, 'Box', 'off')
set(gca, 'XTick', [1:9])
set(gca, 'XTickLabel', {'1-1', '6-1', '11-1', '1-6', '6-6', '11-6', '1-11', '6-11', '11-11'}, 'fontsize', 18)
title('Harmonic')

subplot(1,3,2)
X = squeeze(mean(mean_count.e2(:,:,:,2,:)));
XX = squeeze([X(:,1,:); X(:,2,:); X(:,3,:)])
max(sum(XX(:,2:3)/36,2))
mean(sum(XX(:,2:3)/36,2))

mean(sum(XX([1 9],2:3)/36,2))
mean(sum(XX([3 7],2:3)/36,2))

bar(XX/36)
legend('up', 'both (up dom)', 'both (down dom)', 'down')
set(gca, 'Box', 'off')
set(gca, 'XTick', [1:9])
set(gca, 'XTickLabel', {'1-1', '6-1', '11-1', '1-6', '6-6', '11-6', '1-11', '6-11', '11-11'}, 'fontsize', 18)
title('Inharmonic')

subplot(1,3,3)
X = squeeze(mean(mean_count.e3(:,:,:,:)));
XX = squeeze([X(:,1,:); X(:,2,:); X(:,3,:)])
max(sum(XX(:,2:3)/36,2))
mean(sum(XX(:,2:3)/36,2))
mean(sum(XX(:,[1 4])/36,2))

mean(sum(XX([1 9],2:3)/36,2))
mean(sum(XX([3 7],2:3)/36,2))

bar(XX/36)
set(gca, 'Box', 'off')
set(gca, 'XTick', [1:9])
set(gca, 'XTickLabel', {'1-1', '6-1', '11-1', '1-6', '6-6', '11-6', '1-11', '6-11', '11-11'},'fontsize', 18)
title('Log-equidist.')
