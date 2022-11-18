% evaluate model simulation 

clear all

% go to top level folder 
% cd('/Users/kaisiedenburg/Dropbox (Personal)/2dShepard/ExperimentScriptsData/xScriptsData') % please adjust accordingly! 
cd(fullfile('..', 'data'))

load simu_resp_mat_v8c.mat
load simu_indmat_v8c.mat
load res_mat_v8c.mat

cd(fullfile('..', 'functions'))

%% plot settings 
set(0,'DefaultAxesFontSize',22)
set(0, 'DefaultAxesLineWidth', 2)
set(0, 'DefaultLineLineWidth', 2);
set(0,'DefaultAxesTitleFontWeight','normal');


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


%% specific example of optimization for model schematic plot. 
% [part of preprint Fig 5]

% plot triplets
subplot(1,2,1)
plot3(tripl(:,1), tripl(:,2), tripl(:,3), '.', 'color', [.5 .5 .5])
grid on
xlabel('AC')
ylabel('CCres')
zlabel('CCunres')
text(-.3,1.1, 1.5,'A', 'fontsize', 20, 'FontWeight', 'Bold')


% plot maximization surface
subplot(1,2,2) 
nHarm = 2;
nAtt = 3;
nPart = 1;
nB = 1; 

dt = delaunayTriangulation(alpha_trans(:,1),alpha_trans(:,2)) ;
tri = dt.ConnectivityList ;
%trisurf(tri,alpha_trans(:,1),alpha_trans(:,2),median(res.eval(:,nHarm,nAtt,nLapse,:),5))
trisurf(tri,alpha_trans(:,1),alpha_trans(:,2),(res.eval(:,nHarm,nAtt,nB,nPart)).^2, 'LineStyle','none'); hold on
at = alpha_trans(res.max.ind(nPart,nHarm,nAtt, nB),:);
plot3([at(1), at(1)], [at(2), at(2)], [0 res.max.val(nPart,nHarm,nAtt, nB).^2], ':k')


alph_grid = ... % x1 x2 y1 y2
    [0 1 0 0; 
    1/6 5/6 1/3 1/3; 
    1/3 2/3 2/3 2/3;
    0 1 0 0; 
    1/3 2/3 0 2/3;
    2/3 5/6 0 1/3;
    2/3 1/3 0 2/3;
    1/3 1/6 0 1/3;
    0 1/2 0 1; 
    1/2 1 1 0; ...
    ];
for n = 1:size(alph_grid,1)
    plot3(alph_grid(n, [1 2]), alph_grid(n, [3 4]), [0 0], 'color', [.5 .5 .5], 'linewidth', 2); 
end
text(1/2-.05, 1 + .1, 'AC', 'fontsize', 16)
text(1-.1, -.1, 'CCres',  'fontsize', 16)
text(0-.05, -.1, 'CCunres',  'fontsize', 16)
grid on
axis on
zlabel('R^2')
set(gca, 'XTick', [.5])
%text(-.3,1.1, 1.5,'B', 'fontsize', 20, 'FontWeight', 'Bold')
colormap('winter')
set(gca, 'XTick', [])
set(gca, 'YTick', [])


%% plot model predictions for all conditions 
% [preprint Fig 7]
figure; 
cm = colormap('lines'); 

for n1 = 1:9
    subplot(3,3,n1)
    ylim([0 1])
    set(gca, 'Box', 'off')
end


com = colormap('lines');
com(3,:) = com(4,:); 
cm(1,:) = [0 0 0];% AC 
cm(2,:) = [.5 .5 .5]; % CC 
cm(3,:) = [.5 .5 .5]; % model average
linsty = {'--', ':', '-.', '--'}; 
clear y

% EMPIRICAL DATA 
kk = 0;
for nHarm = 1:3
    for nAtt = 1:3
        kk = kk + 1;
        subplot(3,3,kk); hold on
        yy = squish(mean(res.all.emp(nB).dat(nHarm,nAtt).y(:, :))); % empirical data 
        if nAtt == 3
            plot([1:3], yy(1:3), 'color', [com(nAtt,:), .5], 'linewidth', 8); hold on
            plot([4:6], yy(4:6), 'color', [com(nAtt,:), .5], 'linewidth', 8); hold on
            plot([7:9], yy(7:9), 'color', [com(nAtt,:), .5], 'linewidth', 8); hold on
        else
            plot(yy, 'color', [com(nAtt,:), .5], 'linewidth', 8); hold on
        end
        
        set(gca, 'XTick', [1:11])
        set(gca, 'XTickLabel', {'1', '', '', '', '', '6', '', '', '', '', '11'})
    end
end


% SIMULATION (averaged models)
kk = 0;
for nHarm = 1:3
    for nAtt = 1:3
        kk = kk + 1;
        subplot(3,3,kk); hold on
        yy = []; 
        for nPart = 1:12
            nAlpha = res.max.ind(nPart,nHarm,nAtt,1); % best index 
            yy = [yy; squish((res.all.sim(nAlpha).dat(nHarm,nAtt).x(nPart, :)))'];% simulation 
        end
        yym = mean(yy); 
        if nAtt == 3
            plot([1:3], yym(1:3), 'color', [com(nAtt,:)], 'linewidth', 4); hold on
            plot([4:6], yym(4:6), 'color', [com(nAtt,:)], 'linewidth', 4); hold on
            plot([7:9], yym(7:9), 'color', [com(nAtt,:)], 'linewidth', 4); hold on
        else
            plot(yym, 'color', [com(nAtt,:)], 'linewidth', 4); hold on
        end
        yy = squish(mean(res.all.emp(nB).dat(nHarm,nAtt).y(:, :))); % empirical data 
        rv = corrplot4(yym,yy, com(nAtt,:)) % compare it ! 
        cval(nHarm, nAtt) = rv; 
        set(gca, 'XTick', [1:11])
        set(gca, 'XTickLabel', {'1', '', '', '', '', '6', '', '', '', '', '11'})
    end
end


% RAW FEATURE DATA
% simulated data: AC and CC
nPlot = 0; 
feat = {'AC', 'CC_low', 'CC_high'}; 
nFeat = 0; 

% get special indices: 
iii = [find(sum(abs(alpha_vec-[1 0 0]),2) == 0); 
    find(sum(abs(alpha_vec-[0 1 0]),2) == 0); 
    find(sum(abs(alpha_vec-[0 0 1]),2) == 0)];
    
kk = 0;
for nHarm = 1:3
    for nAtt = 1:3
        kk = kk + 1;
        subplot(3,3,kk); hold on
        for nFeat = 1:3 % plot raw features 
            nPlot = nFeat; 

            yy = []; 
            for nPart = 1:12
                %nAlpha = res.max.ind(nPart,nHarm,nAtt,1); % best index 
                nAlpha = iii(nFeat); 
                yy = [yy; squish((res.all.sim(nAlpha).dat(nHarm,nAtt).x(nPart, :)))'];% simulation 
            end
            yym = mean(yy); 
            if nAtt == 3
                plot([1:3], yym(1:3), linsty{nPlot}, 'color', cm(nPlot, :), 'linewidth', 2); hold on
                plot([4:6], yym(4:6), linsty{nPlot}, 'color', cm(nPlot, :), 'linewidth', 2); hold on
                plot([7:9], yym(7:9), linsty{nPlot}, 'color', cm(nPlot, :), 'linewidth', 2); hold on
            else
                plot(yym, linsty{nPlot}, 'color', cm(nPlot, :), 'linewidth', 2); hold on
            end
            set(gca, 'XTick', [1:11])
            set(gca, 'XTickLabel', {'1', '', '', '', '', '6', '', '', '', '', '11'})
            xlim([0 12])
            ylim([0 1])
        end
    end
end

subplot(3,3,1)
title('SFS')
legend('emp', 'model', 'AC', 'CCres', 'CCunr')
subplot(3,3,2)
title('SE')
subplot(3,3,1)
ylabel('Harmonic')
subplot(3,3,4)
ylabel('Inharmonic')
subplot(3,3,7)
ylabel('Log-Equidist')
xlabel('Shift [st]')

subplot(3,3,3)
set(gca, 'XTick', [1:9])
set(gca, 'XTickLabel', {'1-1', '6-1', '11-1', '1-6', '6-6', '11-6', '1-11', '6-11', '11-11'}, 'fontsize', 12)
xtickangle(90)
set(gca, 'YTick', [0, .5, 1], 'fontsize', 22)

subplot(3,3,6)
set(gca, 'XTick', [1:9])
set(gca, 'XTickLabel', {'1-1', '6-1', '11-1', '1-6', '6-6', '11-6', '1-11', '6-11', '11-11'}, 'fontsize', 12)
xtickangle(90)
set(gca, 'YTick', [0, .5, 1], 'fontsize', 22)

subplot(3,3,9)
set(gca, 'XTick', [1:9])
set(gca, 'XTickLabel', {'1-1', '2-1', '3-1', '1-6', '2-6', '3-6', '1-11', '2-11', '3-11'}, 'fontsize', 12)
xtickangle(90)
set(gca, 'YTick', [0, .5, 1], 'fontsize', 22)

subplot(3,3,3)
title('SFS-SE', 'fontsize', 22)

subplot(3,3,7)
xlim([0 4])
set(gca, 'XTick', [1:3])
set(gca, 'XTickLabel', {'1', '2', '3'})

subplot(3,3,8)
xlim([0 4])
set(gca, 'XTick', [1:3])
set(gca, 'XTickLabel', {'1', '6', '11'})

subplot(3,3,3)
xlim([0 10])

subplot(3,3,6)
xlim([0 10])

subplot(3,3,9)
xlim([0 10])

%% quantitative correlations 
conf_lev = 0.05/9; 
for nFeat = 1:3
    for nHarm = 1:3
        for nAtt = 1:3
            xx = []; 
            for nPart = 1:12
                %nAlpha = res.max.ind(nPart,nHarm,nAtt,1); % best index 
                nAlpha = iii(nFeat); 
                xx = [xx; squish((res.all.sim(nAlpha).dat(nHarm,nAtt).x(nPart, :)))'];% simulation 
            end
            xxm = mean(xx); 
            yy = squish(mean(res.all.emp(nB).dat(nHarm,nAtt).y(:, :))); % empirical data 
            [r,p] = corr(xxm',yy); 
            cormat(nFeat).r(nHarm,nAtt) = r; 
            cormat(nFeat).p(nHarm,nAtt) = p; 
            
        end
    end
end

disp(round([cormat(1).r; cormat(2).r; cormat(3).r].^2,2)) 

disp([cormat(1).p; cormat(2).p; cormat(3).p] < .05/9) 

%% how correlated are the predictors? 

% get data
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

% look into correlations 
cmat = [d.e1s(:, 11),d.e1s(:, 14),d.e1s(:, 15)];
[cr, pmat] = corr(cmat);
cr.*(pmat < .01)
% pretty uncorrelated for exp1

cmat = [d.e2s(:, 14),d.e2s(:, 17),d.e2s(:, 18)];
[cr, pmat] = corr(cmat);
cr.*(pmat < .01)
% and also pretty uncorrelated for exp2

[cmat] = [d.e3bs(:, 11),d.e3bs(:, 14),d.e3bs(:, 15)];
[cr, pmat] = corr(cmat);
cr.*(pmat < .01)
% and also pretty uncorrelated for exp2


% and all in one: 

cmat = [[d.e1s(:, 11);d.e2s(:, 14);d.e3bs(:, 11)], ...
    [d.e1s(:, 14); d.e2s(:, 17);  d.e3bs(:, 14)], ...
    [d.e1s(:, 15); d.e2s(:, 18); d.e3bs(:, 15)]];

[cr, pmat] = corrcoef(cmat);
cr
pmat
cr.*(pmat < .01)

%% INDEX PERMUTATION TEST 
% plot bootstrap distribution
NPERM = 1000; 
figure; 
pps = [1:NPERM];
k = 0;
for nHarm = 1:3
    for nAtt = 1:3
        k = k + 1;
        subplot(3,3,k); hold on
        res.boot_q99(nHarm, nAtt) = quantile(squeeze(res.permcorr(nHarm, nAtt,pps)), .99);
        res.boot_max(nHarm, nAtt) = max(squeeze(res.permcorr(nHarm, nAtt,pps))); 
        res.boot_med(nHarm, nAtt) = median(squeeze(res.permcorr(nHarm, nAtt,pps))); 
        res.boot_p(nHarm, nAtt) = mean(squeeze(res.permcorr(nHarm, nAtt,pps) >= cval(nHarm,nAtt))); 
        histogram(res.permcorr(nHarm, nAtt,pps), 'Normalization', 'probability')
        plot([cval(nHarm, nAtt) cval(nHarm, nAtt)], [0 1], 'r--')
        plot([res.boot_q99(nHarm, nAtt) res.boot_q99(nHarm, nAtt)], [0 1], 'b:')
        ylim([0 .2])
        xlim([-1.1 1.1])
    end
end

subplot(3,3,1)
title('SFS')

subplot(3,3,2)
title('SE')
subplot(3,3,3)
title('SFS-SE')

subplot(3,3,1)
ylabel('Harmonic')

subplot(3,3,4)
ylabel('Inharmonic')

subplot(3,3,7)
ylabel('Log-equidist')
lg = legend('Random indices', 'Actual model', '99th percentile')
lg.FontSize = 12;

subplot(3,3,8)
xlabel('Pearson Correlation')


res.boot_p
res.boot_max
cval
% i.e., our results clearly outperform anything that would have arisen by
% chance [apart from exp3a]


%% plot distrobution of optimal alphas
% [preprint fig 8]

% colormapping 
co = colormap('lines'); 
coma(1,:) = co(1,:);
coma(2,:) = co(2,:);
coma(3,:) = co(4,:);
coma(4:6,:) = co([1 2 4],:);
coma(7:9,:) = co([1 2 4],:);

% settings
plo.lw = 2;
plo.symb = {'s', 's', 's'; 'o', 'o', 'o'; 'd', 'd', 'd'}; 
plo.jit = .025*rand(length(alpha_trans),2)-0.0125;
transpar = .1; 
edgecol = [.1 .1 .1];
%edgecol = [1 1 1];
edgewid = 2; 


% extract maxima for density map 
empdens = zeros(size(res.eval)); 
r2dens = zeros(size(res.eval)); 
for nB = 1:2000 % this should go much better  
    for nP = 1:12
        for nHarm = 1:3
            for nAtt = 1:3
                [~, maxInd] = max(res.eval(:,nHarm, nAtt, nB, nP));
                empdens(maxInd, nHarm,nAtt, nB, nP) = 1; % get empirical density distribution 
                r2dens(:,nHarm,nAtt,nB,nP) = rescale(res.eval(:,nHarm, nAtt, nB, nP).^2);
            end
        end
    end
end
densRes = log10(mean((mean(empdens, 4)),5) + 10^-6); % taking the mean of log(p) [empirical density]
r2Res = mean(mean(r2dens, 4),5);

% plot it
dt = delaunayTriangulation(alpha_trans(:,1),alpha_trans(:,2)) ;
tri = dt.ConnectivityList ;

cmnew = flipud(colormap('gray')); 
%for nP = 1:12
    figure;
    k = 0;
for nHarm = 1:3 % all in one pic
    for nAtt = 1:3
        k = k + 1;
        subplot(3,3,k); hold on
%         for nPart = 1:12  % OLD WAY OF PLOTTING DENSITY, UNSCALED 
%             alpha_opt = alpha_trans(squeeze(res.max.ind(nPart,nHarm,nAtt,:)), :); % indices of all bootstrap samples per participant, fed into alpha 
%             [f,xi] = ksdensity(alpha_opt(:,1:2))
%             plot1 = scatter(alpha_opt(:,1), alpha_opt(:,2), plo.symb{nHarm,nAtt}, ...
%                 'Linewidth', 2, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [coma(k, :)], 'MarkerFaceAlpha',transpar, 'MarkerEdgeAlpha', .1);
%         end
        hold on 
        trisurf(tri,alpha_trans(:,1),alpha_trans(:,2),squeeze(densRes(:,nHarm, nAtt)), 'LineStyle','none'); hold on
        %trisurf(tri,alpha_trans(:,1),alpha_trans(:,2), r2Res(:, nHarm, nAtt), 'LineStyle','none'); hold on
        view([0 90])
        caxis([-6 0])
        set(gca, 'colormap', cmnew)
        %colorbar
        for nPart = 1:12
           ci_dat0 = squeeze(res.max.m(nPart,nHarm,nAtt,:)); % % mean value across bootstrap samples 
           %ci_dat0 = alpha_trans(res.max.ind(nPart,nHarm,nAtt,1),:); % on full data set.
           plot3(ci_dat0(1), ci_dat0(2), 1, plo.symb{nHarm,nAtt}, 'Linewidth', edgewid , ...
               'MarkerEdgeColor',edgecol, 'MarkerFaceColor', coma(k, :), 'MarkerSize', 10)
        end
        triangle_plot3();
    end
end

subplot(3,3,1)
title('SFS')
subplot(3,3,2)
title('SE')
subplot(3,3,3)
title('SFS-SE')

subplot(3,3,1)
t = text(-.2, .1, 'Harmonic', 'fontsize', 22);
set(t, 'Rotation', 90)
subplot(3,3,4)
t = text(-.2, 0, 'Inharmonic', 'fontsize', 22);
set(t, 'Rotation', 90)

subplot(3,3,7)
t = text(-.2, 0, 'Log-equidist', 'fontsize', 22);
set(t, 'Rotation', 90)

subplot(3,3,8)
colorbar
colorbar('XTick', [-6:2:0], 'XTickLabel', {'-6' '10^-4' '10^-2' '10^0'})

%% supplementary figure on necessity of full model. 

% colormapping 
co = colormap('lines'); 
coma(1,:) = co(1,:);
coma(2,:) = co(2,:);
coma(3,:) = co(4,:);

% settings
plo.lw = 0;
plo.symb = {'s', 's', 's'; 'o', 'o', 'o'; 'd', 'd', 'd'}; 
%plo.jit = .015*rand(length(alpha_trans),2)-0.0125;
%transpar = .1; 
edgecol = .7*[1 1 1];
edgewid = 0; 

xpos = 0.2*rand(12,1); 

mod_variants = {'AC', 'CCres', 'CCunr', 'AC+CCres', 'AC+CCunr', 'CCres+CCunr', 'AC+CCres+CCunr'}; % model variants 
colexp = [.5 .5 .5, 1 1 1, 4];
figure;
k = 0; 

for nHarm = 1:3
    for nDim = 1:3
        k = k + 1;
        subplot(3,3,k); hold on  
        for nAl = 1:6
            Y = squeeze(res.r2_check(nDim, nHarm, :,:));
            Y(Y<0) = 0; 
            Y = Y.^2; 
            XX = Y(:,7) - Y(:,nAl);
            %XX = Y(:,nAl); 
            plot(nAl - xpos - .1, XX, '.', 'color', coma(nDim, :).^colexp(nAl), 'linewidth', 2); hold on
            if k == 1 | k == 5 | k == 9
                t = text(nAl + .15, 0.3, mod_variants{nAl}, 'fontsize', 14);
                set(t, 'Rotation', 90, 'color', coma(nDim,:).^colexp(nAl))
            end
            ci = boot_CI(XX);
            plot([nAl + .15, nAl + .15], [ci(1), ci(2)], 'linewidth', 2, 'color', [.5 .5 .5]);
            plot(nAl + .15, nanmean(XX), 's', 'MarkerEdgeColor', coma(nDim, :).^colexp(nAl), 'MarkerFaceColor', coma(nDim, :).^colexp(nAl), 'MarkerSize', 6); hold on
        end
        xlim([0.5 6.5])
        ylim([0 1])
        set(gca, 'XTick', [1:7])
        set(gca, 'XTickLabel', '')
    end
end


subplot(3,3,1)
title('SFS')
subplot(3,3,2)
title('SE')
subplot(3,3,3)
title('SFS-SE')

subplot(3,3,1)
ylabel('Harmonic')

subplot(3,3,4)
ylabel('Inharmonic')

subplot(3,3,7)
ylabel('Log-equidist')

%% plot optimal alphas for four subjects who completed all conditions
fab4ind = ...
[1 5 11 9;
 1 3 4 5; 
 1 3 5 6]

p_map = cat(3, [1 1 1; 1 1 1; 1 1 1], [5 5 3; 5 5 3; 3 3 3], ...
    [11 11 4; 11 11 4; 5 5 5], [9 9 5; 9 9 5; 6 6 6]);

figure
for nPP = 1:4
    k = 0;
    for nHarm = 1:3 % all in one pic
        for nAtt = 1:3
            k = k + 1;
            subplot(3,3,k); hold on
            hold on 
            view([0 90])
            caxis([-6 0])
            set(gca, 'colormap', cmnew)
               ci_dat0 = squeeze(res.max.m(nPP,nHarm,nAtt,:)); % % mean value across bootstrap samples 
               text(ci_dat0(1), ci_dat0(2),  strcat('P', num2str(nPP)), 'fontsize', 12, 'fontweight', 'bold', 'color', coma(k, :))
            triangle_plot3();
        end
    end
end

subplot(3,3,1)
title('SFS')
subplot(3,3,2)
title('SE')
subplot(3,3,3)
title('SFS-SE')

subplot(3,3,1)
t = text(-.2, .1, 'Harmonic', 'fontsize', 22);
set(t, 'Rotation', 90)
subplot(3,3,4)
t = text(-.2, 0, 'Inharmonic', 'fontsize', 22);
set(t, 'Rotation', 90)

subplot(3,3,7)
t = text(-.2, 0, 'Log-equid', 'fontsize', 22);
set(t, 'Rotation', 90)

subplot(3,3,8)
%colorbar 
%colorbar('XTick', [-6:2:0], 'XTickLabel', {'-6' '10^-4' '10^-2' '10^0'})


%% evaluate predictive power ... 


%% GENERALIZATION ACROSS PARTICIPANTS
nB = 1; 
for nPart = 1:12
    for nPartTest = 1:12
        for nHarm = 1:3
            for nAtt = 1:3
                alpha_ind = res.max.ind(nPart,nHarm,nAtt,1);  % best alpha per participant 
                x = squish(res.all.sim(alpha_ind).dat(nHarm,nAtt).x(1,:));
                y = squish(res.all.emp(nB).dat(nHarm,nAtt).y(nPartTest, :));
                res.pred_part(nPart, nPartTest, nHarm,nAtt) = corr(x,y); % predict across participants 
            end
        end
    end
end

%% plotti karotti across parts
figure; 
colormap('parula')
nPl = 0; 
res.pred_part_mean = [];
for nHarm = 1:3
    for nAtt = 1:3
        nPl = nPl + 1;
        subplot(3,3,nPl)
        image(256*(res.pred_part(:, :, nHarm,nAtt)).^2)
        res.pred_part_mean(nHarm, nAtt) = mean(squish(res.pred_part(:, :, nHarm,nAtt).^2)); 
        set(gca, 'XTick', [2:2:12])
        set(gca, 'XTickLabel', {'2', '', '6', '', '10', ''})
        set(gca, 'YTick', [2:2:12])
        set(gca, 'YTickLabel', {'2', '', '6', '', '10', ''})
    end
end
subplot(3,3,9)
colorbar
colorbar('XTick', linspace(1,256, 5), 'XTickLabel', {'0', '.25', '.5', '.75' ,'1'})

subplot(3,3,1)
t = text(-6, -1, 'B', 'fontsize', 20, 'fontweight', 'bold');
title('SFS')
subplot(3,3,2)
title('SE')
subplot(3,3,3)
title('SFS-SE')

subplot(3,3,1)
ylabel('Harm')
subplot(3,3,4)
ylabel('Inharm')
subplot(3,3,7)
ylabel('Log-eq')

subplot(3,3,4)
t = text(-6, 11, 'Fitted', 'fontsize', 22);
set(t, 'Rotation', 90)
subplot(3,3,8)
xlabel('Tested')

% numerically:
res.pred_part_mean


%% GENERALIZATION ACROSS CONDITIONS OF 4 PARTS             

fab4ind = [1 5 11 9;
            1 3 4 5; 
            1 3 5 6]

expmap = [1 1 2; 1 1 2; 3 3 3]; 

nB = 1; 

lin_conds = [ 1 1 1 2 2 2 3 3 3; % harmonicity 
         1 2 3 1 2 3 1 2 3]% attributes 

res.pred_part_acr_con = [];

for nCond = 1:9
    nHarm = lin_conds(1, nCond); 
    nAtt = lin_conds(2, nCond); 
    nExp = expmap(nHarm, nAtt); 
    parts = fab4ind(nExp,:); 

    for nCondTest = 1:9
        nHarmTest = lin_conds(1, nCondTest); 
        nAttTest = lin_conds(2, nCondTest); 
        nExpTest = expmap(nHarmTest, nAttTest); 
        partsTest = fab4ind(nExpTest,:); 

        for nPart = 1:4
        alpha_ind = res.max.ind(parts(nPart),nHarm,nAtt,1);  % best alpha per participant for fitting condition
        x = squish(mean(res.all.sim(alpha_ind).dat(nHarmTest,nAttTest).x(:, :))); % use that weight for model simulation in test condition
        y = squish(res.all.emp(nB).dat(nHarmTest,nAttTest).y(partsTest(nPart), :)); % get empirical data of subject in test condition
        res.pred_part_acr_con(nCond, nCondTest, nPart) = corr(x,y); % compare model prediction and data
        end
    end     
end

% plot it 
figure 
for nPart = 1:4
    subplot(2,2,nPart)
    mm = squeeze(res.pred_part_acr_con(:, :, nPart))
    mm(mm < 0) = nan;
    image(256*(mm.^2))
    title(strcat('P', num2str(nPart)))
    set(gca, 'YTick', [1:9])
    set(gca, 'YTickLabel', {'ha SFS','ha SE', 'ha SFS-SE', 'ih SFS', 'ih SE', 'ih SFS-SE', 'eq SFS', 'eq SE', 'eq SFS-SE'}, 'fontsize',16)
    set(gca, 'XTick', [1:9])
    set(gca, 'XTickLabel', {'ha SFS','ha SE', 'ha SFS-SE', 'ih SFS', 'ih SE', 'ih SFS-SE', 'eq SFS', 'eq SE', 'eq SFS-SE'}, 'fontsize',16)
    xtickangle(45)
end
colorbar
colorbar('XTick', linspace(1,256, 5), 'XTickLabel', {'0', '.25', '.5', '.75' ,'1'})
subplot(2,2,1)
ylabel('Fitted Condition')

subplot(2,2,1)
xlabel('Test Condition')

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

% get sizes of matrices 
clear indsiz
for k = 1:11
    indsiz(k,11) = size(ind.n11(k).nck,1);
end
for k = 1:9
    indsiz(k,9) = size(ind.n9(k).nck,1);
end
for k = 1:3
    indsiz(k,3) = size(ind.n3(k).nck,1);
end
    
% test model on full set of shifts 
nL = 1; 
nB = 1;
res.test_corr = []; 
for nHarm = 1:3
    for nAtt = 1:3
        for nPart = 1:12
            % DIFFERENTIATE 3 CASES!! 
            x = res.all.sim(nAlpha).dat(nHarm,nAtt).x(nPart, :)';
            Lx = length(x);
            for k = 1:Lx % different lengths of the data depending on condition
                for nInd = 1:indsiz(k,Lx) % check all different combinations of index sets given k
                    alpha_ind = res.max_pred.ind(nPart,nHarm,nAtt,k,nInd);  % take alpha that is best on a constrained training set 
                    x = res.all.sim(alpha_ind).dat(nHarm,nAtt).x(nPart, :)';
                    y = res.all.emp(nB).dat(nHarm,nAtt).y(nPart, :)';
                    res.test_corr(nHarm,nAtt, k, nInd, nPart) = corr(x,y); % compare simulation with empirical data from individual participants 
                %res.test_rms(nHarm,nAtt, k, nInd, nPart) = rms(x-y); % compare simulation with empirical data from individual participants 
                end
            end
        end
    end
end
X = res.test_corr.^2; 
X(X == 0) = nan; 


%% plot 
figure; 
nPl = 0; 
for nHarm = 1:3
    for nAtt = 1:3
        nPl = nPl + 1;
        subplot(3,3,nPl)
        plot([2:11]', squeeze(nanmean(X(nHarm,nAtt, 2:end, :, :),4)), 'color', [.5 .5 .5], 'linewidth', 1); hold on
        plot([2:11]', squeeze(mean(nanmean(X(nHarm,nAtt, 2:end, :, :),4),5)), 'color', [.2 .2 .8, .5], 'linewidth', 4)
        ylim([0 1])
        xlim([1 12])
        set(gca, 'Box', 'off')
        set(gca, 'XTick', [2:1:11])
        set(gca, 'XTickLabel', {'2', '', '4', '', '6', '', '8', '', '10', ''})
    end
end
subplot(3,3,8)
xlabel('# Fitted Conditions')

subplot(3,3,7)
legend('Indiv. data', 'Mean')

subplot(3,3,1)
title('SFS')
subplot(3,3,2)
title('SE')
subplot(3,3,3)
title('SFS-SE')

subplot(3,3,1)
ylabel('Harm')
text(-2, 1.2, 'A', 'fontsize', 18, 'fontweight', 'bold')

subplot(3,3,4)
ylabel('Inharm')
subplot(3,3,7)
ylabel('Log-eq')



