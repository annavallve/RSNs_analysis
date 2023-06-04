%% NETS ANALYSIS ON 50 ICs WITH ORDER NETWORKS

% add tolboxes
addpath /Users/anna/Desktop/TFG/FSLNets
addpath (genpath('/Users/anna/Documents/MATLAB/github_repo/cbrewer2'));
addpath ('/Users/anna/Documents/MATLAB/colorspace/colorspace');

% define files & data
drdir = '/Users/anna/Desktop/TFG/Results/dr_50ICs/dr_stage1'; % dual_regression directory (only stage 1 needed)
goodnodes = '/Users/anna/Desktop/TFG/Results/dr_50ICs/goodnodes_50ICs.txt'; % RSN (no noise)
N = 44; % number of subjects

% define colormaps
cmap = cbrewer2 ('div', 'RdBu', 100); 
cmap = cmap(end:-1:1, :);


% =================================================
% ========= TIMESERIES LOADING ====================
% =================================================
ts = nets_load (drdir, 0.8, 1); % 0.8 = time between data acquisition
                                % 1 = normalise whole subject stddev
ts.DD = load(goodnodes);        % remove noise components
ts = nets_tsclean (ts, 1);      % 1 to regress the bad timeseries
tp = ts.NtimepointsPerSubject;  % data timepoints per subject

timeseries_ALL = ts.ts;         % time series of the goodnodes (as a table)
N_goodnodes = ts.Nnodes;        % number of goodnodes


% =================================================
% ========= NETWORKS ORDERING =====================
% =================================================
ts.ts = zeros(tp*N*2, N_goodnodes);
ts.ts(:,1) = timeseries_ALL(:,1);

% ordering the goodnodes in the following groups:
% (1) Visual; (2) Sensory; (3) Parietal; (4) Executive; (5) DMN; (6) Cerebellum
ts.ts(:,1:N_goodnodes) = timeseries_ALL(:,[4, 7, 20, 2, 19, 14, 6, 12, 22, 10, 21, 5, 11, 25, 16, 17, 1, 15, 23, 24, 3, 8, 9, 13, 18, 26]);

% number of subnetworks inside each RSN
net_sum = [3, 6, 5, 6, 4, 2];        % absolute frequency
net_cum = [0, 3, 9, 14, 20, 24, 26]; % cumulative absolute frequency (starting in 0)


% =================================================
% ========= COMPUTE NETS ==========================
% =================================================

netmatAMP = nets_netmats(ts,0,'amp');           % Amplitude data
netmatCORR = nets_netmats(ts,0,'corr');         % Full Correlation data
netmatPCORR = nets_netmats(ts,1,'ridgep', 1);   % Partial Correlation data
netmatCOV = nets_netmats(ts,0,'cov');           % Covariance data


% =============================================================
% ========= PAIRED DIF IN AMPLITUDES ==========================
% =============================================================

netmatAMP_pre = netmatAMP(1:N,:);
netmatAMP_post = netmatAMP(N+1:end,:);

mean_AMP_pre = mean(netmatAMP_pre);
mean_AMP_post = mean(netmatAMP_post);
mean_dif_AMP = (mean_AMP_post)-(mean_AMP_pre);

% plot of ordered amplitudes with NO t-test
figure('Color', 'white', 'Position', [0 0 1000 200])
imagesc(mean_dif_AMP, [-0.1 0.1]); colormap (cmap); title ('Ordered paired differences (POST - PRE) for 50 ICs', 'fontsize', 20)
set(gca, 'YTickLabel', '');

% t-test performance to get only the significant data (p-value < 0.05)
p_AMP = zeros(1,ts.Nnodes);
for i=1:ts.Nnodes
    [h,p] = ttest(netmatAMP_pre(:,i),netmatAMP_post(:,i));
    p_AMP(i)=p;
end

mean_dif_AMP(p_AMP > 0.05 | p_AMP == 0.05) = 0; % Only variables with p-value < 0.05 will
                                                % be significant. Variables with
                                                % p-value => 0.05 will be set to 0

% amplitude plot of ordered networks with t-test
figure('Color', 'white', 'Position', [0 0 1000 200])
imagesc(mean_dif_AMP, [-0.1 0.1]); colormap (cmap); title ('Ordered paired differences (POST - PRE) for 50 ICs', 'fontsize', 20)
set(gca, 'YTickLabel', '');

% calculation of amplitude changes as percentages
for i=1:length(net_sum)
    n = sum(mean_dif_AMP(1+net_cum(i):net_cum(i+1)))/sum(mean_AMP_pre(1+net_cum(i):net_cum(i+1)))*100;
    percent_AMP(i) = n;
end

figure('Position', [300, 300, 300, 700], 'Color', 'White')
imagesc(percent_AMP', [-10 10]); colormap (cmap); title ('Amplitude Changes (%)', 'fontsize', 15)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.55, 1, 'Visual', 'color', 'white', 'fontsize', 15); text(0.55, 2, 'Sensory', 'fontsize', 15, 'color', 'white');
text(0.55, 3, 'Parietal', 'fontsize', 15, 'color', 'white'); text(0.55, 4, 'Executive', 'fontsize', 15, 'color', 'white');
text(0.55, 5, 'DMN', 'fontsize', 15, 'color', 'white'); text(0.55, 6, 'Cerebellum', 'fontsize', 15, 'color', 'white')


% =============================================================
% ========= PAIRED DIF IN CORR ==========================
% =============================================================
netmatCORR_pre = netmatCORR(1:N,:);
netmatCORR_post = netmatCORR(N+1:end,:);
mean_CORR_pre = mean(netmatCORR_pre);
mean_CORR_post = mean(netmatCORR_post);

% correlation and anticorrelation plot of ordered networks 
figure ('Position', [0 0 1200 500], 'Color', 'white')
subplot (1,2,1); imagesc(reshape(mean_CORR_pre,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title ('Mean Correlation 1st scan')
subplot (1,2,2); imagesc(reshape(mean_CORR_post,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title ('Mean Correlation 2nd scan')


% t-test performance to get only the significant data
p_CORR = zeros(1,ts.Nnodes*ts.Nnodes);
for i=1:ts.Nnodes*ts.Nnodes
    [h,p] = ttest(netmatCORR_pre(:,i),netmatCORR_post(:,i));
    p_CORR(i)=p;
end

p_CORR_sig = zeros(1, N_goodnodes*N_goodnodes);
p_CORR_sig(p_CORR<0.05) = 1;

% Differentiation between CORRELATION AND ANTICORRELATION

% ========= CORRELATION =========
% first, get only correlation -> set anticorrelated values to 0
mean_CORR_pre_corr = mean_CORR_pre;
mean_CORR_pre_corr(mean_CORR_pre < 0) = 0;
mean_CORR_pre_corr(mean_CORR_post < 0) = 0;

mean_CORR_post_corr = mean_CORR_post;
mean_CORR_post_corr(mean_CORR_pre < 0) = 0;
mean_CORR_post_corr(mean_CORR_post < 0) = 0;

mean_dif_CORR_corr = mean_CORR_post_corr - mean_CORR_pre_corr;
mat_dif_CORR_corr = reshape(mean_dif_CORR_corr, N_goodnodes, N_goodnodes);
mat_pre_CORR_corr = reshape(mean_CORR_pre_corr, N_goodnodes, N_goodnodes);

% set unsignificant values to 0 for matrices of mean-pre-correlation and
% mean-difference-correlation
p_CORR_sig_corr = p_CORR_sig.*mean_dif_CORR_corr;
mat_dif_CORR_corr(p_CORR>0.05 | p_CORR==0.05) = 0;
mat_dif_CORR_corr(p_CORR>0.05 | p_CORR==0.05) = 0;

% correlation mean diffference plot only with significant values and for
% all voxels
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(reshape(p_CORR_sig_corr,N_goodnodes,N_goodnodes), [-0.1 0.1]), colormap (cmap), title ('CORRELATION PAIRED DIFFERENCES (p<0.05)', 'fontsize', 20)
hold on
xline(3.5); xline(9.5); xline(14.5); xline(20.5); xline(24.5)
yline(3.5); yline(9.5); yline(14.5); yline(20.5); yline(24.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(1.3, 27, 'Visual', 'fontsize', 15); text(5, 27, 'Sensory', 'fontsize', 15); text(11, 27, 'Parietal', 'fontsize', 15)
text(16, 27, 'Executive', 'fontsize', 15); text(21.5, 27, 'DMN', 'fontsize', 15); text(24.7, 27, 'Cerebellum', 'fontsize', 15)
text(-3, 2, 'Visual', 'fontsize', 15); text(-3, 6.5, 'Sensory', 'fontsize', 15); text(-3, 12, 'Parietal', 'fontsize', 15)
text(-3, 17.5, 'Executive', 'fontsize', 15); text(-3, 22.5, 'DMN', 'fontsize', 15); text(-3, 25.5, 'Cerebellum', 'fontsize', 15)

for i=1:length(net_sum)
    for j=1:length(net_sum)
        n = sum(mat_dif_CORR_corr(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')/sum(mat_pre_CORR_corr(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')*100;
        percent_CORR_corr(i, j) = n;
    end
end

percent_CORR_corr(isnan(percent_CORR_corr)) = 0;

% correlation mean difference plot as percentage for each network
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(percent_CORR_corr, [-10 10]); colormap (cmap); title ('Correlation Changes', 'fontsize', 20)
hold on
xline(1.5); xline(2.5); xline(3.5); xline(4.5); xline(5.5)
yline(1.5); yline(2.5); yline(3.5); yline(4.5); yline(5.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.8, 6.7, 'Visual', 'fontsize', 15); text(1.8, 6.7, 'Sensory', 'fontsize', 15); text(2.8, 6.7, 'Parietal', 'fontsize', 15)
text(3.8, 6.7, 'Executive', 'fontsize', 15); text(4.8, 6.7, 'DMN', 'fontsize', 15); text(5.7, 6.7, 'Cerebellum', 'fontsize', 15)
text(-0.3, 1, 'Visual', 'fontsize', 15); text(-0.3, 2, 'Sensory', 'fontsize', 15); text(-0.3, 3, 'Parietal', 'fontsize', 15)
text(-0.3, 4, 'Executive', 'fontsize', 15); text(-0.3, 5, 'DMN', 'fontsize', 15); text(-0.3, 6, 'Cerebellum', 'fontsize', 15)


% ========= ANTICORRELATION =========
% second, get only anticorrelation -> set correlated values to 0
mean_CORR_pre_anticorr = mean_CORR_pre;
mean_CORR_pre_anticorr(mean_CORR_pre > 0) = 0;
mean_CORR_pre_anticorr(mean_CORR_post > 0) = 0;

mean_CORR_post_anticorr = mean_CORR_post;
mean_CORR_post_anticorr(mean_CORR_pre > 0) = 0;
mean_CORR_post_anticorr(mean_CORR_post > 0) = 0;

mean_dif_CORR_anticorr = mean_CORR_post_anticorr - mean_CORR_pre_anticorr;
mat_dif_CORR_anticorr = reshape(mean_dif_CORR_anticorr, N_goodnodes, N_goodnodes);
mat_pre_CORR_anticorr = reshape(mean_CORR_pre_anticorr, N_goodnodes, N_goodnodes);

% set unsignificant values to 0 for matrices of mean-pre-correlation and
% mean-difference-correlation
p_CORR_sig_anticorr = p_CORR_sig.*mean_dif_CORR_anticorr;
mat_dif_CORR_anticorr(p_CORR>0.05 | p_CORR==0.05) = 0;
mat_pre_CORR_anticorr(p_CORR>0.05 | p_CORR==0.05) = 0;

% anticorrelation mean difference plot only with significant values and for
% all voxels
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(reshape(p_CORR_sig_anticorr,N_goodnodes,N_goodnodes), [-0.1 0.1]), colormap (cmap), title ('ANTICORRELATION PAIRED DIFFERENCES (p<0.05)', 'fontsize', 20)
hold on
xline(3.5); xline(9.5); xline(14.5); xline(20.5); xline(24.5)
yline(3.5); yline(9.5); yline(14.5); yline(20.5); yline(24.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(1.3, 27, 'Visual', 'fontsize', 15); text(5, 27, 'Sensory', 'fontsize', 15); text(11, 27, 'Parietal', 'fontsize', 15)
text(16, 27, 'Executive', 'fontsize', 15); text(21.5, 27, 'DMN', 'fontsize', 15); text(24.7, 27, 'Cerebellum', 'fontsize', 15)
text(-3, 2, 'Visual', 'fontsize', 15); text(-3, 6.5, 'Sensory', 'fontsize', 15); text(-3, 12, 'Parietal', 'fontsize', 15)
text(-3, 17.5, 'Executive', 'fontsize', 15); text(-3, 22.5, 'DMN', 'fontsize', 15); text(-3, 25.5, 'Cerebellum', 'fontsize', 15)

for i=1:length(net_sum)
    for j=1:length(net_sum)
        n = sum(mat_dif_CORR_anticorr(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')/abs(sum(mat_pre_CORR_anticorr(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all'))*100;
        percent_CORR_anticorr(i, j) = n;
    end
end

percent_CORR_anticorr(isnan(percent_CORR_anticorr)) = 0;

% anticorrelation mean difference plot as percentage for each network
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(percent_CORR_anticorr, [-350 350]); colormap (cmap); title ('Anticorrelation Changes', 'fontsize', 20)
hold on
xline(1.5); xline(2.5); xline(3.5); xline(4.5); xline(5.5)
yline(1.5); yline(2.5); yline(3.5); yline(4.5); yline(5.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.8, 6.7, 'Visual', 'fontsize', 15); text(1.8, 6.7, 'Sensory', 'fontsize', 15); text(2.8, 6.7, 'Parietal', 'fontsize', 15)
text(3.8, 6.7, 'Executive', 'fontsize', 15); text(4.8, 6.7, 'DMN', 'fontsize', 15); text(5.7, 6.7, 'Cerebellum', 'fontsize', 15)
text(-0.3, 1, 'Visual', 'fontsize', 15); text(-0.3, 2, 'Sensory', 'fontsize', 15); text(-0.3, 3, 'Parietal', 'fontsize', 15)
text(-0.3, 4, 'Executive', 'fontsize', 15); text(-0.3, 5, 'DMN', 'fontsize', 15); text(-0.3, 6, 'Cerebellum', 'fontsize', 15)


% =============================================================
% ========= PAIRED DIF IN PCORR ==========================
% =============================================================
netmatPCORR_pre = netmatPCORR(1:N,:);
netmatPCORR_post = netmatPCORR(N+1:end,:);
mean_PCORR_pre = mean(netmatPCORR_pre);
mean_PCORR_post = mean(netmatPCORR_post);

% correlation and anticorrelation plot of ordered networks 
figure ('Position', [0 0 1200 500], 'Color', 'white')
subplot (1,2,1); imagesc(reshape(mean_PCORR_pre,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title ('Mean Partial Correlation 1st scan')
subplot (1,2,2); imagesc(reshape(mean_PCORR_post,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title ('Mean Partial Correlation 2nd scan')

% t-test performance to get only the significant data
p_PCORR = zeros(1,ts.Nnodes*ts.Nnodes);
for i=1:ts.Nnodes*ts.Nnodes
    [h,p] = ttest(netmatPCORR_pre(:,i),netmatPCORR_post(:,i));
    p_PCORR(i)=p;
end

p_PCORR_sig = zeros(1, N_goodnodes*N_goodnodes);
p_PCORR_sig(p_PCORR<0.05) = 1;

% Differentiation between PARTIAL CORRELATION AND PARTIAL ANTICORRELATION 

% ========= PARTIAL CORRELATION =========
% first, get only partial correlation -> set partial anticorrelated values to 0
mean_PCORR_pre_corr = mean_PCORR_pre;
mean_PCORR_pre_corr(mean_PCORR_pre < 0) = 0;
mean_PCORR_pre_corr(mean_PCORR_post < 0) = 0;

mean_PCORR_post_corr = mean_PCORR_post;
mean_PCORR_post_corr(mean_PCORR_pre < 0) = 0;
mean_PCORR_post_corr(mean_PCORR_post < 0) = 0;

mean_dif_PCORR_corr = mean_PCORR_post_corr - mean_PCORR_pre_corr;
mat_dif_PCORR_corr = reshape(mean_dif_PCORR_corr, N_goodnodes, N_goodnodes);
mat_pre_PCORR_corr = reshape(mean_PCORR_pre_corr, N_goodnodes, N_goodnodes);

% set unsignificant values to 0 for matrices of mean-pre-correlation and
% mean-difference-correlation
p_PCORR_sig_corr = p_PCORR_sig.*mean_dif_PCORR_corr;
mat_dif_PCORR_corr(p_PCORR>0.05 | p_PCORR==0.05) = 0;
mat_dif_PCORR_corr(p_PCORR>0.05 | p_PCORR==0.05) = 0;

% correlation mean diffference plot only with significant values and for
% all voxels
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(reshape(p_PCORR_sig_corr,N_goodnodes,N_goodnodes), [-0.5 0.5]), colormap (cmap), title ('PARTIAL CORRELATION PAIRED DIFFERENCES (p<0.05)', 'fontsize', 20)
hold on
xline(3.5); xline(9.5); xline(14.5); xline(20.5); xline(24.5)
yline(3.5); yline(9.5); yline(14.5); yline(20.5); yline(24.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(1.3, 27, 'Visual', 'fontsize', 15); text(5, 27, 'Sensory', 'fontsize', 15); text(11, 27, 'Parietal', 'fontsize', 15)
text(16, 27, 'Executive', 'fontsize', 15); text(21.5, 27, 'DMN', 'fontsize', 15); text(24.7, 27, 'Cerebellum', 'fontsize', 15)
text(-3, 2, 'Visual', 'fontsize', 15); text(-3, 6.5, 'Sensory', 'fontsize', 15); text(-3, 12, 'Parietal', 'fontsize', 15)
text(-3, 17.5, 'Executive', 'fontsize', 15); text(-3, 22.5, 'DMN', 'fontsize', 15); text(-3, 25.5, 'Cerebellum', 'fontsize', 15)

for i=1:length(net_sum)
    for j=1:length(net_sum)
        n = sum(mat_dif_PCORR_corr(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')/sum(mat_pre_PCORR_corr(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')*100;
        percent_PCORR_corr(i, j) = n;
    end
end

percent_PCORR_corr(isnan(percent_PCORR_corr)) = 0;

% partial correlation mean difference plot as percentage for each network
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(percent_PCORR_corr, [-10 10]); colormap (cmap); title ('Partial Correlation Changes', 'fontsize', 20)
hold on
xline(1.5); xline(2.5); xline(3.5); xline(4.5); xline(5.5)
yline(1.5); yline(2.5); yline(3.5); yline(4.5); yline(5.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.8, 6.7, 'Visual', 'fontsize', 15); text(1.8, 6.7, 'Sensory', 'fontsize', 15); text(2.8, 6.7, 'Parietal', 'fontsize', 15)
text(3.8, 6.7, 'Executive', 'fontsize', 15); text(4.8, 6.7, 'DMN', 'fontsize', 15); text(5.7, 6.7, 'Cerebellum', 'fontsize', 15)
text(-0.3, 1, 'Visual', 'fontsize', 15); text(-0.3, 2, 'Sensory', 'fontsize', 15); text(-0.3, 3, 'Parietal', 'fontsize', 15)
text(-0.3, 4, 'Executive', 'fontsize', 15); text(-0.3, 5, 'DMN', 'fontsize', 15); text(-0.3, 6, 'Cerebellum', 'fontsize', 15)


% ========= PARTIAL ANTICORRELATION =========
% second, get only partial anticorrelation -> set partial correlated values to 0
mean_PCORR_pre_anticorr = mean_PCORR_pre;
mean_PCORR_pre_anticorr(mean_PCORR_pre > 0) = 0;
mean_PCORR_pre_anticorr(mean_PCORR_post > 0) = 0;

mean_PCORR_post_anticorr = mean_PCORR_post;
mean_PCORR_post_anticorr(mean_PCORR_pre > 0) = 0;
mean_PCORR_post_anticorr(mean_PCORR_post > 0) = 0;

mean_dif_PCORR_anticorr = mean_PCORR_post_anticorr - mean_PCORR_pre_anticorr;
mat_dif_PCORR_anticorr = reshape(mean_dif_PCORR_anticorr, N_goodnodes, N_goodnodes);
mat_pre_PCORR_anticorr = reshape(mean_PCORR_pre_anticorr, N_goodnodes, N_goodnodes);

% set unsignificant values to 0 for matrices of mean-pre-partial anticorrelation
% and mean-difference-partial anticorrelation
p_PCORR_sig_anticorr = p_PCORR_sig.*mean_dif_PCORR_anticorr;
mat_dif_PCORR_anticorr(p_PCORR>0.05 | p_PCORR==0.05) = 0;
mat_pre_PCORR_anticorr(p_PCORR>0.05 | p_PCORR==0.05) = 0;

% partial anticorrelation mean difference plot only with significant values
% and for all voxels
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(reshape(p_PCORR_sig_anticorr,N_goodnodes,N_goodnodes), [-0.5 0.5]), colormap (cmap), title ('PARTIAL ANTICORRELATION PAIRED DIFFERENCES (p<0.05)', 'fontsize', 20)
hold on
xline(3.5); xline(9.5); xline(14.5); xline(20.5); xline(24.5)
yline(3.5); yline(9.5); yline(14.5); yline(20.5); yline(24.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(1.3, 27, 'Visual', 'fontsize', 15); text(5, 27, 'Sensory', 'fontsize', 15); text(11, 27, 'Parietal', 'fontsize', 15)
text(16, 27, 'Executive', 'fontsize', 15); text(21.5, 27, 'DMN', 'fontsize', 15); text(24.7, 27, 'Cerebellum', 'fontsize', 15)
text(-3, 2, 'Visual', 'fontsize', 15); text(-3, 6.5, 'Sensory', 'fontsize', 15); text(-3, 12, 'Parietal', 'fontsize', 15)
text(-3, 17.5, 'Executive', 'fontsize', 15); text(-3, 22.5, 'DMN', 'fontsize', 15); text(-3, 25.5, 'Cerebellum', 'fontsize', 15)

for i=1:length(net_sum)
    for j=1:length(net_sum)
        n = sum(mat_dif_PCORR_anticorr(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')/abs(sum(mat_pre_PCORR_anticorr(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all'))*100;
        percent_PCORR_anticorr(i, j) = n;
    end
end

percent_PCORR_anticorr(isnan(percent_PCORR_anticorr)) = 0;

% partial anticorrelation mean difference plot as percentage for each network
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(percent_PCORR_anticorr, [-350 350]); colormap (cmap); title ('Partial Anticorrelation Changes', 'fontsize', 20)
hold on
xline(1.5); xline(2.5); xline(3.5); xline(4.5); xline(5.5)
yline(1.5); yline(2.5); yline(3.5); yline(4.5); yline(5.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.8, 6.7, 'Visual', 'fontsize', 15); text(1.8, 6.7, 'Sensory', 'fontsize', 15); text(2.8, 6.7, 'Parietal', 'fontsize', 15)
text(3.8, 6.7, 'Executive', 'fontsize', 15); text(4.8, 6.7, 'DMN', 'fontsize', 15); text(5.7, 6.7, 'Cerebellum', 'fontsize', 15)
text(-0.3, 1, 'Visual', 'fontsize', 15); text(-0.3, 2, 'Sensory', 'fontsize', 15); text(-0.3, 3, 'Parietal', 'fontsize', 15)
text(-0.3, 4, 'Executive', 'fontsize', 15); text(-0.3, 5, 'DMN', 'fontsize', 15); text(-0.3, 6, 'Cerebellum', 'fontsize', 15)


% =============================================================
% ========= PAIRED DIF IN COV ==========================
% =============================================================

netmatCOV_pre = netmatCOV(1:N,:);
netmatCOV_post = netmatCOV(N+1:end,:);
mean_COV_pre = mean(netmatCOV_pre);
mean_COV_post = mean(netmatCOV_post);

% positive and negative covariance plot of ordered networks 
figure ('Position', [0 0 1200 500], 'Color', 'white')
subplot (1,2,1); imagesc(reshape(mean_COV_pre,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title ('Mean Covariance 1st scan')
subplot (1,2,2); imagesc(reshape(mean_COV_post,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title ('Mean Covariance 2nd scan')


% t-test performance to get only the significant data
p_COV = zeros(1,ts.Nnodes*ts.Nnodes);
for i=1:ts.Nnodes*ts.Nnodes
    [h,p] = ttest(netmatCOV_pre(:,i),netmatCOV_post(:,i));
    p_COV(i)=p;
end

p_COV_sig = zeros(1, N_goodnodes*N_goodnodes);
p_COV_sig(p_COV<0.05) = 1;

% Differentiation between POSITIVE COVARIANCE AND NEGATIVE COVARIANCE

% ========= POSITIVE COVARIANCE =========
% first, get only positive covariance values -> set negative covariance values to 0
mean_COV_pre_p = mean_COV_pre;
mean_COV_pre_p(mean_COV_pre < 0) = 0;
mean_COV_pre_p(mean_COV_post < 0) = 0;

mean_COV_post_p = mean_COV_post;
mean_COV_post_p(mean_COV_pre < 0) = 0;
mean_COV_post_p(mean_COV_post < 0) = 0;

mean_dif_COV_p = mean_COV_post_p - mean_COV_pre_p;

% differentiation between increased and decreased positive covariance
mean_dif_COV_p(mean_dif_COV_p<0) = -0.5;
mean_dif_COV_p(mean_dif_COV_p>0) = 0.5;

mat_dif_COV_p = reshape(mean_dif_COV_p, N_goodnodes, N_goodnodes);
mat_pre_COV_p = reshape(mean_COV_pre_p, N_goodnodes, N_goodnodes);

% set unsignificant values to 0 for matrices of mean-pre-positive-covariance and
% mean-difference-positive-covariance
p_COV_sig_p = p_COV_sig.*mean_dif_COV_p;
mat_dif_COV_p(p_COV>0.05 | p_COV==0.05) = 0;
mat_dif_COV_p(p_COV>0.05 | p_COV==0.05) = 0;

% positive covarianve mean diffference plot only with significant values
% and for all voxels
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(reshape(p_COV_sig_p,N_goodnodes,N_goodnodes), [-1 1]), colormap (cmap), title ('POSITIVE COVARIANCE PAIRED DIFFERENCES (p<0.05)', 'fontsize', 20)
hold on
xline(3.5); xline(9.5); xline(14.5); xline(20.5); xline(24.5)
yline(3.5); yline(9.5); yline(14.5); yline(20.5); yline(24.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(1.3, 27, 'Visual', 'fontsize', 15); text(5, 27, 'Sensory', 'fontsize', 15); text(11, 27, 'Parietal', 'fontsize', 15)
text(16, 27, 'Executive', 'fontsize', 15); text(21.5, 27, 'DMN', 'fontsize', 15); text(24.7, 27, 'Cerebellum', 'fontsize', 15)
text(-3, 2, 'Visual', 'fontsize', 15); text(-3, 6.5, 'Sensory', 'fontsize', 15); text(-3, 12, 'Parietal', 'fontsize', 15)
text(-3, 17.5, 'Executive', 'fontsize', 15); text(-3, 22.5, 'DMN', 'fontsize', 15); text(-3, 25.5, 'Cerebellum', 'fontsize', 15)

% Percentage of increased positive covariance
p_COV_pos_sig_inc(p_COV_sig_p==0.5) = 1;
p_COV_pos_sig_inc(p_COV_sig_p==-0.5) = 0;
p_COV_pos_sig_inc(p_COV_sig_p==0) = 0;
mat_COV_pos_sig_inc = reshape(p_COV_pos_sig_inc, N_goodnodes, N_goodnodes);

for i=1:length(net_sum)
    for j=1:length(net_sum)
        n = sum(mat_COV_pos_sig_inc(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')/(net_sum(j)*net_sum(i))*100;
        percent_COV_p_inc(i, j) = n;
    end
end

percent_COV_p_inc(isnan(percent_COV_p_inc)) = 0;

figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(percent_COV_p_inc, [-50 50]); colormap (cmap); title ('Positive Covariance Changes', 'fontsize', 20)
hold on
xline(1.5); xline(2.5); xline(3.5); xline(4.5); xline(5.5)
yline(1.5); yline(2.5); yline(3.5); yline(4.5); yline(5.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.8, 6.7, 'Visual', 'fontsize', 15); text(1.8, 6.7, 'Sensory', 'fontsize', 15); text(2.8, 6.7, 'Parietal', 'fontsize', 15)
text(3.8, 6.7, 'Executive', 'fontsize', 15); text(4.8, 6.7, 'DMN', 'fontsize', 15); text(5.7, 6.7, 'Cerebellum', 'fontsize', 15)
text(-0.3, 1, 'Visual', 'fontsize', 15); text(-0.3, 2, 'Sensory', 'fontsize', 15); text(-0.3, 3, 'Parietal', 'fontsize', 15)
text(-0.3, 4, 'Executive', 'fontsize', 15); text(-0.3, 5, 'DMN', 'fontsize', 15); text(-0.3, 6, 'Cerebellum', 'fontsize', 15)

% Percentage of decreased positive covariance
p_COV_pos_sig_dec(p_COV_sig_p==0.5) = 0;
p_COV_pos_sig_dec(p_COV_sig_p==-0.5) = -1;
p_COV_pos_sig_dec(p_COV_sig_p==0) = 0;
mat_COV_pos_sig_dec = reshape(p_COV_pos_sig_dec, N_goodnodes, N_goodnodes);

for i=1:length(net_sum)
    for j=1:length(net_sum)
        n = sum(mat_COV_pos_sig_dec(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')/(net_sum(j)*net_sum(i))*100;
        percent_COV_p_dec(i, j) = n;
    end
end

percent_COV_p_dec(isnan(percent_COV_p_dec)) = 0;

figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(percent_COV_p_dec, [-50 50]); colormap (cmap); title ('Positive Covariance Changes', 'fontsize', 20)
hold on
xline(1.5); xline(2.5); xline(3.5); xline(4.5); xline(5.5)
yline(1.5); yline(2.5); yline(3.5); yline(4.5); yline(5.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.8, 6.7, 'Visual', 'fontsize', 15); text(1.8, 6.7, 'Sensory', 'fontsize', 15); text(2.8, 6.7, 'Parietal', 'fontsize', 15)
text(3.8, 6.7, 'Executive', 'fontsize', 15); text(4.8, 6.7, 'DMN', 'fontsize', 15); text(5.7, 6.7, 'Cerebellum', 'fontsize', 15)
text(-0.3, 1, 'Visual', 'fontsize', 15); text(-0.3, 2, 'Sensory', 'fontsize', 15); text(-0.3, 3, 'Parietal', 'fontsize', 15)
text(-0.3, 4, 'Executive', 'fontsize', 15); text(-0.3, 5, 'DMN', 'fontsize', 15); text(-0.3, 6, 'Cerebellum', 'fontsize', 15)


% ========= NEGATIVE COVARIANCE =========
% second, get only negative covariance values -> set positive covariance values to 0
mean_COV_pre_n = mean_COV_pre;
mean_COV_pre_n(mean_COV_pre > 0) = 0;
mean_COV_pre_n(mean_COV_post > 0) = 0;

mean_COV_post_n = mean_COV_post;
mean_COV_post_n(mean_COV_pre > 0) = 0;
mean_COV_post_n(mean_COV_post > 0) = 0;

mean_dif_COV_n = mean_COV_post_n - mean_COV_pre_n;

% differentiation between increased and decreased negative covariance

mean_dif_COV_n(mean_dif_COV_n<0) = -0.5;
mean_dif_COV_n(mean_dif_COV_n>0) = 0.5;

mat_dif_COV_n = reshape(mean_dif_COV_n, N_goodnodes, N_goodnodes);
mat_pre_COV_n = reshape(mean_COV_pre_n, N_goodnodes, N_goodnodes);

% set unsignificant values to 0 for matrices of mean-pre-negative-covariance
% and mean-difference-negative-covariance
p_COV_sig_n = p_COV_sig.*mean_dif_COV_n;
mat_dif_COV_n(p_COV>0.05 | p_COV==0.05) = 0;
mat_pre_COV_n(p_COV>0.05 | p_COV==0.05) = 0;

% negative covariance mean difference plot only with significant values and for
% all voxels
figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(reshape(p_COV_sig_n,N_goodnodes,N_goodnodes), [-1 1]), colormap (cmap), title ('NEGATIVE COVARIANCE PAIRED DIFFERENCES (p<0.05)', 'fontsize', 20)
hold on
xline(3.5); xline(9.5); xline(14.5); xline(20.5); xline(24.5)
yline(3.5); yline(9.5); yline(14.5); yline(20.5); yline(24.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(1.3, 27, 'Visual', 'fontsize', 15); text(5, 27, 'Sensory', 'fontsize', 15); text(11, 27, 'Parietal', 'fontsize', 15)
text(16, 27, 'Executive', 'fontsize', 15); text(21.5, 27, 'DMN', 'fontsize', 15); text(24.7, 27, 'Cerebellum', 'fontsize', 15)
text(-3, 2, 'Visual', 'fontsize', 15); text(-3, 6.5, 'Sensory', 'fontsize', 15); text(-3, 12, 'Parietal', 'fontsize', 15)
text(-3, 17.5, 'Executive', 'fontsize', 15); text(-3, 22.5, 'DMN', 'fontsize', 15); text(-3, 25.5, 'Cerebellum', 'fontsize', 15)

% Percentage of decreased negative covariance (less negative)
p_COV_neg_sig_inc(p_COV_sig_n==0.5) = 1;
p_COV_neg_sig_inc(p_COV_sig_n==-0.5) = 0;
p_COV_neg_sig_inc(p_COV_sig_n==0) = 0;
mat_COV_neg_sig_inc = reshape(p_COV_neg_sig_inc, N_goodnodes, N_goodnodes);

for i=1:length(net_sum)
    for j=1:length(net_sum)
        n = sum(mat_COV_neg_sig_inc(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')/(net_sum(j)*net_sum(i))*100;
        percent_COV_n_inc(i, j) = n;
    end
end

percent_COV_n_inc(isnan(percent_COV_n_inc)) = 0;

figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(percent_COV_n_inc, [-50 50]); colormap (cmap); title ('Negative Covariance Changes (less negative)', 'fontsize', 20)
hold on
xline(1.5); xline(2.5); xline(3.5); xline(4.5); xline(5.5)
yline(1.5); yline(2.5); yline(3.5); yline(4.5); yline(5.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.8, 6.7, 'Visual', 'fontsize', 15); text(1.8, 6.7, 'Sensory', 'fontsize', 15); text(2.8, 6.7, 'Parietal', 'fontsize', 15)
text(3.8, 6.7, 'Executive', 'fontsize', 15); text(4.8, 6.7, 'DMN', 'fontsize', 15); text(5.7, 6.7, 'Cerebellum', 'fontsize', 15)
text(-0.3, 1, 'Visual', 'fontsize', 15); text(-0.3, 2, 'Sensory', 'fontsize', 15); text(-0.3, 3, 'Parietal', 'fontsize', 15)
text(-0.3, 4, 'Executive', 'fontsize', 15); text(-0.3, 5, 'DMN', 'fontsize', 15); text(-0.3, 6, 'Cerebellum', 'fontsize', 15)

% Percentage of increased negative covariance (more negative)
p_COV_neg_sig_dec(p_COV_sig_n==0.5) = 0;
p_COV_neg_sig_dec(p_COV_sig_n==-0.5) = -1;
p_COV_neg_sig_dec(p_COV_sig_n==0) = 0;
mat_COV_neg_sig_dec = reshape(p_COV_neg_sig_dec, N_goodnodes, N_goodnodes);

for i=1:length(net_sum)
    for j=1:length(net_sum)
        n = sum(mat_COV_neg_sig_dec(1+net_cum(i):net_cum(i+1),1+net_cum(j):net_cum(j+1)), 'all')/(net_sum(j)*net_sum(i))*100;
        percent_COV_n_dec(i, j) = n;
    end
end

percent_COV_n_dec(isnan(percent_COV_n_dec)) = 0;

figure('Position', [300, 300, 800, 800], 'Color', 'White')
imagesc(percent_COV_n_dec, [-50 50]); colormap (cmap); title ('Negative Covariance Changes (more negative)', 'fontsize', 20)
hold on
xline(1.5); xline(2.5); xline(3.5); xline(4.5); xline(5.5)
yline(1.5); yline(2.5); yline(3.5); yline(4.5); yline(5.5)
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
text(0.8, 6.7, 'Visual', 'fontsize', 15); text(1.8, 6.7, 'Sensory', 'fontsize', 15); text(2.8, 6.7, 'Parietal', 'fontsize', 15)
text(3.8, 6.7, 'Executive', 'fontsize', 15); text(4.8, 6.7, 'DMN', 'fontsize', 15); text(5.7, 6.7, 'Cerebellum', 'fontsize', 15)
text(-0.3, 1, 'Visual', 'fontsize', 15); text(-0.3, 2, 'Sensory', 'fontsize', 15); text(-0.3, 3, 'Parietal', 'fontsize', 15)
text(-0.3, 4, 'Executive', 'fontsize', 15); text(-0.3, 5, 'DMN', 'fontsize', 15); text(-0.3, 6, 'Cerebellum', 'fontsize', 15)
