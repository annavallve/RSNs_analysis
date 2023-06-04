%% NETS ANALYSIS ON 100 ICs

% add tolboxes
addpath /Users/anna/Desktop/TFG/FSLNets
addpath (genpath('/Users/anna/Documents/MATLAB/github_repo/cbrewer2'));
addpath ('/Users/anna/Documents/MATLAB/colorspace/colorspace');

% define files & data
drdir = '/Users/anna/Desktop/TFG/Results/dr_100ICs/dr_stage1'; % dual_regression directory (only stage 1 needed)
goodnodes='/Users/anna/Desktop/TFG/Results/dr_100ICs/goodnodes_100ICs.txt'; % RSN (no noise)
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
% ========= NETWORK MATRICES COMPUTING ============
% =================================================

netmatAMP = nets_netmats(ts,0,'amp');           % Amplitude data
netmatCORR = nets_netmats(ts,0,'corr');         % Full Correlation data
netmatPCORR = nets_netmats(ts,1,'ridgep', 1);   % Partial Correlation data
netmatCOV = nets_netmats(ts,0,'cov');           % Covariance data


% =============================================================
% ========= PAIRED DIF IN AMPLITUDES ==========================
% =============================================================

% separation of scan 1 (pre) and scan 2 (post) of each subject
% separation of scan 1 (pre) and scan 2 (post) of each subject
netmatAMP_pre = netmatAMP(1:N,:);
netmatAMP_post = netmatAMP(N+1:end,:);

% mean value of scan 1 amplitudes (pre) and scan 2 amplitudes (post)
mean_AMP_pre = mean(netmatAMP_pre);
mean_AMP_post = mean(netmatAMP_post);

mean_dif_AMP = (mean_AMP_post)-(mean_AMP_pre); % comparison scan 1 and scan 2

% evaluation of ALL amplitudes (NO t-test)
figure('Color', 'white', 'Position', [0 0 1000 200])
imagesc(mean_dif_AMP, [-0.1 0.1]); colormap (cmap); title ('Paired differences 100 ICs (AMPpost - AMPpre)', 'fontsize', 20)
set(gca, 'YTickLabel', '');


% =============================================================
% ========= PAIRED DIF IN CORR ==========================
% =============================================================

% separation of scan 1 (pre) and scan 2 (post) of each subject
netmatCORR_pre = netmatCORR(1:N,:);
netmatCORR_post = netmatCORR(N+1:end,:);

mean_CORR_pre = mean(netmatCORR_pre);
mean_CORR_post = mean(netmatCORR_post);

figure ('Position', [0 0 1200 500], 'Color', 'white')
subplot (1,2,1); imagesc(reshape(mean_CORR_pre,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title (' Mean Correlation scan 1; 100 ICs')
subplot (1,2,2); imagesc(reshape(mean_CORR_post,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title (' Mean Correlation scan 2; 100 ICs')


% =============================================================
% ========= PAIRED DIF IN PCORR ==========================
% =============================================================

% separation of scan 1 (pre) and scan 2 (post) of each subject
netmatPCORR_pre = netmatPCORR(1:N,:);
netmatPCORR_post = netmatPCORR(N+1:end,:);

mean_PCORR_pre = mean(netmatPCORR_pre);
mean_PCORR_post = mean(netmatPCORR_post);

figure ('Position', [0 0 1200 500], 'Color', 'white')
subplot (1,2,1); imagesc(reshape(mean_PCORR_pre,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title (' Mean Partial Correlation scan 1; 100 ICs')
subplot (1,2,2); imagesc(reshape(mean_PCORR_post,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title (' Mean Partial Correlation scan 2; 100 ICs')


% =============================================================
% ========= PAIRED DIF IN COV ==========================
% =============================================================

% separation of scan 1 (pre) and scan 2 (post) of each subject
netmatCOV_pre = netmatCOV(1:N,:);
netmatCOV_post = netmatCOV(N+1:end,:);

mean_COV_pre = mean(netmatCOV_pre);
mean_COV_post = mean(netmatCOV_post);

figure ('Position', [0 0 1200 500], 'Color', 'white')
subplot (1,2,1); imagesc(reshape(mean_COV_pre,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title (' Mean Covariance scan 1; 100 ICs')
subplot (1,2,2); imagesc(reshape(mean_COV_post,N_goodnodes,N_goodnodes), [-1 1]); colormap (cmap); title (' Mean Covariance scan 2; 100 ICs')

