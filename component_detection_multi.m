%% Extract temporal components from spatial map time series
%	This script extracts the temporally recurrent network states from
% spatial or 


%%	SET UP FILE SYSTEM

% Clear workspace
clear; close all; clc

% Shuffle random seed.  Necessary to avoid repeating random seeds across parallel computations.
rng('shuffle');

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{2,1} = strjoin(path{1},'/');
path{1,1} = strjoin(path{1}(1:end-1),'/');

% Set data-specific subdirectories
path{3,1} = fullfile(path{2}, 'Data');
path{4,1} = fullfile(path{2}, 'Results');

% Add relevant paths
fpath{1,1} = fullfile(path{1}, 'MATLAB','spm12');
fpath{2,1} = fullfile(path{1}, 'MATLAB','gift','GroupICAT','icatb');
fpath{3,1} = fullfile(path{1}, 'MATLAB','FastICA');
fpath{4,1} = fullfile(path{1}, 'MATLAB','permutationTest');
fpath{5,1} = fullfile(path{2}, 'Functions');
fpath{6,1} = fullfile(path{1}, 'MATLAB','BCT');
% fpath{5,1} = fullfile(path{2}, 'LEICA', 'Functions');
addpath(fpath{1});
for k = 2:numel(fpath)
	addpath(genpath(fpath{k}));
end
clear fpath k


%% Load data

% Load formatted dFNC data
load(fullfile(path{3}, 'FBIRN_DFNC_table.mat'));

% Load functional network labels
labels = readtable(fullfile(path{3}, 'NeuroMark_FNC_labels.xlsx')); % NeuroMark functional network labels & locations

% Remove borders between functional domains
[r,~] = find(strcmpi(labels{:,:}, ""));   % find rows which separate functional domains
r = unique(r);
coords = array2table(cellfun(@str2num, labels{:,2:end}, 'UniformOutput', false), 'RowNames',labels{:,1}, 'VariableNames',labels.Properties.VariableNames(2:end));
coords(r,:) = [];

% Concatenate dFNC time series (space = rows, time = columns)
ts = cell2mat(DFNC_FBIRN)';


%% Rename & convert table variables (should be placed in separate script)

% Convert group diagnosis code to string(s) & logical array
groups = analysis_data{:,'diagnosis(1:sz; 2:hc)'};
N.conditions = numel(unique(groups));
group = string(groups);
group(groups==2) = 'HC'; group(groups==1) = 'SZ'; clear groups
I = [strcmpi(group, 'SZ'), strcmpi(group, 'HC')];
I = array2table(I, "VariableNames",["SZ", "HC"]);
groups = unique(group, 'stable');

% Convert gender code to string(s)
g = analysis_data{:,'gender(1:male; 2:female)'};
gender = string(g);
gender(g==1) = 'M'; gender(g==2) = 'F';
clear g

% Rename & convert table variables
analysis_data = renamevars(analysis_data, ["diagnosis(1:sz; 2:hc)","gender(1:male; 2:female)"], ["Diagnosis","Gender"]);
analysis_data = convertvars(analysis_data, ["Diagnosis","Gender"], 'string');

% Replace numeric codes with strings
analysis_data{:,'Diagnosis'} = group;
analysis_data{:,'Gender'} = gender;
clear group gender

% Replace numeric missing data code with NaN
analysis_data{:,5:end}(analysis_data{:,5:end} == -9999) = NaN;

% Set figure counter
N.fig = 1;

% Index time, subjects
N.TR = size(DFNC_FBIRN{1},1);
N.subjects = sum(I{:,:},1);

% Set number of neighbors to search for in KNN
co = HShannon_kNN_k_initialization(1);

% Set number of ROIs
N.ROI = size(coords,1);

% Determine pairwise comparisons to make
comps = "all";
if strcmpi(comps, "all")
    comps = nchoosek(1:N.conditions, 2);    % all comparisons
else
    comps = find(matches(groups, comps))';
end
N.comp = size(comps,1);

% Analysis Settings
intercept = false;	% determines whether to use intercept term in NBS

% Set parameter sweeps
tstat = 3:0.5:5;


%%  Find number of ICs

% find maximum number(s) of components to search for
disp("Identifying number independent components from Marcenko-Pasteur distribution.");
N.IC = NumberofIC(ts);

% Set IC counts to use (if preset)
N.IC = [5; 10; 22; N.IC];


%% Define filename based on parameters

% Get number of ICs
fileName = strjoin([strtrim(string(num2str(N.IC))); "ICs"],'_');

% Set iteration number
fList = dir(fullfile(path{4}, strcat(strjoin([fileName, "iteration"], '_'), '*.mat')));	% Get file list
a = false(numel(fList),1);
for n = 1:numel(fList)
    a(n) = matches("iteration", strsplit(fList(n).name, '_'));
end
nIter = numel(fList)-sum(a)+1;

% Set full filename
fileName = strjoin([fileName, strcat("iteration", num2str(nIter))], '_');
clear fList nIter a
clear a k n


%% Visualise static FNC matrix with indices

% subject-level sFNC
sFNC{1} = cell2mat(cellfun(@mean, DFNC_FBIRN, 'UniformOutput', false));

% group-level sFNC
sFNC{2} = nan(numel(I.Properties.VariableNames), size(ts,1));
for k = 1:numel(I.Properties.VariableNames)
    sFNC{2}(k,:) = mean(sFNC{1}(I{:,'SZ'}==(k-1),:));
end

% Convert to square matrices
sFNC{1} = icatb_vec2mat(sFNC{1});
sFNC{2} = icatb_vec2mat(sFNC{2});
N.ROI = size(sFNC{1},2);

% Check if need intercept term
if intercept == true
	I = horzcat(ones(size(I,1),1), I);
end

% Allocate contrast matrices
if size(comps,1) > 1
    cont = zeros((size(comps,1)+size(I,2))*2, size(I,2));
    strcont = cell((size(comps,1)+size(I,2))*2, 1);

    % Set all-way contrasts
    c = 2*(size(comps,1)+1:size(comps,1)+N.conditions);
    cont(c-1, :) = -ones(size(I,2), size(I,2)) + diag(2*ones(N.conditions,1));
    cont(c, :) = ones(size(I,2), size(I,2)) - diag(2*ones(N.conditions,1));
    for c = 2*(size(comps,1)+1:size(comps,1)+N.conditions)
        strcont{c-1} = strjoin([groups((c-2*size(comps,1))/2), '> ALL']);
        strcont{c} = strjoin([groups((c-2*size(comps,1))/2), '< ALL']);
    end
else
    cont = zeros(size(comps,1)*2, size(I,2));
    strcont = cell(size(comps,1)*2, 1);
end

% Set pairwise contrasts
for c = 2*(1:size(comps,1))
    cont(c-1, comps(c/2,:)) = [1 -1];
    cont(c, comps(c/2,:)) = [-1 1];
    strcont{c-1} = strjoin([groups(comps(c/2,1)), '>', groups(comps(c/2,2))]);
    strcont{c} = strjoin([groups(comps(c/2,1)), '<', groups(comps(c/2,2))]);
end
clear c

% Run NBS analysis
sFNC{1} = permute(sFNC{1}, [2 3 1]);
[nbs, STATS, GLM, storarray] = runNBS(sFNC{1}, cont, I, N, tstat);

% Visualise group-level sFNC
F(N.fig) = figure; N.fig = N.fig + 1;
for k = 1:numel(I.Properties.VariableNames)
	ax(1,k) = subplot(1,numel(I.Properties.VariableNames),k);
    imagesc(squeeze(sFNC{2}(k,:,:))); colormap jet;
    clim([-max(abs(sFNC{2}),[],'all'), max(abs(sFNC{2}),[],'all')]); colorbar;
    set(ax(1,k), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;
    pbaspect([1 1 1]);
	title(strjoin(["Mean sFNC of", I.Properties.VariableNames{k}], " "), 'FontSize',16);
end
clear k


%% Extract components & activities from dFNC

% Run PCA decomposition
[~,~,pc.ltnt, pc.tsqrd, pc.expl, pc.mu] = pca(ts');

% extract components & activities from dFNC
clust = nan(N.TR*sum(N.subjects), numel(N.IC));
activity = cell(numel(N.IC),1);
sources = cell(numel(N.IC),1);
A = cell(numel(N.IC),1);
for k = 1:numel(N.IC)
    [whitesig, dewhiteM] = icatb_calculate_pca(ts, N.IC(k));
    [~, ~, A{k}, sources{k}] = icatb_icaAlgorithm(1, whitesig');
    A{k} = dewhiteM*A{k};
    [~, clust(:,k)] = max(A,[],2);
    
    activity{k} = reshape(A{k}, [N.IC(k), N.TR, sum(N.subjects)]);
end
clust = clust';
clear k whitesig dewhiteM

% Split activity, time series by subject

% reshape activities to allow easy plotting
ts_ic = nan(N.IC, N.TR, max(N.subjects), 2);
for c = 1:2
	ts_ic(:,:, 1:nnz(I{:,'HC'} == c-1), c) = activity(:,:,I{:,'HC'} == c-1);
end
ts_ic = reshape(ts_ic, [N.IC, N.TR*max(N.subjects), 2]);

% Cluster evaluations
ev.sil = evalclusters(ts', clust', 'silhouette');
ev.var = evalclusters(ts', clust', 'CalinskiHarabasz');
ev.db = evalclusters(ts', clust', 'DaviesBouldin');

% Plot variance explained per component
F(N.fig) = figure; N.fig = N.fig + 1;
ax(1,1) = subplot(2, 17, 1:8);
bar(pc.expl);
ylim([0 max(pc.expl)]); xlim([0 N.ROI*(N.ROI-1)/2]); hold on;
xlabel("Principal Components");
ylabel("% Variance Captured");
title('% Variance per tIC');

% Plot cumulative variance explained as function of component count
ax(1,2) = subplot(2, 17, 10:17);
plot(cumsum(pc.expl));
ylim([0 100]); xlim([0 N.ROI*(N.ROI-1)/2]); hold on;
t = plot([N.IC N.IC]', repmat(ax(1,2).YLim,[numel(N.IC) 1])', '--', 'LineWidth',0.25);
xlabel("Principal Components");
ylabel("% Cumulative Variance Captured");
legend(t, strcat([repmat('t = ',[numel(N.IC),1]), num2str(N.IC)]));
title("% Variance as Function of tIC Count");

% Plot cluster evaluations
ax(2,1) = subplot(2, 17, 17+[1:5]); plot(ev.sil); title("Silhouette Scores");
ax(2,2) = subplot(2, 17, 17+[7:11]); plot(ev.var); title("Calinski-Harabasz Scores");
ax(2,3) = subplot(2, 17, 17+[13:17]); plot(ev.db); title("Davies-Bouldin Scores");
clear ax t k clust activations


%% Isolate components & activity from dFNC

% Plot sample connectivity maps
F(N.fig) = figure; N.fig = N.fig + 1;
N.samplot = 4;
C = nchoosek(1:N.IC, N.samplot);
C = C(ceil(rand*size(C,1)),:);
fncmat = nan(N.ROI, N.ROI, N.samplot);
for j = 1:N.samplot
    % Plot connectivity maps
	ax(j,1) = subplot(2, N.samplot, j); hold on
    fncmat(:,:,j) = icatb_vec2mat(sources(C(j),:)');
	imagesc(fncmat(:,:,j)); colormap jet;
    clim([-max(abs(sources),[],'all'), max(abs(sources),[],'all')]); colorbar;
    set(ax(j,1), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;
    pbaspect([1 1 1]);
	title(strjoin(["Connectivity Matrix of Component", num2str(C(j))], " "), 'FontSize',16);
	ylabel("Spatial Networks");
    xlabel("Spatial Networks");

    % Plot time courses
    ax(j,2) = subplot(2, N.samplot, j+N.samplot); hold on
    plot(1:N.TR*max(N.subjects), squeeze(ts_ic(C(j),:,:)));
    title(strjoin(["Time Course of Component", num2str(C(j))], " "), 'FontSize',16);
    set(ax(j,2), {'XTick', 'XTickLabels', 'XLim'}, {0:N.TR:max(N.subjects)*N.TR, [], [0 max(N.subjects)*N.TR]}); hold on;
    legend(I.Properties.VariableNames);
    ylabel("Activity");
    xlabel("Subjects");
end
clear j c


%% Calculate subject-level entropy & joint entropy

% Find component entropies
e.comp = nan(N.IC, max(sum(I{:,'SZ'}), sum(N.subjects)-sum(I{:,'SZ'})), 2);
for c = 1:2
    a = activity(:,:,I{:,'SZ'} == c-1);
	for s = 1:size(a,3)
        for ass = 1:N.IC
            e.comp(ass, s, c) = HShannon_kNN_k_estimation(a(ass,:,s), co);
        end
	end
end
clear a s ass c

% Joint entropies
e.joint = squeeze(sum(e.comp,1));


%% Test for group-level changes

% joint entropy
isnorm = ~(jbtest(e.joint(:,1))) || jbtest(squeeze(e.joint(:,1)));
if isnorm
    [h_joint.t, p_joint.t] = ttest2(squeeze(e.joint(:,1)), squeeze(e.joint(:,2)));
    if h_joint.t
        disp(strjoin(["Student's two-tailed t-test shows significant difference between patient and control joint entropy (p = .", num2str(p_joint.t), ")."], " "));
    end
end
[h_joint.ks, p_joint.ks] = kstest2(squeeze(e.joint(:,1)), squeeze(e.joint(:,2)));
if h_joint.ks
    disp(strjoin(["Kolmogorov-Smirnov two-tailed test shows significant difference between patient and control joint entropy (p = .", num2str(p_joint.ks), ")."], " "));
end
p_joint.perm = permutationTest(squeeze(e.joint(:,1)), squeeze(e.joint(:,2)), 10000, 'exact',0, 'sidedness','both');
h_joint.perm = (p_joint.perm < 0.05);
if h_joint.perm
    disp(strjoin(["Difference-of-means permutation test shows significant difference between patient and control joint entropy (p = .", num2str(p_joint.perm), ")."], " "));
end

% subject entropy
p.ks = nan(N.IC,1);
p.perm = p.ks;
for r = 1:N.IC
    [~, p.ks(r,:)] = kstest2(e.comp(r,:,1), e.comp(r,:,2));
    p.perm(r,:) = permutationTest(e.comp(r,:,1), e.comp(r,:,2), 10000, 'exact',0, 'sidedness','both');
end
clear r

% multiple-comparison correction
FDR = nan(N.IC,2);
Bonferroni = FDR; Sidak = FDR;
[FDR(:,1), Bonferroni(:,1), Sidak(:,1)] = mCompCorr(N.IC, p.ks, 0.05);
[FDR(:,2), Bonferroni(:,2), Sidak(:,2)] = mCompCorr(N.IC, p.perm, 0.05);
Bonferroni = array2table(Bonferroni, 'VariableNames',["Kolmogorov-Smirnov", "Permutation"]);
Sidak = array2table(Sidak, 'VariableNames',["Kolmogorov-Smirnov", "Permutation"]);
FDR = array2table(FDR, 'VariableNames',["Kolmogorov-Smirnov", "Permutation"]);

% Find indices of significantly different entropies
ind.e = find(FDR{:,"Kolmogorov-Smirnov"});

% Display number of significant differences detected
disp(strjoin(["Kolmogorov-Smirnov two-tailed test detects", num2str(nnz(FDR{:,"Kolmogorov-Smirnov"})), "significant ICs"], " "));
disp(strjoin(["Difference-of-means permutation test detects", num2str(nnz(FDR{:,"Permutation"})), "significant ICs"], " "));

% Compile tables: entropy means, standard deviations
e.jointmean = mean(e.joint,'omitnan'); e.jointmean = array2table(e.jointmean, 'VariableNames',I.Properties.VariableNames);
e.compmean = squeeze(mean(e.comp,2,'omitnan')); e.compmean = array2table(e.compmean, 'VariableNames',I.Properties.VariableNames, 'RowNames',string(1:size(e.comp,1)));
e.jointstd = std(e.joint,0,'omitnan'); e.jointstd = array2table(e.jointstd, 'VariableNames',I.Properties.VariableNames);
e.compstd = squeeze(std(e.comp,0,2,'omitnan')); e.compstd = array2table(e.compstd, 'VariableNames',I.Properties.VariableNames, 'RowNames',string(1:size(e.comp,1)));
e.sigcompmean = e.compmean(FDR{:,"Kolmogorov-Smirnov"},:);
e.sigcompstd = e.compstd(FDR{:,"Kolmogorov-Smirnov"},:);


%% Visualise components with group-level changes

% Visualise joint entropy of both conditions
F(N.fig) = figure; N.fig = N.fig + 1; hold on;
boxplot(e.joint, I.Properties.VariableNames, 'Notch','on');
ylim([min(e.joint,[],'all','omitnan')-10, max(e.joint,[],'all')+10]);
title("Joint Entropy", 'FontSize',16); ylabel("Joint Entropy");
if h_joint.ks
    axes = gca;
    plot(1:2, [max(e.joint,[],'all','omitnan')+3, max(e.joint,[],'all','omitnan')+3], 'k-', 'LineWidth',1);
    plot(1.5, max(e.joint,[],'all','omitnan')+5, 'k*', 'MarkerSize',10);
end

% Plot entropies of significantly differing ICs
F(N.fig) = figure; N.fig = N.fig + 1;
for j = 1:numel(ind.e)
	subplot(ceil(numel(ind.e)/3), 3, j); hold on
	boxplot(squeeze(e.comp(ind.e(j),:,:)), I.Properties.VariableNames, 'Notch','on');
	ylim([min(e.comp(ind.e(j),:,:),[],'all')-2, max(e.comp(ind.e(j),:,:),[],'all','omitnan')+2]);
	title(strjoin(["Component", num2str(ind.e(j))], " "), 'FontSize',16);
	ylabel("Shannon Entropy");
	plot(1:2, [max(e.comp(ind.e(j),:,:),[],'all','omitnan')+0.5, max(e.comp(ind.e(j),:,:),[],'all','omitnan')+0.5], 'k-', 'LineWidth',1);
	plot(1.5, max(e.comp(ind.e(j),:,:),[],'all','omitnan')+1, 'k*', 'MarkerSize',10);
end

% Plot connectivity maps of significantly differing ICs
F(N.fig) = figure; N.fig = N.fig + 1;
fncmat = nan(N.ROI, N.ROI, numel(ind.e));
for j = 1:numel(ind.e)
	ax(j,1) = subplot(ceil(numel(ind.e)/3), 3, j); hold on
    fncmat(:,:,j) = icatb_vec2mat(sources(ind.e(j),:)');  % fncmat(:,:,j) = icatb_vec2mat(W(ind.e(j),:)');
	imagesc(fncmat(:,:,j)); colormap jet;
    clim([-max(abs(sources),[],'all'), max(abs(sources),[],'all')]); colorbar;
    set(ax(j,1), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;
    pbaspect([1 1 1]);
	title(strjoin(["Connectivity Matrix of Component", num2str(ind.e(j))], " "), 'FontSize',16);
	ylabel("Spatial Networks");
    xlabel("Spatial Networks");
end
clear j c

% Plot time courses of significantly differing ICs
ts_ic = nan(nnz(FDR{:,"Kolmogorov-Smirnov"}), size(activity,2), size(e.comp,2), 2);
for c = 1:2
	ts_ic(:,:, 1:nnz(I{:,'SZ'} == c-1), c) = activity(FDR{:,"Kolmogorov-Smirnov"},:,I{:,'SZ'} == c-1);
end
ts_ic = reshape(ts_ic, [nnz(FDR{:,"Kolmogorov-Smirnov"}), N.TR*max(N.subjects), 2]);
F(N.fig) = figure; N.fig = N.fig+1;
for j = 1:numel(ind.e)
	ax(j,2) = subplot(ceil(numel(ind.e)/3), 3, j); hold on
    plot(1:size(ts_ic,2), squeeze(ts_ic(j,:,:)));
    title(strjoin(["Time Course of Component", num2str(ind.e(j))], " "), 'FontSize',16);
    set(ax(j,2), {'XTick', 'XTickLabels', 'XLim'}, {0:N.TR:max(N.subjects)*N.TR, [], [0 max(N.subjects)*N.TR]}); hold on;
    legend(I.Properties.VariableNames);
    ylabel("Activity");
    xlabel("Subjects");
end
clear j c


%% Correlate entropy score against clinical variables

% Collate matrix (table) of joint entropy scores & 
edata = nan(nnz(isfinite(e.joint)), size(e.comp,1)+1);
edata(~I{:,'SZ'},1) = e.joint(1:nnz(~I{:,'SZ'}), 1);
edata(I{:,'SZ'},1) = e.joint(1:nnz(I{:,'SZ'}), 2);
for i = 2:size(e.comp,1)+1
    edata(~I{:,'SZ'},i) = e.comp(i-1, 1:nnz(~I{:,'SZ'}), 1);
    edata(I{:,'SZ'},i) = e.comp(i-1, 1:nnz(I{:,'SZ'}), 2);
end
edata = array2table(edata, 'VariableNames',["joint"; strcat(repmat("Component ",[N.IC 1]), string(1:size(e.comp,1))')], 'RowNames',analysis_data.Properties.RowNames);
clear i

% numerize clinical information
adata = analysis_data{:,:};
adata(:,2) = single(strcmpi(analysis_data{:,2}, "F"));
adata(:,3) = single(strcmpi(analysis_data{:,3}, "SZ"));
adata = double(adata);
% adata = array2table(adata, 'VariableNames',analysis_data.Properties.VariableNames, 'RowNames',analysis_data.Properties.RowNames);

% Compute correlation matrix with joint entropy
[rho, p.corr] = corr(adata, edata{:,:}, 'Rows','pairwise');

% Correct p-values for multiple comparison (family-wise error)
p_bin = zeros(numel(p.corr), 1);
for s = 1:length(p_bin)
	p_bin(sort(FDR_benjHoch(reshape(p.corr, [numel(p.corr),1]), 0.05)), 1) = 1;
end; clear s
p_bin = logical(squeeze(p_bin));
p_bin = reshape(p_bin, size(p.corr));
[y,x] = find(p_bin);


%% Plot correlation values & significance

% Plot correlation values (matrix)
F(N.fig) = figure; N.fig = N.fig+1;
l = numel(ind.e);
ax(l+1,1) = subplot(2,1,1);
imagesc(rho); colormap jet;
clim([-max(abs(rho),[],'all'), max(abs(rho),[],'all')]); colorbar; hold on;
xticks(1:3:size(edata,2));
xticklabels(edata.Properties.VariableNames(1:3:size(edata,2)));
yticks(1:size(adata,2));
yticklabels(analysis_data.Properties.VariableNames);
title("Correlation of Entropy to Clinical Variables", 'FontSize',18);

% Plot correlation values (bar chart)
ax(l+1,2) = subplot(2,1,2);
b = bar3(rho); hold on;
xticks(1:3:size(edata,2));
yticks(1:3:size(adata,2));
xticklabels(num2str(0:size(edata,2)));
yticklabels(num2str(0:size(edata,2)));
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
title("Correlation Coefficients", 'FontSize',18);
xlabel("Component Numbers");
clear k

% Add significance markers to correlation bar chart
plot3(ax(l+1,2), x,y, ones(numel(x),1)*1.5*max(abs(rho),[],'all'), '*k');

% Plot p-values
F(N.fig) = figure; N.fig = N.fig+1;
ax(l+2,1) = subplot(2,1,1);
imagesc(p.corr); colormap(ax(l+2,1), jet);
clim([-max(abs(p.corr),[],'all'), max(abs(p.corr),[],'all')]); colorbar; hold on;
xticks(1:3:size(edata,2));
xticklabels(edata.Properties.VariableNames(1:3:size(edata,2)));
yticks(1:size(adata,2));
yticklabels(analysis_data.Properties.VariableNames);
title("Correlation p-values of Entropy to Clinical Values", 'FontSize',18);
xlabel("Component Numbers");

% Plot binarized p-values
ax(l+2,2) = subplot(2,1,2);
imagesc(p_bin); colormap(ax(l+2,2), gray); hold on
xticks(1:3:size(edata,2));
xticklabels(edata.Properties.VariableNames(1:3:size(edata,2)));
yticks(1:size(adata,2));
yticklabels(analysis_data.Properties.VariableNames);
title("Significant Correlations of Entropy to Clinical Values", 'FontSize',18);
xlabel("Component Numbers");
% Ax = gca;
% xt = Ax.XTick';
% xtv = compose('%.0f',xt)';
% xt = xt(ind.e);
% xtv = xtv(ind.e);
% text(xt,zeros(size(xt)), xtv, 'Color','r', 'Horiz','center', 'Vert','top')

% Get signficant correlation coefficients & variables
corrvals = nan(numel(find(p_bin)),4);
[corrvals(:,1), corrvals(:,2)] = find(p_bin);
corrvals(:,3) = rho(p_bin);
corrvals(:,2) = corrvals(:,2)-1;
corrvals(:,4) = p.corr(p_bin);
corrvals = num2cell(corrvals);
corrvals(:,1) = analysis_data.Properties.VariableNames(cell2mat(corrvals(:,1)));
corrvals = cell2table(corrvals, 'VariableNames',["Variable", "Component", "rho", "p"]);


%% Save results & figure(s)

% Save figures
savefig(F, fullfile(path{4}, fileName), 'compact');
for i = 1:numel(F)
    saveas(F(i), fullfile(path{4}, "Figures", strjoin([fileName, num2str(i)], '_')), 'svg');
end
clear i F

% Save files
N.fig = N.fig - 1;
save(fullfile(path{4}, fileName));