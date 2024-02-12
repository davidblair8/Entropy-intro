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


%% Load data and convert table variables (should be placed in separate script)

% Load formatted dFNC data
load(fullfile(path{3}, 'FBIRN_DFNC_table.mat'));

% Sort data by site
[~, I(:,1)] = sort(analysis_data{:,"Site"});
analysis_data = analysis_data(I,:);
analysis_ID = analysis_ID(I,:);
analysis_SCORE = analysis_SCORE(I,:);
DFNC_FBIRN = DFNC_FBIRN(I,:);
clear I

% Set diagnosis, gender labels
labels.diagnosis = ["SZ";"HC"];
labels.gender = ["M";"F"];
labels.data = ["Diagnosis", "Gender"];
labels.I = ["Site", "Diagnosis"];

% Identify table variables
i(1,:) = contains(analysis_data.Properties.VariableNames, 'diagnosis');
i(2,:) = contains(analysis_data.Properties.VariableNames, 'gender');

% Replace table variables
analysis_data{:,i(1,:)} = labels.diagnosis(analysis_data{:,i(1,:)});
analysis_data.Properties.VariableNames(i(1,:)) = "Diagnosis";
analysis_data{:,i(2,:)} = labels.diagnosis(analysis_data{:,i(2,:)});
analysis_data.Properties.VariableNames(i(2,:)) = "Gender";

% Set figure counter
N.fig = 1;

% Index time
N.TR = size(DFNC_FBIRN{1},1);

% Set diagnosis, gender labels
labels.diagnosis = ["SZ","HC"];
labels.gender = ["M","F"];
labels.I = ["Site", "Diagnosis"];

% Locate & isolate site indices
ind.site = unique(analysis_data{:,'Site'});


%% Extract and index time series 



%% Remove site mean from time series


% get sitewise dFNC values
A = cell(numel(ind.site),1);
ts.subj = cell(numel(ind.site), numel(unique(analysis_data{:,3})));
ts.site = ts.subj;
c = nan(numel(ind.site),1);
for k = 1:numel(ind.site)
    A{k} = DFNC_FBIRN(analysis_data{:,'Site'}==ind.site(k));
    % I(analysis_data{:,'Site'}==ind.site(k), 2) = ind.site(k);   % record relation of site to index
    for j = 1:numel(unique(analysis_data{:,3}))
        ts.subj{k,j} = DFNC_FBIRN(analysis_data{:,'Site'}==ind.site(k) & analysis_data{:,3}==j); % organise time series by site, group
        I.subj{k,j} = analysis_data(analysis_data{:,'Site'}==ind.site(k) & analysis_data{:,3}==j, :);
        ts.site{k,j} = cell2mat(ts.subj{k,j});
        I.site{k,j} = repmat([ind.site(k), j], [numel(ts.subj{k,j})*N.TR,1]);
    end
    A{k} = cell2mat(A{k});
    c(k) = max(abs(cat(1,ts.site{k,:})), [], 'all');
end
c = max(c,[],'all');

% plot mean FNC matrix for each site
F(N.fig) = figure; N.fig = N.fig + 1;
F(N.fig-1).Position = get(0,'screensize');
for k = 1:numel(ind.site)
    ax(1,k) = subplot(3, numel(ind.site), k);
    imagesc(icatb_vec2mat(mean(A{k}))); colormap jet; clim([-c c]);
    pbaspect([1 1 1]); colorbar;
    ylabel('Functional Networks'); xlabel('Functional Networks');
    title(strjoin(["Mean FNC of Site", num2str(k)]));
end

% subtract mean site FNC matrix from time-resolved dFNC values
for k = 1:numel(ind.site)
    for j = 1:numel(unique(analysis_data{:,3}))
        ts.site{k,j} = ts.site{k,j} - mean(A{k});
        for q = 1:numel(ts.subj{k,j})
            ts.subj{k,j}{q} = ts.subj{k,j}{q} - mean(A{k});
        end
    end
    A{k} = A{k} - mean(A{k});
end

% plot mean group FNC after removing site means
for k = 1:numel(ind.site)
    for j = 1:numel(unique(analysis_data{:,3}))
        ax(j+1,1) = subplot(3, numel(ind.site), k+(numel(ind.site)*j));
        imagesc(icatb_vec2mat(mean(ts.site{k,j}))); colormap jet; clim([-c c]);
        pbaspect([1 1 1]); colorbar;
        ylabel('Functional Networks'); xlabel('Functional Networks');
        title([strjoin(["Mean FNC of Group", num2str(j), ", Site", num2str(k)]), strjoin(["- Mean of Site", num2str(k)])]);
    end
end
clear k j q c

% Confirm that dFNC ordering maintained
% assert(all(I(:,1) == I(:,2)), "dFNC ordering has been lost!")

% Set site index
S = I(:,1);
s = repmat(S', [N.TR 1]);
s = reshape(s, [],1);

% Concatenate dFNC time series (space = rows, time = columns)
ts.all = cell(1,numel(unique(analysis_data{:,3})));
for j = 1:numel(unique(analysis_data{:,3}))
    ts.all{j} = cell2mat(ts.site(:,j));
    I.all{j} = cell2mat(I.site(:,j));
end
ts.concat = cat(1, ts.all{:});
I.concat = cat(1, I.all{:});
A = cell2mat(A);

% Confirm that dFNC ordering maintained
assert(all(ts.concat == A, 'all'), "dFNC ordering has been lost!")

% Convert I to table
for j = 1:numel(unique(analysis_data{:,3}))
    I.all{j} = array2table(I.all{j}, 'VariableNames',labels.I);
    for k = 1:numel(ind.site)
        I.site{k,j} = array2table(I.site{k,j}, 'VariableNames',labels.I);
    end
end


%% index diagnosis, gender

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

% Index subjects
N.TR = size(DFNC_FBIRN{1},1);
N.subjects = sum(I{:,:},1);

% Set number of neighbors to search for in KNN
co = HShannon_kNN_k_initialization(1);

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


%% Set region labels

% Load functional network labels
labels.FNC = readtable(fullfile(path{3}, 'NeuroMark_FNC_labels.xlsx')); % NeuroMark functional network labels & locations

% Remove borders between functional domains
[r,~] = find(strcmpi(labels.FNC{:,:}, ""));   % find rows which separate functional domains
r = unique(r);
coords = array2table(cellfun(@str2num, labels.FNC{:,2:end}, 'UniformOutput', false), 'RowNames',labels.FNC{:,1}, 'VariableNames',labels.FNC.Properties.VariableNames(2:end));
coords(r,:) = [];

% Set number of ROIs
N.ROI = size(coords,1);

% % Set number of ROIs
% [r,~] = find(strcmpi(labels{:,:}, ""));   % find rows which separate functional domains
% r = numel(unique(r));
% N.ROI = size(coords,1) - r;
% clear r


%% Set number of ICs

% % Find maximum number of components
% disp("Identifying number independent components from Marcenko-Pasteur distribution.");
% N.IC = NumberofIC(ts);

% % Evaluate IC counts vs. captured variance
% d = factor(N.IC);
% N.IC = d(1)*d(2):d(1)*d(2):N.IC; clear d
% [ev, F(N.fig)] = evaluateICnumbers(N, ts);
% N.fig = N.fig + 1;

% Set number of ICs
N.IC = 22;


%% Define filename based on parameters

% Get number of ICs
fileName = strcat(num2str(N.IC), "ICs-SiteMeanFNC");

% Set iteration number
fList = dir(fullfile(path{4}, strcat(strjoin([fileName, "iteration"], '-'), '*.mat')));	% Get file list
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
sFNC{2} = nan(numel(I.Properties.VariableNames), size(ts.all,1));
for k = 1:numel(I.Properties.VariableNames)
    sFNC{2}(k,:) = mean(sFNC{1}(I{:,'SZ'}==(k-1),:));
end

% Convert to square matrices
sFNC{1} = icatb_vec2mat(sFNC{1});
sFNC{2} = icatb_vec2mat(sFNC{2});
N.ROI = size(sFNC{1},2);

% % Check if need intercept term
% if intercept == true
% 	I = horzcat(ones(size(I,1),1), I);
% end
% 
% % Allocate contrast matrices
% if size(comps,1) > 1
%     cont = zeros((size(comps,1)+size(I,2))*2, size(I,2));
%     strcont = cell((size(comps,1)+size(I,2))*2, 1);
% 
%     % Set all-way contrasts
%     c = 2*(size(comps,1)+1:size(comps,1)+N.conditions);
%     cont(c-1, :) = -ones(size(I,2), size(I,2)) + diag(2*ones(N.conditions,1));
%     cont(c, :) = ones(size(I,2), size(I,2)) - diag(2*ones(N.conditions,1));
%     for c = 2*(size(comps,1)+1:size(comps,1)+N.conditions)
%         strcont{c-1} = strjoin([groups((c-2*size(comps,1))/2), '> ALL']);
%         strcont{c} = strjoin([groups((c-2*size(comps,1))/2), '< ALL']);
%     end
% else
%     cont = zeros(size(comps,1)*2, size(I,2));
%     strcont = cell(size(comps,1)*2, 1);
% end
% 
% % Set pairwise contrasts
% for c = 2*(1:size(comps,1))
%     cont(c-1, comps(c/2,:)) = [1 -1];
%     cont(c, comps(c/2,:)) = [-1 1];
%     strcont{c-1} = strjoin([groups(comps(c/2,1)), '>', groups(comps(c/2,2))]);
%     strcont{c} = strjoin([groups(comps(c/2,1)), '<', groups(comps(c/2,2))]);
% end
% clear c
% 
% % Run NBS analysis
% sFNC{1} = permute(sFNC{1}, [2 3 1]);
% [nbs, STATS, GLM, storarray] = runNBS(sFNC{1}, cont, I, N, tstat);
% 
% % Visualise group-level sFNC
% F(N.fig) = figure; N.fig = N.fig + 1;
% for k = 1:numel(I.Properties.VariableNames)
% 	ax(1,k) = subplot(1,numel(I.Properties.VariableNames),k);
%     imagesc(squeeze(sFNC{2}(k,:,:))); colormap jet;
%     clim([-max(abs(sFNC{2}),[],'all'), max(abs(sFNC{2}),[],'all')]); colorbar;
%     set(ax(1,k), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;
%     pbaspect([1 1 1]);
% 	title(strjoin(["Mean sFNC of", I.Properties.VariableNames{k}], " "), 'FontSize',16);
% end
% clear k


%% Isolate components & activity from dFNC

% extract components & activities from dFNC
% [weights,sphere,activations] = icatb_runica(ts, 'ncomps',N.IC);
% weights = weights';
% [activities, memberships, W] = fastica(ts, 'numOfIC', N.IC, 'verbose','off'); % time series input should be in format space x time

% extract components & activities from dFNC
[whitesig, dewhiteM] = icatb_calculate_pca(A', N.IC);
[~, ~, A, sources] = icatb_icaAlgorithm(1, whitesig');
A = dewhiteM*A;
% W = pinv(A);

% Split activity, time series by subject
activity = reshape(A, [N.IC, N.TR, sum(N.subjects)]);
S = repmat(S', [N.TR,1]);

% reshape activities to allow easy plotting
ts.ic = nan(N.IC, N.TR, max(N.subjects), 2);
s = nan(N.TR, max(N.subjects), 2);
for c = 1:2
	ts.ic(:,:, 1:nnz(I{:,'HC'} == c-1), c) = activity(:,:,I{:,'HC'} == c-1);
    s(:, 1:nnz(I{:,'HC'} == c-1), c) = S(:,I{:,'HC'} == c-1);
end
ts.ic = reshape(ts.ic, [N.IC, N.TR*max(N.subjects), 2]);
s = reshape(s, [N.TR*max(N.subjects), 2]);

% Select sample connectivity maps
N.samplot = min([5, N.IC]);
C = nchoosek(1:N.IC, N.samplot);
C = C(ceil(rand*size(C,1)),:);

% Plot sample connectivity maps
F(N.fig) = figure; N.fig = N.fig + 1;
F(N.fig-1).Position = get(0,'screensize');
fncmat = nan(N.ROI, N.ROI, N.samplot);
for j = 1:N.samplot
    % Plot connectivity maps
	ax(j,1) = subplot(3, N.samplot, j); hold on
    fncmat(:,:,j) = icatb_vec2mat(sources(C(j),:)');
	imagesc(fncmat(:,:,j)); colormap jet;
    clim([-max(abs(sources),[],'all'), max(abs(sources),[],'all')]); colorbar;
    set(ax(j,1), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;
    pbaspect([1 1 1]);
	title(strjoin(["FNC of Component", num2str(C(j))], " "), 'FontSize',16);
	ylabel("Spatial Networks");
    xlabel("Spatial Networks");

    % Plot time courses
    ax(j,2) = subplot(3, N.samplot, j+N.samplot); hold on
    plot(1:N.TR*max(N.subjects), squeeze(ts.ic(C(j),:,:)));
    title(strjoin(["Time Course of Component", num2str(C(j))], " "), 'FontSize',16);
    set(ax(j,2), {'XTick', 'XTickLabels', 'XLim'}, {0:N.TR:max(N.subjects)*N.TR, [], [0 max(N.subjects)*N.TR]}); hold on;
    legend(I.Properties.VariableNames);
    ylabel("Activity");
    xlabel("Subjects");

    % Plot site of recording
    ax(j,3) = subplot(3, N.samplot, j+(2*N.samplot));
    plot(s); hold on
    title(strjoin("Site of Recording"), 'FontSize',16);
    set(ax(j,3), {'XTick', 'XTickLabels', 'XLim'}, {0:N.TR:max(N.subjects)*N.TR, [], [0 max(N.subjects)*N.TR]}); hold on;
    legend(I.Properties.VariableNames);
    ylabel("Site");
    xlabel("Subjects");
end
clear j c s


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
clear isnorm

% subject entropy
p.ks = nan(N.IC,1);
p.perm = p.ks;
p.t = p.ks;
for r = 1:N.IC
    [~, p.ks(r,:)] = kstest2(e.comp(r,:,1), e.comp(r,:,2));
    p.perm(r,:) = permutationTest(e.comp(r,:,1), e.comp(r,:,2), 10000, 'exact',0, 'sidedness','both');
    isnorm = ~(jbtest(e.comp(r,:,1))) || jbtest(squeeze(e.comp(r,:,2)));
    if isnorm
        [~, p.t(r,:)] = ttest2(squeeze(e.comp(r,:,1)), squeeze(e.comp(r,:,2)));
    end
end
clear r isnorm

% multiple-comparison correction
FDR = nan(N.IC,2);
Bonferroni = FDR; Sidak = FDR;
[FDR(:,2), Bonferroni(:,2), Sidak(:,1)] = mCompCorr(N.IC, p.t, 0.05);
[FDR(:,2), Bonferroni(:,2), Sidak(:,2)] = mCompCorr(N.IC, p.ks, 0.05);
[FDR(:,3), Bonferroni(:,3), Sidak(:,3)] = mCompCorr(N.IC, p.perm, 0.05);

% convert NaNs to false
FDR(isnan(FDR)) = false;
Bonferroni(isnan(Bonferroni)) = false;
Sidak(isnan(Sidak)) = false;

% convert to logical arrays
Bonferroni = array2table(logical(Bonferroni), 'VariableNames',["Student's t", "Kolmogorov-Smirnov", "Permutation"]);
Sidak = array2table(logical(Sidak), 'VariableNames',["Student's t", "Kolmogorov-Smirnov", "Permutation"]);
FDR = array2table(logical(FDR), 'VariableNames',["Student's t", "Kolmogorov-Smirnov", "Permutation"]);

% Find indices of significantly different entropies
if nnz(FDR{:,"Student's t"}) > 0
    ind.e = find(FDR{:,"Student's t"});
elseif nnz(FDR{:,"Permutation"}) > 0
    ind.e = find(FDR{:,"Permutation"});
else
    ind.e = find(FDR{:,"Kolmogorov-Smirnov"});
end

% Display number of significant differences detected
disp(strjoin(["Student's two-sample t-test detects", num2str(nnz(FDR{:,"Student's t"})), "significant ICs"], " "));
disp(strjoin(["Kolmogorov-Smirnov two-tailed test detects", num2str(nnz(FDR{:,"Kolmogorov-Smirnov"})), "significant ICs"], " "));
disp(strjoin(["Difference-of-means permutation test detects", num2str(nnz(FDR{:,"Permutation"})), "significant ICs"], " "));

% Compile tables: entropy means, standard deviations
e.jointmean = mean(e.joint,'omitnan'); e.jointmean = array2table(e.jointmean, 'VariableNames',I.Properties.VariableNames);
e.compmean = squeeze(mean(e.comp,2,'omitnan')); e.compmean = array2table(e.compmean, 'VariableNames',I.Properties.VariableNames, 'RowNames',string(1:size(e.comp,1)));
e.jointstd = std(e.joint,0,'omitnan'); e.jointstd = array2table(e.jointstd, 'VariableNames',I.Properties.VariableNames);
e.compstd = squeeze(std(e.comp,0,2,'omitnan')); e.compstd = array2table(e.compstd, 'VariableNames',I.Properties.VariableNames, 'RowNames',string(1:size(e.comp,1)));
if nnz(FDR{:,"Kolmogorov-Smirnov"}) > 0
    e.sigcompmean = e.compmean(FDR{:,"Kolmogorov-Smirnov"},:);
    e.sigcompstd = e.compstd(FDR{:,"Kolmogorov-Smirnov"},:);
end


%% Visualise components with group-level changes

% Visualise joint entropy of both conditions
F(N.fig) = figure; N.fig = N.fig + 1;
F(N.fig-1).Position = get(0,'screensize'); hold on;
boxplot(e.joint, I.Properties.VariableNames, 'Notch','on');
ylim([min(e.joint,[],'all','omitnan')-10, max(e.joint,[],'all')+10]);
title("Joint Entropy", 'FontSize',16); ylabel("Joint Entropy");
if isfield(h_joint, "t") && h_joint.t
    axes = gca;
    plot(1:2, [max(e.joint,[],'all','omitnan')+3, max(e.joint,[],'all','omitnan')+3], 'k-', 'LineWidth',1);
    plot(1.5, max(e.joint,[],'all','omitnan')+5, 'k*', 'MarkerSize',10);
elseif h_joint.ks
    axes = gca;
    plot(1:2, [max(e.joint,[],'all','omitnan')+3, max(e.joint,[],'all','omitnan')+3], 'k-', 'LineWidth',1);
    plot(1.5, max(e.joint,[],'all','omitnan')+5, 'k*', 'MarkerSize',10);
end

% Plot entropies of significantly differing ICs
F(N.fig) = figure; N.fig = N.fig + 1;
F(N.fig-1).Position = get(0,'screensize');
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
F(N.fig-1).Position = get(0,'screensize');
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
if nnz(FDR{:,"Kolmogorov-Smirnov"}) > 0
    ts.ic = nan(nnz(FDR{:,"Kolmogorov-Smirnov"}), size(activity,2), size(e.comp,2), 2);
    for c = 1:2
	    ts.ic(:,:, 1:nnz(I{:,'SZ'} == c-1), c) = activity(FDR{:,"Kolmogorov-Smirnov"},:,I{:,'SZ'} == c-1);
    end
    ts.ic = reshape(ts.ic, [nnz(FDR{:,"Kolmogorov-Smirnov"}), N.TR*max(N.subjects), 2]);
end
F(N.fig) = figure; N.fig = N.fig+1;
F(N.fig-1).Position = get(0,'screensize');
for j = 1:numel(ind.e)
	ax(j,2) = subplot(ceil(numel(ind.e)/3), 3, j); hold on
    plot(1:size(ts.ic,2), squeeze(ts.ic(j,:,:)));
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
F(N.fig-1).Position = get(0,'screensize');
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
F(N.fig-1).Position = get(0,'screensize');
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
for c = 1:numel(F)
    saveas(F(c), fullfile(path{4}, "Figures", strjoin([fileName, num2str(c)], '_')), 'svg');
end
clear c F

% Save files
N.fig = N.fig - 1;
save(fullfile(path{4}, fileName));