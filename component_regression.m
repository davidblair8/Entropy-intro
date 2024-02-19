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
path{4,1} = fullfile(path{2}, 'Results','Regression');

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


%% Load and sort data

% Load formatted dFNC data
load(fullfile(path{3}, 'FBIRN_DFNC_table.mat'));
load(fullfile(path{3}, 'head_motion_meanFD.mat'));

% Confirm that IDs, data are properly indexed
assert(all(str2double(string(cell2mat(analysis_ID))) == str2double(analysis_data.Properties.RowNames)), "Data labels are not properly ordered!");
clear analysis_ID analysis_SCORE
assert(all(strcmpi(string(FILE_ID), string(analysis_data.Properties.VariableNames))), "Clinical variables are not properly ordered!");
clear FILE_ID

% add head motion to data array
analysis_data = [table(head_motion_meanFD, 'VariableNames',"Mean Head Motion"), analysis_data];
clear head_motion_meanFD

% % Sort data by site
% [~, I(:,1)] = sort(analysis_data{:,"Site"});
% analysis_data = analysis_data(I,:);
% DFNC_FBIRN = DFNC_FBIRN(I,:);
% clear I

% Set diagnosis, gender labels
labels.diagnosis = ["SZ";"HC"];
labels.gender = ["M";"F"];
labels.data = ["Diagnosis", "Gender"];


%% Index counters

% Set figure counter
N.fig = 1;

% Index time
N.TR = size(DFNC_FBIRN{1},1);

% identify diagnosis column
i(1,:) = contains(analysis_data.Properties.VariableNames, 'diagnosis');

% Index conditions
d = unique(analysis_data{:,i});
N.conditions = numel(d);

% Index subjects
I = cell(1,N.conditions);
for j = 1:N.conditions
    I{j} = nnz(analysis_data{:,i} == d(j));
end
N.subjects = cell2table(I, 'VariableNames',labels.diagnosis);

% Locate & isolate site indices
ind.site = unique(analysis_data{:,'Site'});
clear I d j i


%% Set analysis parameters

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

% Determine whether to use intercept term in NBS
intercept = false;

% Set parameter sweeps
tstat = 3:0.5:5;


%% Convert table variables (should be placed in separate script)

% Replace numeric missing data code with NaN
analysis_data{:,:}(analysis_data{:,:} == -9999) = NaN;

% Identify table variables to change
i(1,:) = contains(analysis_data.Properties.VariableNames, 'diagnosis');
i(2,:) = contains(analysis_data.Properties.VariableNames, 'gender');

% generate string arrays
groups = labels.diagnosis(analysis_data{:,i(1,:)});
gender = labels.gender(analysis_data{:,i(2,:)});

% Convert variable type
analysis_data = convertvars(analysis_data, ["diagnosis(1:sz; 2:hc)","gender(1:male; 2:female)"], "string");

% Replace table numeric indices with strings
analysis_data{:,i(1,:)} = groups;
analysis_data{:,i(2,:)} = gender;
clear groups gender i

% Rename table variables
analysis_data = renamevars(analysis_data, ["age" "diagnosis(1:sz; 2:hc)" "gender(1:male; 2:female)"], ["Age" "Diagnosis" "Gender"]);


%% Set region labels & maps

% Load functional network labels
labels.FNC = readtable(fullfile(path{3}, 'NeuroMark_FNC_labels.xlsx')); % NeuroMark functional network labels & locations
labels.FNC = renamevars(labels.FNC, "SelectedComponentsAsRegionsOfInterest", "Functional Networks");

% Remove borders between functional domains
[r,~] = find(strcmpi(labels.FNC{:,:}, ""));   % find rows which separate functional domains
r = unique(r);
labels.ROI = array2table(cellfun(@str2num, labels.FNC{:,2:end}, 'UniformOutput', false), 'RowNames',labels.FNC{:,1}, 'VariableNames',labels.FNC.Properties.VariableNames(2:end));
labels.ROI(r,:) = [];

% Set functional domain labels
labels.FDs = labels.FNC(r,"Functional Networks");
labels.FDs = renamevars(labels.FDs, "Functional Networks", "Functional Domains");

% Set number of ROIs, FDs
N.ROI = size(labels.ROI,1);
N.FD = size(labels.FDs,1);

% Establish FN-level map of FDs
r = [r; size(labels.FNC,1)+1];
labels.FND = labels.FNC;
for i = 2:numel(r)
    labels.FND{r(i-1)+1:r(i)-1,"Functional Networks"} = repmat(labels.FDs{i-1,"Functional Domains"}, [r(i)-1-r(i-1) ,1]);
end
labels.FND(r(1:end-1),:) = [];

% Establish FN-level map of FDs
ind.FND = zeros(N.ROI, N.FD);
for d = 1:N.FD
    ind.FND(:,d) = strcmpi(labels.FND{:,"Functional Networks"}, labels.FDs{d,"Functional Domains"});
end
ind.FND = array2table(ind.FND, 'RowNames',labels.ROI.Properties.RowNames, 'VariableNames',labels.FDs{:,"Functional Domains"});
clear d i r


%% Concatenate and index time series

% rename FNC data
FNC.subj = DFNC_FBIRN;
FNC.full = cell2mat(DFNC_FBIRN);

% convert variables to row form
I.subject = str2double(string(analysis_data.Properties.RowNames)');
I.diagnosis = analysis_data{:,"Diagnosis"}';
I.gender = analysis_data{:,"Gender"}';
I.site = analysis_data{:,"Site"}';
I.age = analysis_data{:,"Age"}';

% fill in data for each timepoint
f = fieldnames(I);
for j = 1:numel(f)
    I.(f{j}) = repmat(I.(f{j}), [N.TR 1]);
    I.(f{j}) = reshape(I.(f{j}), [sum(N.subjects{:,:})*N.TR 1]);    % reshape indices to column form
end
clear j f k

% Confirm that indices are in proper order
assert(all(unique(I.subject)==str2double(string(analysis_data.Properties.RowNames))), "Indices are out of order!");

% convert index to table
I = struct2table(I);


%% Set number of ICs

% % Find maximum number of components
% disp("Identifying number independent components from Marcenko-Pasteur distribution.");
% N.IC = NumberofIC(FNC.full);

% % Evaluate IC counts vs. captured variance
% d = factor(N.IC);
% N.IC = d(1)*d(2):d(1)*d(2):N.IC; clear d
% [ev, F(N.fig)] = evaluateICnumbers(N, ts);
% N.fig = N.fig + 1;

% Set number of ICs
N.IC = 9;


%% Define filename based on parameters

% Get number of ICs
fileName = fullfile(strcat(num2str(N.IC), "ICs"));

% Set iteration number
fList = dir(fullfile(path{4}, strcat(strjoin([fileName, "iteration"], '-'), '*.mat')));	% Get file list
a = false(numel(fList),1);
for n = 1:numel(fList)
    a(n) = matches("iteration", strsplit(fLis.name, '_'));
end
nIter = numel(fList)-sum(a)+1;

% Set full filename
fileName = strjoin([fileName, strcat("iteration", num2str(nIter))], '_');
clear fList nIter a
clear a k n


%% Visualise static FNC matrix with indices

% subject-level sFNC
sFNC{1} = cell2mat(cellfun(@mean, DFNC_FBIRN, 'UniformOutput', false));
sFNC{1} = icatb_vec2mat(sFNC{1});       % Convert to square matrices

% group-level sFNC
sFNC{2} = nan(N.conditions, N.ROI, N.ROI);
F(N.fig) = figure; F(N.fig).Position = get(0,'screensize');
N.fig = N.fig + 1;
lim.c = max(abs(FNC.full), [], 'all');
for c = 1:N.conditions
    % compute mean group-level FNC for condition
    sFNC{2}(c,:,:) = mean(sFNC{1}(strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c)),:,:),1,"omitmissing");

    % Visualize group sFNC
    a(1,c) = subplot(1, N.conditions, c);
    display_FNC(squeeze(sFNC{2}(c,:,:)), 0.8);    % imagesc(squeeze(sFNC{2}(c,:,:)));
    % colormap jet; clim([-lim.c lim.c]);
    pbaspect([1 1 1]); colorbar;
    ylabel('Functional Networks'); xlabel('Functional Networks');
    title(strjoin(["Mean FNC of", labels.diagnosis(c)]));
end
clear c


%% Run network-based statistic

% % Check if need NBS intercept term
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
% for k = 1:numel(labels.diagnosis)
% 	ax(1,k) = subplot(1,numel(labels.diagnosis),k);
%     imagesc(squeeze(sFNC{2}(k,:,:))); colormap jet;
%     clim([-max(abs(sFNC{2}),[],'all'), max(abs(sFNC{2}),[],'all')]); colorbar;
%     set(ax(1,k), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;
%     pbaspect([1 1 1]);
% 	title(strjoin(["Mean sFNC of", labels.diagnosis{k}], " "), 'FontSize',16);
% end
% clear k


%% Isolate components & activity from dFNC

% extract components & activities from dFNC
[whitesig, dewhiteM] = icatb_calculate_pca(FNC.full, N.IC);     % ICA runs along the first (column) dimension.
[~, ~, A, sources] = icatb_icaAlgorithm(1, whitesig');     %   for spatial ICA, orient data as space x time
A = dewhiteM*A;                                                 %   for temporal ICA, orient data as time x space
W = pinv(A);

% Set concatenated time series structure
sources = sources';

% Set colorbar, y-limits
lim.color = max(abs(W),[],'all');        % A = mixing matrix; W = unmixing matrix
lim.y = max(abs(sources),[],'all');   % ts.sources.full = component time courses

% Select sample connectivity maps
N.samplot = min([5, N.IC]);
C = nchoosek(1:N.IC, N.samplot);
C = C(ceil(rand*size(C,1)),:);


%% Calculate subject-level entropy & joint entropy

% Find component entropies
entropy = nan(sum(N.subjects{:,:}), N.IC+1);
E = strcat(repmat("Comp", [N.IC+1, 1]), string(1:N.IC+1)')';

for s = 1:sum(N.subjects{:,:})
    ts = sources(I{:,"subject"} == str2double(analysis_data.Properties.RowNames{s}), :);
    for k = 1:N.IC
        entropy(s, k) = HShannon_kNN_k_estimation(ts(:,k)', co);
    end
end
clear s ts k si

% Joint entropies
entropy(:,N.IC+1) = squeeze(sum(entropy(:,1:N.IC), 2));
E(N.IC+1) = "Joint";

% Convert entropy to table
entropy = array2table(entropy, "RowNames",analysis_data.Properties.RowNames, "VariableNames",E);


%% Regress confounder variables against entropy

% Format confounder variables
X = analysis_data(:, ["Diagnosis" "Age" "Gender"]);
X = convertvars(X, labels.data, 'double');
X{:,"Diagnosis"} = strcmpi(labels.diagnosis(1), analysis_data{:,"Diagnosis"});
X{:,"Gender"} = strcmpi(labels.gender(1), analysis_data{:,"Gender"});
X = renamevars(X, labels.data, ["Diagnosis" "Gender (1=M, 0=F)"]);

% add site information
S = zeros(sum(N.subjects{:,:}), numel(ind.site));
for s = 1:numel(ind.site)
    S(:,s) = analysis_data{:,"Site"} == ind.site(s);
end
S = array2table(S, "VariableNames",strcat(repmat("Site", [numel(ind.site) 1]), num2str(ind.site)));
X = horzcat(X,S); clear S

% fit linear model
diagnosis = table('Size',[N.IC+1, 4], 'VariableTypes',repmat("double",[1 4]));
for c = 1:N.IC+1
    x = horzcat(X, entropy(:,c));
    mdl.(x.Properties.VariableNames{end}) = fitlm(x);
    diagnosis(c,:) = mdl.(x.Properties.VariableNames{end}).Coefficients("Diagnosis",:);
end
diagnosis.Properties.RowNames = entropy.Properties.VariableNames;
diagnosis.Properties.VariableNames = mdl.(x.Properties.VariableNames{end}).Coefficients.Properties.VariableNames;
clear c s x S


%% Regress clinical variables against entropy (patients only)

% Format confounder variables
X = analysis_data(:, ["Diagnosis" "Age" "Gender"]);
X = convertvars(X, labels.data, 'double');
X{:,"Diagnosis"} = strcmpi(labels.diagnosis(1), analysis_data{:,"Diagnosis"});
X{:,"Gender"} = strcmpi(labels.gender(1), analysis_data{:,"Gender"});
X = renamevars(X, labels.data, ["Diagnosis" "Gender (1=M, 0=F)"]);

% Format site information
S = zeros(sum(N.subjects{:,:}), numel(ind.site));
for s = 1:numel(ind.site)
    S(:,s) = analysis_data{:,"Site"} == ind.site(s);
end
S = array2table(S, "VariableNames",strcat(repmat("Site", [numel(ind.site) 1]), num2str(ind.site)));
X = horzcat(X,S);

% Add clinical variables
X = horzcat(analysis_data(:,contains(analysis_data.Properties.VariableNames, "PANSS")), X); % fitlm ignores rows with missing values

% trim entropy & confounder matrices: patients only
e = entropy(strcmpi("SZ", analysis_data{:,"Diagnosis"}), :);
X = X(strcmpi("SZ", analysis_data{:,"Diagnosis"}), ~strcmpi("Diagnosis", X.Properties.VariableNames));

% trim entropy & confounder matrices: remove rows with NaN values
e = e(~(isnan(X{:,"PANSS(positive)"}) | isnan(X{:,"PANSS(negative)"})), :);
X = X(~(isnan(X{:,"PANSS(positive)"}) | isnan(X{:,"PANSS(negative)"})), :);

% fit linear model
clinical = table('Size',[2*(N.IC+1), 4], 'VariableTypes',repmat("double",[1 4]));
for c = 1:N.IC+1
    x = horzcat(X, e(:,c));
    mdl.(x.Properties.VariableNames{end}) = fitlm(x);
    clinical(2*c-1,:) = mdl.(x.Properties.VariableNames{end}).Coefficients("PANSS(positive)",:);
    clinical(2*c,:) = mdl.(x.Properties.VariableNames{end}).Coefficients("PANSS(negative)",:);
end
clinical.Properties.RowNames = join([entropy.Properties.VariableNames(floor(1:0.5:N.IC+1.5))' repmat(["PANSS(positive)";"PANSS(negative)"], [N.IC+1,1])]);
clinical.Properties.VariableNames = mdl.(x.Properties.VariableNames{end}).Coefficients.Properties.VariableNames;
clear c s x S


%% Test cognitive variables for cross-correlation

% Correlate analysis data variables against each other
datcorr = corr(analysis_data{:, ~contains(analysis_data.Properties.VariableNames, ["Diagnosis" "Gender" "Site"])}, 'rows','pairwise');
i = analysis_data.Properties.VariableNames(~contains(analysis_data.Properties.VariableNames, ["Diagnosis" "Gender" "Site"]));

% Visualize covariance matrix
sFNC{2} = nan(N.conditions, N.ROI, N.ROI);
F(N.fig) = figure; F(N.fig).Position = get(0,'screensize');
N.fig = N.fig + 1;
lim.c = max(abs(datcorr), [], 'all');
imagesc(datcorr); colormap jet; pbaspect([1 1 1]); colorbar;
clim([-max(abs(datcorr), [], 'all') max(abs(datcorr), [], 'all')]);
ylabel('Cognitive Variables'); xlabel('Cognitive Variables');
xticks(1:numel(i)); yticks(1:numel(i));
xticklabels(i); yticklabels(i);
title("Correlation Matrix of Cognitive Variables");
clear i


%% Regress cognitive variables against entropy

% Format confounder variables
% X = analysis_data(:, ~contains(analysis_data.Properties.VariableNames, ["PANSS" "Site"]));
X = analysis_data(:, ["Diagnosis" "Age" "Gender" "CMINDS_composite"]);
X = convertvars(X, "Gender", 'double');
X{:,"Gender"} = strcmpi(labels.gender(1), analysis_data{:,"Gender"});
X = renamevars(X, "Gender", "Gender (1=M, 0=F)");

% add site information
S = zeros(sum(N.subjects{:,:}), numel(ind.site));
for s = 1:numel(ind.site)
    S(:,s) = analysis_data{:,"Site"} == ind.site(s);
end
S = array2table(S, "VariableNames",strcat(repmat("Site", [numel(ind.site) 1]), num2str(ind.site)));
X = horzcat(X,S); clear S

% fit linear model
cognitive = table('Size',[N.IC+1, 4], 'VariableTypes',repmat("double",[1 4]));
for c = 1:N.IC+1
    x = horzcat(X, entropy(:,c));
    mdl.(x.Properties.VariableNames{end}) = fitlm(x);
    cognitive(c,:) = mdl.(x.Properties.VariableNames{end}).Coefficients("CMINDS_composite",:);
end
cognitive.Properties.RowNames = entropy.Properties.VariableNames;
cognitive.Properties.VariableNames = mdl.(x.Properties.VariableNames{end}).Coefficients.Properties.VariableNames;
clear c s x S


%% Test for group-level changes

% Compile entropy array: subject x condition x component
E = nan(max(N.subjects{:,:}), N.conditions, N.IC+1);
for k = 1:N.conditions
    for c = 1:N.IC
        E(1:N.subjects{:,labels.diagnosis(k)},k,c) = entropy{strcmpi(analysis_data{:,"Diagnosis"}, labels.diagnosis(k)), strcat("Comp",num2str(c))};
    end
    E(1:N.subjects{:,labels.diagnosis(k)},k, N.IC+1) = entropy{strcmpi(analysis_data{:,"Diagnosis"}, labels.diagnosis(k)), "Joint"};
end

% component entropy comparisons
p = table('Size',[N.IC+1,3], 'VariableTypes',{'double','double','double'}, 'VariableNames',["Student's t-test" "Kolmogorov-Smirnov" "Permutation"], 'RowNames',entropy.Properties.VariableNames);
h = table('Size',[N.IC+1,3], 'VariableTypes',{'logical','logical','logical'}, 'VariableNames',["Student's t-test" "Kolmogorov-Smirnov" "Permutation"], 'RowNames',entropy.Properties.VariableNames);
for c = 1:N.IC+1
    [h(c,:), p(c,:)] = sigtest(squeeze(E(:,:,c)));
end

% add regression results to p, h
p = horzcat(diagnosis(:,'pValue'), p);
p = renamevars(p, "pValue", "Regression");
t = table(p{:,"Regression"} < 0.05, 'VariableNames',"regression");
h = horzcat(t, h);
clear t

% display result of joint entropy comparisons
if h{"Joint","Student's t-test"}
    disp(strjoin(["Student's two-sample t-test shows significant difference between patient and control joint entropy (p = .", num2str(p{"Joint","Student's t-test"}), ")."], " "));
end
if h{"Joint",'Kolmogorov-Smirnov'}
    disp(strjoin(["Kolmogorov-Smirnov two-tailed test shows significant difference between patient and control joint entropy (p = .", num2str(p{"Joint","Kolmogorov-Smirnov"}), ")."], " "));
end
if h{"Joint",'Permutation'}
    disp(strjoin(["Difference-of-means permutation test shows significant difference between patient and control joint entropy (p = .", num2str(p{"Joint","Permutation"}), ")."], " "));
end

% multiple-comparison correction
multcorr.FDR = table('Size',[N.IC+1,numel(p.Properties.VariableNames)], 'VariableTypes',repmat("double",[numel(p.Properties.VariableNames) 1]), ...
    'VariableNames',p.Properties.VariableNames, 'RowNames',entropy.Properties.VariableNames);
multcorr.Bonferroni = multcorr.FDR;
multcorr.Sidak = multcorr.FDR;
for t = 1:numel(p.Properties.VariableNames)
    c = p.Properties.VariableNames(t);
    [multcorr.FDR{:,c}, multcorr.Bonferroni{:,c}, multcorr.Sidak{:,c}] = mCompCorr(N.IC+1, p{:,c}, 0.05);
end
clear c k t

% convert NaNs to false
S = string(fieldnames(multcorr));
for s = 1:numel(S)
    c = logical(multcorr.(S(s)){:,:});
    c(isnan(table2array(multcorr.(S(s))))) = false;
    multcorr.(S(s)){:,:} = c;
end
clear c k t s

% Find indices of individual components with significantly different entropies
i = contains(string(entropy.Properties.VariableNames), "Comp");
if nnz(multcorr.FDR{i,"Student's t-test"}) > 0
    ind.e = find(multcorr.FDR{i,"Student's t-test"});
elseif nnz(multcorr.FDR{i,"Permutation"}) > 0
    ind.e = find(multcorr.FDR{i,"Permutation"});
else
    ind.e = find(multcorr.FDR{i,"Kolmogorov-Smirnov"});
end

% Display number of significant differences detected
disp(strjoin(["Student's two-sample t-test detects", num2str(nnz(multcorr.FDR{i,"Student's t-test"})), "significant ICs"], " "));
disp(strjoin(["Kolmogorov-Smirnov two-tailed test detects", num2str(nnz(multcorr.FDR{i,"Kolmogorov-Smirnov"})), "significant ICs"], " "));
disp(strjoin(["Difference-of-means permutation test detects", num2str(nnz(multcorr.FDR{i,"Permutation"})), "significant ICs"], " "));
clear i

% Compile tables: entropy group mean & variances
stats = nan(2*numel(labels.diagnosis), N.IC+1);
i = [repmat("Mean", [numel(labels.diagnosis) 1]); repmat("Variance", [numel(labels.diagnosis) 1])];
for k = 1:numel(labels.diagnosis)
    stats(k,:) = mean(entropy{strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(k)),:},'omitnan');
    stats(k+numel(labels.diagnosis),:) = var(entropy{strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(k)),:},'omitnan');
    i(k) = strjoin([labels.diagnosis(k), i(k)]);
    i(k+numel(labels.diagnosis)) = strjoin([labels.diagnosis(k), i(k+numel(labels.diagnosis))]);
end
stats = array2table(stats, "RowNames",i, "VariableNames",entropy.Properties.VariableNames);


%% Visualise components with group-level changes

% Visualise joint entropy of both conditions
F(N.fig) = figure; N.fig = N.fig + 1;
F(N.fig-1).Position = get(0,'screensize'); hold on;
boxplot(squeeze(E(:,:,N.IC+1)), labels.diagnosis, 'Notch','on');
ylim([min(squeeze(E(:,:,N.IC+1)),[],'all',"omitmissing")-10, max(squeeze(E(:,:,N.IC+1)),[],'all',"omitmissing")+10]);
title("Joint Entropy", 'FontSize',16); ylabel("Joint Entropy");
if h{"Joint",'Student''s t-test'}
    axes = gca;
    plot(1:2, [max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+3, max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+3], 'k-', 'LineWidth',1);
    plot(1.5, max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+5, 'k*', 'MarkerSize',10);
elseif h{"Joint",'Kolmogorov-Smirnov'}
    axes = gca;
    plot(1:2, [max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+3, max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+3], 'k-', 'LineWidth',1);
    plot(1.5, max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+5, 'k*', 'MarkerSize',10);
end

% Plot significantly differing ICs
if isfield(ind, 'e')
    F(N.fig) = figure; N.fig = N.fig + 1;
    F(N.fig-1).Position = get(0,'screensize');
    for j = 1:numel(ind.e)
        % compile time series
        ts = nan(N.TR*max(N.subjects{:,:}), N.conditions);
        for k = 1:N.conditions
            ts(1:N.TR*N.subjects{:,labels.diagnosis(k)},k) = sources(strcmpi(I{:,"diagnosis"}, labels.diagnosis(k)), ind.e(j));
        end

        % Plot ICA source matrices
        ax(j,1) = subplot(3, numel(ind.e), j);
        display_FNC(icatb_vec2mat(W(ind.e(j),:)), 0.8);     % imagesc(icatb_vec2mat(W(ind.e(j),:)));
        % colormap jet; clim([-lim.color lim.color]);
        pbaspect([1 1 1]); colorbar; hold on
        ylabel('Functional Networks'); xlabel('Functional Networks');
        title(strjoin(["FNC of Component", num2str(ind.e(j))]), 'FontSize',16);
        set(ax(j,1), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;

        % Plot time series
        ax(j,2) = subplot(3, numel(ind.e), j+numel(ind.e)); hold on
        plot(1:N.TR*max(N.subjects{:,:}), ts); hold on
        title({"Time Course of ", strjoin(["Component", num2str(ind.e(j))], " ")}, 'FontSize',16);
        set(ax(j,2), {'XTick', 'XTickLabels', 'XLim', 'YLim'}, {0:N.TR:max(N.subjects{:,:})*N.TR, [], [0 max(N.subjects{:,:})*N.TR], [-lim.y lim.y]});
        legend(labels.diagnosis);
        ylabel("Activity");
        xlabel("Subjects");

        % Plot entropies
        ax(j,3) = subplot(3, numel(ind.e), j+(2*numel(ind.e))); hold on
        boxplot(squeeze(E(:,:,ind.e(j))), labels.diagnosis, 'Notch','on');
        ylim([min(E(:,:,ind.e(j)),[],'all')-2, max(E(:,:,ind.e(j)),[],'all','omitnan')+2]);
        title({"Group Entropies of", strjoin(["Component", num2str(ind.e(j))], " ")}, 'FontSize',16);
        ylabel("Shannon Entropy");
        plot(1:2, [max(E(:,:,ind.e(j)),[],'all','omitnan')+0.5, max(E(:,:,ind.e(j)),[],'all','omitnan')+0.5], 'k-', 'LineWidth',1);
        plot(1.5, max(E(:,:,ind.e(j)),[],'all','omitnan')+1, 'k*', 'MarkerSize',10);
    end
end
clear k ts E


%% Correlate entropy score against clinical variables

% % Collate matrix (table) of joint entropy scores & component entropeis
% edata = nan(nnz(isfinite(e.joint)), size(e.comp,1)+1);      % set table
% for j = 1:N.conditions
%     edata() =  e.joint();
% edata(~I{:,'SZ'},1) = e.joint(1:nnz(~I{:,'SZ'}), 1);        % select healthy controls
% edata(I{:,'SZ'},1) = e.joint(1:nnz(I{:,'SZ'}), 2);          % select patients
% end
% for i = 2:size(e.comp,1)+1
%     edata(~I{:,'SZ'},i) = e.comp(i-1, 1:nnz(~I{:,'SZ'}), 1);
%     edata(I{:,'SZ'},i) = e.comp(i-1, 1:nnz(I{:,'SZ'}), 2);
% end
% edata = array2table(edata, 'VariableNames',["joint"; strcat(repmat("Component ",[N.IC 1]), string(1:size(e.comp,1))')], 'RowNames',analysis_data.Properties.RowNames);
% clear i
% 
% % numerize clinical information
% adata = analysis_data{:,:};
% adata(:,2) = single(strcmpi(analysis_data{:,2}, "F"));
% adata(:,3) = single(strcmpi(analysis_data{:,3}, "SZ"));
% adata = double(adata);
% % adata = array2table(adata, 'VariableNames',analysis_data.Properties.VariableNames, 'RowNames',analysis_data.Properties.RowNames);
% 
% % Compute correlation matrix with joint entropy
% [rho, p.corr] = corr(adata, edata{:,:}, 'Rows','pairwise');
% 
% % Correct p-values for multiple comparison (family-wise error)
% p_bin = zeros(numel(p.corr), 1);
% for s = 1:length(p_bin)
% 	p_bin(sort(FDR_benjHoch(reshape(p.corr, [numel(p.corr),1]), 0.05)), 1) = 1;
% end; clear s
% p_bin = logical(squeeze(p_bin));
% p_bin = reshape(p_bin, size(p.corr));
% [y,x] = find(p_bin);
% 

% %% Plot correlation values & significance
% 
% % Plot correlation values (matrix)
% F(N.fig) = figure; N.fig = N.fig+1;
% F(N.fig-1).Position = get(0,'screensize');
% l = numel(ind.e);
% ax(l+1,1) = subplot(2,1,1);
% imagesc(rho); colormap jet;
% clim([-max(abs(rho),[],'all'), max(abs(rho),[],'all')]); colorbar; hold on;
% xticks(1:3:size(edata,2));
% xticklabels(edata.Properties.VariableNames(1:3:size(edata,2)));
% yticks(1:size(adata,2));
% yticklabels(analysis_data.Properties.VariableNames);
% title("Correlation of Entropy to Clinical Variables", 'FontSize',18);
% 
% % Plot correlation values (bar chart)
% ax(l+1,2) = subplot(2,1,2);
% b = bar3(rho); hold on;
% xticks(1:3:size(edata,2));
% yticks(1:3:size(adata,2));
% xticklabels(num2str(0:size(edata,2)));
% yticklabels(num2str(0:size(edata,2)));
% for k = 1:length(b)
%     zdata = b(k).ZData;
%     b(k).CData = zdata;
%     b(k).FaceColor = 'interp';
% end
% title("Correlation Coefficients", 'FontSize',18);
% xlabel("Component Numbers");
% clear k
% 
% % Add significance markers to correlation bar chart
% plot3(ax(l+1,2), x,y, ones(numel(x),1)*1.5*max(abs(rho),[],'all'), '*k');
% 
% % Plot p-values
% F(N.fig) = figure; N.fig = N.fig+1;
% F(N.fig-1).Position = get(0,'screensize');
% ax(l+2,1) = subplot(2,1,1);
% imagesc(p.corr); colormap(ax(l+2,1), jet);
% clim([-max(abs(p.corr),[],'all'), max(abs(p.corr),[],'all')]); colorbar; hold on;
% xticks(1:3:size(edata,2));
% xticklabels(edata.Properties.VariableNames(1:3:size(edata,2)));
% yticks(1:size(adata,2));
% yticklabels(analysis_data.Properties.VariableNames);
% title("Correlation p-values of Entropy to Clinical Values", 'FontSize',18);
% xlabel("Component Numbers");
% 
% % Plot binarized p-values
% ax(l+2,2) = subplot(2,1,2);
% imagesc(p_bin); colormap(ax(l+2,2), gray); hold on
% xticks(1:3:size(edata,2));
% xticklabels(edata.Properties.VariableNames(1:3:size(edata,2)));
% yticks(1:size(adata,2));
% yticklabels(analysis_data.Properties.VariableNames);
% title("Significant Correlations of Entropy to Clinical Values", 'FontSize',18);
% xlabel("Component Numbers");
% % Ax = gca;
% % xt = Ax.XTick';
% % xtv = compose('%.0f',xt)';
% % xt = xt(ind.e);
% % xtv = xtv(ind.e);
% % text(xt,zeros(size(xt)), xtv, 'Color','r', 'Horiz','center', 'Vert','top')
% 
% % Get signficant correlation coefficients & variables
% corrvals = nan(numel(find(p_bin)),4);
% [corrvals(:,1), corrvals(:,2)] = find(p_bin);
% corrvals(:,3) = rho(p_bin);
% corrvals(:,2) = corrvals(:,2)-1;
% corrvals(:,4) = p.corr(p_bin);
% corrvals = num2cell(corrvals);
% corrvals(:,1) = analysis_data.Properties.VariableNames(cell2mat(corrvals(:,1)));
% corrvals = cell2table(corrvals, 'VariableNames',["Variable", "Component", "rho", "p"]);
clear n sm s c k j ts E


%% Save results & figure(s)

% Save figures
savefig(F, fullfile(path{4}, fileName), 'compact');
for c = 1:numel(F)
    saveas(F(c), fullfile(path{4}, "Figures", strjoin([fileName, num2str(c)], '-')), 'svg');
end
clear c F

% Save files
N.fig = N.fig - 1;
clear ts;
save(fullfile(path{4}, fileName));