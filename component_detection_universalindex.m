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


%% Load and sort data

% Load formatted dFNC data
load(fullfile(path{3}, 'FBIRN_DFNC_table.mat'));

% Confirm that IDs, data are properly indexed
assert(all(str2double(string(cell2mat(analysis_ID))) == str2double(analysis_data.Properties.RowNames)), "Data labels are not properly ordered!");
clear analysis_ID analysis_SCORE
assert(all(strcmpi(string(FILE_ID), string(analysis_data.Properties.VariableNames))), "Clinical variables are not properly ordered!");
clear FILE_ID

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
analysis_data = renamevars(analysis_data, ["diagnosis(1:sz; 2:hc)","gender(1:male; 2:female)"], ["Diagnosis","Gender"]);

% Replace numeric missing data code with NaN
analysis_data{:,5:end}(analysis_data{:,5:end} == -9999) = NaN;


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
N.IC = 8;


%% Define filename based on parameters

% Get number of ICs
fileName = fullfile("Domains", strcat(num2str(N.IC), "ICs"));

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


%% Concatenate and index time series

% rename FNC data
FNC.subj = DFNC_FBIRN;
FNC.full = cell2mat(DFNC_FBIRN);

% convert variables to row form
I.subject = str2double(string(analysis_data.Properties.RowNames)');
I.diagnosis = analysis_data{:,"Diagnosis"}';
I.gender = analysis_data{:,"Gender"}';
I.site = analysis_data{:,"Site"}';
I.age = analysis_data{:,"age"}';

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


%% Remove site mean from time series

% Preallocate arrays
A = cell(numel(ind.site),1);
lim.c = max(abs(FNC.full), [], 'all');
ax = gobjects(numel(ind.site), 1+2*N.conditions);

% plot mean group FNC matrices for each site
F(N.fig) = figure; F(N.fig).Position = get(0,'screensize');
% T(N.fig) = tiledlayout(F(N.fig), N.conditions, numel(ind.site));
N.fig = N.fig + 1;
for k = 1:numel(ind.site)
    for j = 1:N.conditions
        ax(k,j) = subplot(N.conditions, numel(ind.site), k+numel(ind.site)*(j-1));
        imagesc(ax(k,j), icatb_vec2mat(mean(FNC.full(I{:,'site'}==ind.site(k) & strcmpi(labels.diagnosis(j),I{:,'diagnosis'}),:), 1, "omitmissing")));
        colormap jet; clim([-lim.c lim.c]); pbaspect([1 1 1]); colorbar;
        ylabel('Functional Networks'); xlabel('Functional Networks');
        title([strcat(["Mean FNC of ", labels.diagnosis(j), ", Site ", num2str(k)])]);
    end
end

% remove sitewise means from FNC
F(N.fig) = figure; F(N.fig).Position = get(0,'screensize');
% T(N.fig) = tiledlayout(F(N.fig), 1+2*N.conditions, numel(ind.site));
N.fig = N.fig + 1;
for k = 1:numel(ind.site)
    % get sitewise dFNC values
    A{k} = FNC.full(I{:,'site'}==ind.site(k), :);
    A{k} = mean(A{k}, 1, "omitmissing");

    % get number of subjects per site
    n = nnz(analysis_data{:,"Site"}==ind.site(k));

    % plot mean site-wise FNC matrix for each site
    ax(k,N.conditions+1) = subplot(N.conditions+1, numel(ind.site), k);
    imagesc(icatb_vec2mat(A{k})); colormap jet;
    clim([-lim.c lim.c]); pbaspect([1 1 1]); colorbar;
    ylabel('Functional Networks'); xlabel('Functional Networks');
    title(strjoin(["Mean FNC of Site", num2str(k)]));

    % subtract mean site FNC matrix from time-resolved dFNC values
    FNC.full(I{:,'site'}==ind.site(k),:) = FNC.full(I{:,'site'}==ind.site(k),:) - A{k};
    FNC.subj(analysis_data{:,"Site"}==ind.site(k)) = cellfun(@minus, FNC.subj(analysis_data{:,"Site"}==ind.site(k)), repmat({repmat(A{k}, [N.TR 1])}, [n 1]), 'UniformOutput', false);
    
    % plot mean site-wise group FNC after removing site means
    for j = 1:N.conditions
        ax(k,j+N.conditions+1) = subplot(N.conditions+1, numel(ind.site), k+numel(ind.site)*j);
        imagesc(icatb_vec2mat(mean(FNC.full(I{:,'site'}==ind.site(k) & strcmpi(labels.diagnosis(j),I{:,'diagnosis'}),:), 1, "omitmissing")));
        colormap jet; clim([-lim.c lim.c]); pbaspect([1 1 1]); colorbar;
        ylabel('Functional Networks'); xlabel('Functional Networks');
        title([strcat(["Mean FNC of ", labels.diagnosis(j), ", Site ", num2str(k)]), strjoin([" - Mean of Site", num2str(k)])]);
    end
end
clear k j A n


%% Visualise static FNC matrix with indices

% subject-level sFNC
sFNC{1} = cell2mat(cellfun(@mean, DFNC_FBIRN, 'UniformOutput', false));
sFNC{1} = icatb_vec2mat(sFNC{1});       % Convert to square matrices

% group-level sFNC
sFNC{2} = nan(N.conditions, N.ROI, N.ROI);
F(N.fig) = figure; F(N.fig).Position = get(0,'screensize');
N.fig = N.fig + 1;
for c = 1:N.conditions
    % compute mean group-level FNC for condition
    sFNC{2}(c,:,:) = mean(sFNC{1}(strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c)),:,:),1,"omitmissing");

    % Visualize group sFNC
    a(1,c) = subplot(1, N.conditions, c);
    imagesc(squeeze(sFNC{2}(c,:,:))); colormap jet;
    clim([-lim.c lim.c]); pbaspect([1 1 1]); colorbar;
    ylabel('Functional Networks'); xlabel('Functional Networks');
    title(strjoin(["Mean FNC of", labels.diagnosis(c)]));
end
clear c

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


%% Remove site means from ICA components
% 
% F(N.fig) = figure; F(N.fig).Position = get(0,'screensize');
% N.fig = N.fig + 1;
% for k = 1:N.samplot
%     % Plot ICA source matrices
%     ax(k,1) = subplot(3, N.samplot, k);
%     imagesc(icatb_vec2mat(W(C(k),:))); colormap jet; clim([-limits.color limits.color]);
%     pbaspect([1 1 1]); colorbar; hold on
%     ylabel('Functional Networks'); xlabel('Functional Networks');
%     title(strjoin(["FNC of Component", num2str(C(k))]), 'FontSize',16);
% 
%     % Plot sample ICA time series pre-correction
%     ax(k,2) = subplot(3, N.samplot, N.samplot+k);
%     plot(1:N.TR*max(N.subjects{:,:}), ts.sources.ic{C(k)}); hold on
%     title({"Uncorrected Time Course of", strjoin(["Component", num2str(C(k))], " ")}, 'FontSize',16);
%     set(ax(k,2), {'XTick', 'XTickLabels', 'XLim', 'YLim'}, {0:N.TR:max(N.subjects{:,:})*N.TR, [], [0 max(N.subjects{:,:})*N.TR], [-limits.y limits.y]});
%     legend(labels.diagnosis);
%     ylabel("Activity");
%     xlabel("Subjects");
% end
% 
% % Remove site means from ICA components
% sitemean = nan(numel(ind.site), N.IC);
% sm = sources;
% for k = 1:numel(ind.site)
%     sitemean(k,:) = mean(sources(I.FNC.concat{:,"Site"}==ind.site(k),:), 1,'omitnan');
%     sm(I.FNC.concat{:,"Site"}==ind.site(k),:) = sm(I.FNC.concat{:,"Site"}==ind.site(k),:) - repmat(sitemean(k,:), [nnz(I.FNC.concat{:,"Site"}==ind.site(k)) 1]);
% end
% 
% % Recollate time series by IC
% for k = 1:N.IC
%     for c = 1:N.conditions
%         ts.sources.ic{k}(c,1:N.TR*N.subjects{:,labels.diagnosis(c)}) = sm(strcmpi(I.FNC.concat{:,"Diagnosis"},labels.diagnosis(c)),k)';
%     end
% end
% 
% % Recollate time series by subject
% for c = 1:N.conditions
%     for s = 1:N.subjects{:,c}
%         ts.sources.subj{s,c} = sm(strcmpi(I.FNC.concat{:,"Diagnosis"},labels.diagnosis(c)) & I.FNC.concat{:,"Subject"}==ind.subject(s,c), :);
%     end
% end
% 
% % Plot ICA time series post-correction
% for k = 1:N.samplot
%     ax(k,3) = subplot(3, N.samplot, (2*N.samplot)+k);
%     plot(1:N.TR*max(N.subjects{:,:}), ts.sources.ic{C(k)}); hold on
%     title({"Corrected Time Course of", strjoin(["Component", num2str(C(k))], " ")}, 'FontSize',16);
%     set(ax(k,3), {'XTick', 'XTickLabels', 'XLim', 'YLim'}, {0:N.TR:max(N.subjects{:,:})*N.TR, [], [0 max(N.subjects{:,:})*N.TR], [-limits.y limits.y]});
%     legend(labels.diagnosis);
%     ylabel("Activity");
%     xlabel("Subjects");
% end


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
p = table('Size',[N.IC+1,3], 'VariableTypes',{'double','double','double'}, 'VariableNames',["t" "ks" "perm"], 'RowNames',entropy.Properties.VariableNames);
h = table('Size',[N.IC+1,3], 'VariableTypes',{'logical','logical','logical'}, 'VariableNames',["t" "ks" "perm"], 'RowNames',entropy.Properties.VariableNames);
for c = 1:N.IC
    [h(c,:), p(c,:)] = sigtest(squeeze(E(:,:,c)));
end

% multiple-comparison correction
FDR = nan(N.IC,2);
Bonferroni = FDR;
Sidak = FDR;
[FDR(:,1), Bonferroni(:,1), Sidak(:,1)] = mCompCorr(N.IC, p{1:N.IC,"t"}, 0.05);
[FDR(:,2), Bonferroni(:,2), Sidak(:,2)] = mCompCorr(N.IC, p{1:N.IC,"ks"}, 0.05);
[FDR(:,3), Bonferroni(:,3), Sidak(:,3)] = mCompCorr(N.IC, p{1:N.IC,"perm"}, 0.05);

% convert NaNs to false
FDR(isnan(FDR)) = false;
Bonferroni(isnan(Bonferroni)) = false;
Sidak(isnan(Sidak)) = false;

% convert to logical arrays
Bonferroni = array2table(logical(Bonferroni), 'VariableNames',["Student's t", "Kolmogorov-Smirnov", "Permutation"], 'RowNames',entropy.Properties.VariableNames(1:N.IC));
Sidak = array2table(logical(Sidak), 'VariableNames',["Student's t", "Kolmogorov-Smirnov", "Permutation"], 'RowNames',entropy.Properties.VariableNames(1:N.IC));
FDR = array2table(logical(FDR), 'VariableNames',["Student's t", "Kolmogorov-Smirnov", "Permutation"], 'RowNames',entropy.Properties.VariableNames(1:N.IC));

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

% joint entropy comparisons
[h("Joint",:), p("Joint",:)] = sigtest(squeeze(E(:,:,N.IC+1)));
if h{"Joint",'t'}
    disp(strjoin(["Student's two-sample t-test shows significant difference between patient and control joint entropy (p = .", num2str(p{"Joint",'t'}), ")."], " "));
end
if h{"Joint",'ks'}
    disp(strjoin(["Kolmogorov-Smirnov two-tailed test shows significant difference between patient and control joint entropy (p = .", num2str(p{"Joint",'ks'}), ")."], " "));
end
if h{"Joint",'perm'}
    disp(strjoin(["Difference-of-means permutation test shows significant difference between patient and control joint entropy (p = .", num2str(p{"Joint",'perm'}), ")."], " "));
end
clear c k


%% Compile tables: entropy means, standard deviations

% estats.joint.mean = mean(entropy{:,"Joint"},'omitnan');
% estats.joint.mean = array2table(estats.joint.mean, 'VariableNames',labels.diagnosis);
% estats.comp.mean = squeeze(mean(e.comp,2,'omitnan'));
% estats.comp.mean = array2table(estats.comp.mean, 'VariableNames',labels.diagnosis);
% estats.joint.std = std(e.joint,0,'omitnan');
% estats.joint.std = array2table(estats.joint.std, 'VariableNames',labels.diagnosis);
% estats.comp.std = squeeze(std(e.comp,0,2,'omitnan'));
% estats.comp.std = array2table(estats.comp.std, 'VariableNames',labels.diagnosis);


%% Visualise components with group-level changes

% Visualise joint entropy of both conditions
F(N.fig) = figure; N.fig = N.fig + 1;
F(N.fig-1).Position = get(0,'screensize'); hold on;
boxplot(squeeze(E(:,:,N.IC+1)), labels.diagnosis, 'Notch','on');
ylim([min(squeeze(E(:,:,N.IC+1)),[],'all',"omitmissing")-10, max(squeeze(E(:,:,N.IC+1)),[],'all',"omitmissing")+10]);
title("Joint Entropy", 'FontSize',16); ylabel("Joint Entropy");
if h{"Joint",'t'}
    axes = gca;
    plot(1:2, [max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+3, max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+3], 'k-', 'LineWidth',1);
    plot(1.5, max(squeeze(E(:,:,N.IC+1)),[],'all','omitnan')+5, 'k*', 'MarkerSize',10);
elseif h{"Joint",'ks'}
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
        imagesc(icatb_vec2mat(W(ind.e(j),:))); colormap jet; clim([-lim.color lim.color]);
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


%% Regress clinical variables against entropy

% Format clinical variables

% Set up multiple linear regression
e.joint;
e.comp;
I.sources.subj;
[b, stats] = robustfit();
[b,bint, r,rint, stats] = regress();


% %% Correlate entropy score against clinical variables
% 
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

% % Save figures
% savefig(F, fullfile(path{4}, fileName), 'compact');
% for c = 1:numel(F)
%     saveas(F(c), fullfile(path{4}, "Figures", strjoin([fileName, num2str(c)], '_')), 'svg');
% end
% clear c F
% 
% % Save files
% N.fig = N.fig - 1;
% clear ts;
% save(fullfile(path{4}, fileName));