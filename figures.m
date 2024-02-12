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


%% ICs

% reshape activities to allow easy plotting
ts_ic = nan(N.IC, N.TR, max(N.subjects), 2);
for c = 1:2
	ts_ic(:,:, 1:nnz(I{:,'HC'} == c-1), c) = activity(:,:,I{:,'HC'} == c-1);
end
ts_ic = reshape(ts_ic, [N.IC, N.TR*max(N.subjects), 2]);

% Plot connectivity maps
F(N.fig) = figure; N.fig = N.fig + 1;
fncmat = nan(N.ROI, N.ROI, N.IC);
for j = 1:N.IC
	ax(j,1) = subplot(2, 5, j); hold on
    fncmat(:,:,j) = icatb_vec2mat(W(j,:)');
	imagesc(fncmat(:,:,j)); colormap jet; c = colorbar;
    c.Limits = [-max(abs(c.Limits)), max(abs(c.Limits))];
    set(ax(j,1), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;
    pbaspect([1 1 1]);
	title(strjoin(["Connectivity Matrix of Component", num2str(j)], " "), 'FontSize',16);
	ylabel("Spatial Networks");
    xlabel("Spatial Networks");
end
clear j c

% Plot time courses
F(N.fig) = figure; N.fig = N.fig + 1;
fncmat = nan(N.ROI, N.ROI, N.IC);
for j = 1:N.IC
    ax(j,2) = subplot(2, 5, j); hold on
    plot(1:N.TR*max(N.subjects), squeeze(ts_ic(j,:,:)));
    title(strjoin(["Time Course of Component", num2str(j)], " "), 'FontSize',16);
    set(ax(j,2), {'XTick', 'XTickLabels', 'XLim'}, {0:N.TR:max(N.subjects)*N.TR, [], [0 max(N.subjects)*N.TR]}); hold on;
    legend(I.Properties.VariableNames);
    ylabel("Activity");
    xlabel("Subjects");
end
clear j c


%% Entropy

% Plot entropies of ICs
ind.e = find(FDR{:,"Kolmogorov-Smirnov"});
F(N.fig) = figure; N.fig = N.fig + 1;
for j = 1:N.IC
	a(j) = subplot(2, N.IC/2, j); hold on
	boxplot(squeeze(e.comp(j,:,:)), I.Properties.VariableNames, 'Notch','on');
	ylim([min(e.comp(j,:,:),[],'all')-2, max(e.comp(j,:,:),[],'all','omitnan')+2]);
	title(strjoin(["Component", num2str(j)], " "), 'FontSize',16);
	ylabel("Shannon Entropy");
end

for j = 1:numel(ind.e)
	plot(a(ind.e(j)), 1:2, [max(e.comp(ind.e(j),:,:),[],'all','omitnan')+0.5, max(e.comp(ind.e(j),:,:),[],'all','omitnan')+0.5], 'k-', 'LineWidth',1);
	plot(a(ind.e(j)), 1.5, max(e.comp(ind.e(j),:,:),[],'all','omitnan')+1, 'k*', 'MarkerSize',10);
end


%% Spatial Correlations

% Load relevant mapping matrices
w = load('/Users/David/Library/CloudStorage/GoogleDrive-dblair@gsu.edu/My Drive/Calhoun/Results/5ICs_iteration1.mat', 'W'); W{5} = w.W;
w = load('/Users/David/Library/CloudStorage/GoogleDrive-dblair@gsu.edu/My Drive/Calhoun/Results/10ICs_iteration1.mat', 'W'); W{10} = w.W;
w = load('/Users/David/Library/CloudStorage/GoogleDrive-dblair@gsu.edu/My Drive/Calhoun/Results/22ICs_iteration1.mat', 'W'); W{22} = w.W;
w = load('/Users/David/Library/CloudStorage/GoogleDrive-dblair@gsu.edu/My Drive/Calhoun/Results/171ICs_iteration1.mat', 'W'); W{171} = w.W;
clear w

% fill correlation & significance matrices
ind.w = find(~cellfun(@isempty, W));
nk = nchoosek(ind.w, 2);
C = cell(numel(unique(nk(:,1))), numel(unique(nk(:,2)))); p = C;
% for j = 1:size(nk,1)
%     [C{nk(j,:)}, p{nk(j,:)}] = corr(W{nk(j,1)}', W{nk(j,2)}');
% end
% clear j
[C{1,1}, p{1,1}] = corr(W{nk(1,1)}', W{nk(1,2)}');
[C{1,2}, p{1,2}] = corr(W{nk(2,1)}', W{nk(2,2)}');
[C{1,3}, p{1,3}] = corr(W{nk(3,1)}', W{nk(3,2)}');
[C{2,2}, p{2,2}] = corr(W{nk(4,1)}', W{nk(4,2)}');
[C{2,3}, p{2,3}] = corr(W{nk(5,1)}', W{nk(5,2)}');
[C{3,3}, p{3,3}] = corr(W{nk(6,1)}', W{nk(6,2)}');

% Get limits for color maps
l.c = cellfun(@abs,C, 'UniformOutput',false);
l.c = cellfun(@max, l.c, 'UniformOutput',false);
l.c = cellfun(@max, l.c, 'UniformOutput',false);
l.c = max(cell2mat(l.c(~cellfun(@isempty, l.c))));
l.p = cellfun(@abs,p, 'UniformOutput',false);
l.p = cellfun(@max, l.p, 'UniformOutput',false);
l.p = cellfun(@max, l.p, 'UniformOutput',false);
l.p = max(cell2mat(l.p(~cellfun(@isempty, l.p))));

% Plot correlation matrices
figure;
ind.c = find(~cellfun(@isempty, C));
for j = 1:nnz(~cellfun(@isempty, C))
    subplot(3, nnz(~cellfun(@isempty, C)), j);
    imagesc(C{ind.c(j)});
    colormap jet; c(j,1) = colorbar;
    clim([-l.c, l.c]); hold on;
    title(strjoin(["Correlation:", num2str(nk(j,1)), "tICs vs.", num2str(nk(j,2)), "tICs"]));
end
clear j

% Plot correlation p-values
ind.p = find(~cellfun(@isempty, p));
for j = 1:nnz(~cellfun(@isempty, p))
    subplot(3, nnz(~cellfun(@isempty, p)), j+nnz(~cellfun(@isempty, p)));
    imagesc(p{ind.p(j)});
    colormap jet; c(j,2) = colorbar;
    clim([0, l.p]); hold on;
    title(strjoin(["p-values:", num2str(nk(j,1)), "tICs vs.", num2str(nk(j,2)), "tICs"]));
end
clear j

% Correct p-values for multiple comparison (family-wise error)
for j = 1:nnz(~cellfun(@isempty, p))
    p_bin{ind.p(j)} = zeros(numel(p{ind.p(j)}), 1);
    for s = 1:length(p_bin)
	    p_bin{ind.p(j)}(sort(FDR_benjHoch(reshape(p{ind.p(j)}, [numel(p{ind.p(j)}),1]), 0.05)), 1) = 1;
    end; clear s
    p_bin{ind.p(j)} = logical(squeeze(p_bin{ind.p(j)}));
    p_bin{ind.p(j)} = reshape(p_bin{ind.p(j)}, size(p{ind.p(j)}));
    [y{ind.p(j)},x{ind.p(j)}] = find(p_bin{ind.p(j)});
    
    a = subplot(3, nnz(~cellfun(@isempty, p)), j+2*nnz(~cellfun(@isempty, p)));
    imagesc(p_bin{ind.p(j)}); colormap(a,'gray'); hold on
    title(strjoin(["FDR-Corrected p-values:", num2str(nk(j,1)), "tICs vs.", num2str(nk(j,2)), "tICs"]));
end
clear j
