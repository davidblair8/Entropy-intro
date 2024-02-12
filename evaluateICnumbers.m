function [ev, F] = evaluateICnumbers(N, ts)
% UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Run PCA decomposition
[~,~,pc.ltnt, pc.tsqrd, pc.expl, pc.mu] = pca(ts');

% Search around N.IC
clust = nan(N.TR*sum(N.subjects), numel(N.IC));
for k = 1:numel(N.IC)
    [whitesig, dewhiteM] = icatb_calculate_pca(ts, N.IC(k));
    [~, ~, A, ~] = icatb_icaAlgorithm(1, whitesig');
    A = dewhiteM*A;
    [~, clust(:,k)] = max(A,[],2);
%     activations = fastica(ts, 'numOfIC',N.IC(k), 'verbose','off');
%     [~, clust(k,:)] = max(activations);
end
clust = clust';
clear k

% Cluster evaluations
ev.sil = evalclusters(ts', clust', 'silhouette');
ev.var = evalclusters(ts', clust', 'CalinskiHarabasz');
ev.db = evalclusters(ts', clust', 'DaviesBouldin');

% Plot variance explained per component
F = figure;
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

end