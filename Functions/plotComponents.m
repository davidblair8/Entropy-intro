function [F] = plotComponents(memberships, thresh, ROI, cortex, origin, cind, N, fDims)
%	PLOTCOMPONENTS plots the spatial structure of each component 
% 

% Storage arrays
m = nan(N.ROI, N.ROI, N.IC);


%% Plot combined figures

% Open combined figures
F(1) = figure('Position', fDims(1,:));
F(2) = figure('Position', fDims(1,:)); tiledlayout(3, N.IC, 'TileSpacing','none', 'Padding','tight');

% Plot index labels on bar chart
figure(F(1));
subplot(1, N.IC+1, 1);
ylim([1 N.ROI]); yticks([1 3:2:N.ROI/2 1+N.ROI/2:2:N.ROI]); yticklabels(ROI{[1 3:2:N.ROI/2 1+N.ROI/2:2:N.ROI],'Label'});

% Bar charts
for k = 1:N.IC
    % Enforce consensus: smaller community should be "positive"
    if sum(sign(memberships)) > 0
        memberships(:,k) = -memberships(:,k);
    end
    
    % Plot component bar charts
    figure(F(1));
    ax(k,1) = subplot(1, N.IC+1, k+1); hold on;
    ind(:,1) = memberships(:,k) < -thresh;	% select node weights which surpass threshold
    ind(:,2) = memberships(:,k) > thresh; 	% select node weights which surpass threshold
    ind(:,3) = memberships(:,k) < 0;
    ind(:,4) = memberships(:,k) > 0;
    a = memberships(:,k); a(~ind(:,1)) = 0; barh(1:N.ROI, a, 'b');
    a = memberships(:,k); a(~ind(:,2)) = 0; barh(1:N.ROI, a, 'r');
    a = memberships(:,k); a(~ind(:,3)) = 0; barh(1:N.ROI, a, 'b', 'FaceAlpha',0.3);
    a = memberships(:,k); a(~ind(:,4)) = 0; barh(1:N.ROI, a, 'r', 'FaceAlpha',0.3);
    set(ax(k,1), {'YTick'}, {[]});
    title(num2str(k));
end

% Brain space
for k = 1:N.IC
    % Plot component connectivity matrix
    ax(k,3) = subplot(3, N.IC, k+N.IC);
    m(:,:,k) = zscore(squeeze(memberships(:,k))*squeeze(memberships(:,k))', 0, 'all');
    imagesc(m(:,:,k)); colormap jet;
    set(ax(k,3), {'XLim','YLim','XTick','YTick'}, {[0.5 N.ROI+0.5], [0.6 N.ROI+0.5], [], []}); hold on;
    pbaspect([1 1 1]);
    ax(k,3).CLim = [min(m,[],'all') max(m,[],'all')];
    
    % Plot component in brain space
    figure(F(2));
    ax(k,2) = subplot(3, N.IC, k); hold on;
    plot_nodes_in_cortex(cortex, memberships(:,k), ROI{:,{'x','y','z'}}, origin, thresh, [], cind, [], []);
end

% Label component bar charts
sgtitle(F(1), "Component Spatial Maps", 'FontSize',20);
% xlabel(ax(floor(N.IC/2)+1, 1), 'z-score', 'FontSize',15, 'HorizontalAlignment','right');
% ylabel(ax(1,1), 'Brain Areas', 'FontSize',15);
% set(ax(1,1), {'YTick','YTickLabel','FontSize'}, {1:N.ROI, ROI{:,"Label"}, 6});

% Set component map color range
for k = 1:N.IC
    ax(k,3).CLim = [min(m,[],'all') max(m,[],'all')];
end

% Add connectivity matrix colorbar
figure(F(2));
cax = subplot(3, N.IC, [1+2*N.IC 3*N.IC]); axis off
cax.CLim = [min(m,[],'all') max(m,[],'all')];
c = colorbar(cax, 'northoutside');
c.Label.String = 'Outer Product';

% Add connectivity matrix labels
% ylabel(ax(1,3), 'Brain Areas', 'FontSize',12);
% xlabel(ax(floor(N.IC/2)+1, 3), 'Brain Areas', 'FontSize',12, 'HorizontalAlignment','right');


%% Plot individual components

% Plot brain maps and connectivity matrices
for k = 1:N.IC
    F(k+2) = figure('Position', fDims(2,:));
    
    % Brain maps
    ax(k,4) = subplot(2, 1, 1); hold on;
	plot_nodes_in_cortex(cortex, memberships(:,k), ROI{:,{'x','y','z'}}, origin, thresh, [], cind, [], []);
    title("Cortical Connectivity");
    
    % Connectivity matrices
    ax(k,5) = subplot(2, 1, 2);
    ax(k,5).CLim = [min(m,[],'all') max(m,[],'all')]; hold on;
    imagesc(m(:,:,k)); colormap jet; colorbar;
    set(ax(k,5), {'XLim','XTick','XTickLabel','FontSize','XTickLabelRotation'}, {[0.5 N.ROI+0.5], 5:5:N.ROI, ROI{5:5:N.ROI,"Label"}, 7, -45});
    set(ax(k,5), {'YLim','YTick','YTickLabel','FontSize'}, {[0.6 N.ROI+0.5], 5:5:N.ROI, ROI{5:5:N.ROI,"Label"}, 7});
    title("Outer Product");
    pbaspect([1 1 1]);
    text('Position',N.ROI*[1.11 1.03], 'String',"z-score", 'Rotation',0, 'FontSize',7);
    
    % sgtitle(strjoin("Component", num2str(k)));
end

end

