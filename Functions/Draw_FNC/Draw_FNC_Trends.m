%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Draw pair FNC
%%%% Daniel
%%%% Trends
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = Draw_FNC_Trends(FC_input, domain, domain_ICN, lw, Bar_range)
% FC_input: is a N x N matrix, where N is the number of ROI
% domain: is the label of ICN, if drawing ROI based matrix, just set it as [1:N];
% domain_ICN: 
% lw: line widths for interregional and inter-domain borders.
%   First element is interregional width
%   Second element is inter-domain width
% Bar_range: the max and min of the colorbar. For example, [-0.5 0.5].

if nargin < 4
    Bar_range = [-1 1];
end

ICN_label = [];
line_idx = zeros(length(domain)-1,1);
for i = 1:length(domain)
    ICN_label = [ICN_label;domain_ICN{i}];
    if i ~= length(domain)
        line_idx(i) = length(ICN_label);
    end
end
%% draw figure

%figure

% load FC matrix
FC_draw = FC_input;
[X_range, ~] = size(FC_draw);

% figure size
%set(gcf,'Position',[0 0 850 800]);
%set(gca,'Position',[.06 .06 .90 .90]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% previous codes 
% set XY axis range
% imagesc(FC_draw)
% colorbar;
% set(gca,'Clim', Bar_range);
% set(gca,'ytick',[],'yticklabel',[])
% set(gca,'xtick',[],'xticklabel',[])

% save and load axis
% SaveAxis = axes('Position',get(gca,'Position'),'XAxisLocation','top','YAxisLocation','right');
% set(SaveAxis,'XTick',[],'YTick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% draw figure
figure_handle = imagesc(FC_draw);
colorbar;
set(gca,'Fontsize', 14)
set(gca,'Clim', Bar_range);

%% intersect colmap
% if opt_inter == 1
%     inter_num = 4;
%     new_colmap = zeros(64+63*(inter_num-1),3);
%     for i = 1:(size(colmap,1)-1)
%         diff_map = colmap(i+1,:) - colmap(i,:);
%         for j = 1:4
%         new_colmap((i-1)*4+j,:) = colmap(i,:) + (j-1)*diff_map/4;
%         end
%     end
%     new_colmap(end,:) = colmap(end,:);
% end
% colmap = new_colmap;

%% set the zero value to gray
% if (0 >= Bar_range(1)) && (0 <= Bar_range(2))
% %     color_idx = round(((0-Bar_range(1))./(Bar_range(2) - Bar_range(1)))*size(colmap,1));
% %     colmap(color_idx,:) = [1 1 1];
% %     colmap(color_idx+1,:) = [1 1 1];
% %     colmap(color_idx-1,:) = [1 1 1];
%       inte_bar = (Bar_range(2)-Bar_range(1))./size(colmap,1);
%       color_idx = round((0-Bar_range(1))/inte_bar)+1;
%       colmap(color_idx,:) = [1 1 1];
%       colmap(color_idx+1,:) = [1 1 1];
%       colmap(color_idx-1,:) = [1 1 1];
% end

%% gradient map
ini_idx = 20; % small:edge close to black;
zero_modified_idx = 2;
ratio = (Bar_range(2) - 0)/(0 - Bar_range(1));
num_hot = 100;
num_col = round((num_hot-zero_modified_idx)/ratio);
% hot map
colormap(rand(num_hot,3))
hot_map = colormap(hot);
hot_map_use = hot_map(end-zero_modified_idx:-1:(1+ini_idx),:);
% col map
colormap(rand(num_col,3))
hot_map = colormap(hot);
col_map_use = zeros(size(hot_map,1)-round(ini_idx/ratio),size(hot_map,2));
col_map_use(:,1) = hot_map((1+round(ini_idx/ratio):end),3);
col_map_use(:,2) = hot_map((1+round(ini_idx/ratio):end),2);
col_map_use(:,3) = hot_map((1+round(ini_idx/ratio):end),1);

colmap = [col_map_use;hot_map_use];
%%
colormap(colmap);
set(gca,'YTick',[])
set(gca,'XTick',[])

% diagonal entries
axis_ratio = (1:1:X_range);
for i = 1:X_range
    tx = num2str(ICN_label(i));
    text(axis_ratio(i),axis_ratio(i),tx,'Fontsize',8,'HorizontalAlignment','center');
end

% draw lines
box off
axis off
draw_idx = [line_idx+1;1;X_range+1];
for i = 1:(X_range+1)
    hold on
    Linerange = ((0.5):1:(X_range+0.5));
    axisposi  = (i-0.5).*ones(length(Linerange),1);
    if (i ~= draw_idx) & lw(1) > 0
        plot(Linerange, axisposi, 'b', 'LineWidth', lw(1));
        plot(axisposi, Linerange, 'b', 'LineWidth', lw(1));
    elseif (i ~= draw_idx)
    else
        plot(Linerange, axisposi, 'k', 'LineWidth', lw(2));
        plot(axisposi, Linerange, 'k', 'LineWidth', lw(2));
    end
end

% draw domain
for i = 1:length(domain)
    tx = domain{i};
    if i == 1
        axis_domain = length(domain_ICN{i})/2 + 0.5;
    else
        axis_domain = line_idx(i-1) + length(domain_ICN{i})/2 + 0.5;
    end
    text(-0.5, axis_domain, tx, 'Fontsize',14,'HorizontalAlignment','center', 'Rotation', 90);
    text(axis_domain, length(ICN_label) + 1.8, tx, 'Fontsize',14,'HorizontalAlignment','center');
end

