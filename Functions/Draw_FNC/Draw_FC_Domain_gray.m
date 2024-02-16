%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Draw FC figures
%%%% Daniel
%%%% MRN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = Draw_FC_Domain_gray(FC_input, domain, domain_ICN, Bar_range, opt_inter, colmap)
% FC_input: is a N x N matrix, where N is the number of ROI
% ICN_label: is the label of ICN, if drawing ROI based matrix, just set it as [1:N];
% Bar_range: the max and min of the colorbar. For example, [-0.5 0.5].
if nargin < 5
    Bar_range = [-1 1];
    colmap = colormap(jet);
    close all
elseif nargin < 6
    colmap = colormap(jet);
    close all
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
figure
% load FC matrix
FC_draw = zeros(size(FC_input,1)+1,size(FC_input,2)+1);
FC_draw(1:size(FC_input,1),1:size(FC_input,2)) = FC_input;
[X_range Y_range] = size(FC_draw(1:end-1,1:end-1));

% figure size
set(gcf,'Position',[50 50 950 900]);
set(gca,'Position',[.06 .06 .90 .90]);

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
figure_handle = pcolor(FC_draw);
set(figure_handle, 'linestyle','none')
colorbar;
set(gca,'Clim', Bar_range);
%% intersect colmap
if opt_inter == 1
    inter_num = 4;
    new_colmap = zeros(64+63*(inter_num-1),3);
    for i = 1:(size(colmap,1)-1)
        diff_map = colmap(i+1,:) - colmap(i,:);
        for j = 1:4
        new_colmap((i-1)*4+j,:) = colmap(i,:) + (j-1)*diff_map/4;
        end
    end
    new_colmap(end,:) = colmap(end,:);
end
colmap = new_colmap;

%% set the zero value to gray
if (0 >= Bar_range(1)) && (0 <= Bar_range(2))
%     color_idx = round(((0-Bar_range(1))./(Bar_range(2) - Bar_range(1)))*size(colmap,1));
%     colmap(color_idx,:) = [1 1 1];
%     colmap(color_idx+1,:) = [1 1 1];
%     colmap(color_idx-1,:) = [1 1 1];
      inte_bar = (Bar_range(2)-Bar_range(1))./size(colmap,1);
      color_idx = round((0-Bar_range(1))/inte_bar)+1;
      colmap(color_idx,:) = [1 1 1];
      colmap(color_idx+1,:) = [1 1 1];
      colmap(color_idx-1,:) = [1 1 1];
end

%%
colormap(colmap);
set(gca,'YTick',[])
set(gca,'XTick',[])

% diagnol entries
axis_ratio = (1:1:X_range);
for i = 1:X_range
    tx = num2str(ICN_label(i));
    text(axis_ratio(i)+0.5,axis_ratio(i)+0.5,tx,'Fontsize',8,'HorizontalAlignment','center');
end

% draw lines
box off
axis off
draw_idx = [line_idx+1;1;X_range+1];
for i = 1:(X_range+1)
    hold on
    Linerange = ((1):1:(X_range+1));
    axisposi  = (i).*ones(length(Linerange),1);
    if (i ~= draw_idx)
        plot(Linerange, axisposi, 'b', 'LineWidth', 0.5);
        plot(axisposi, Linerange, 'b', 'LineWidth', 0.5);
    else
        plot(Linerange, axisposi, 'k', 'LineWidth', 1.5);
        plot(axisposi, Linerange, 'k', 'LineWidth', 1.5);
    end
end

% draw domain
for i = 1:length(domain)
    tx = domain{i};
    if i == 1
        axis_domain = length(domain_ICN{i})/2 + 1;
    else
        axis_domain = line_idx(i-1) + length(domain_ICN{i})/2 + 1;
    end
    text(-0.5, axis_domain, tx, 'Fontsize',14,'HorizontalAlignment','center', 'Rotation', 90);
    text(axis_domain, length(ICN_label) + 2, tx, 'Fontsize',14,'HorizontalAlignment','center');
end

