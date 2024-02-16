%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Draw FC figures
%%%% Daniel
%%%% MRN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = Draw_FC(FC_input, ICN_label, Bar_range, colmap)
% FC_input: is a N x N matrix, where N is the number of ROI
% ICN_label: is the label of ICN, if drawing ROI based matrix, just set it as [1:N];
% Bar_range: the max and min of the colorbar. For example, [-0.5 0.5].
if nargin < 2
    ICN_label = [1:size(FC_input),1];
    Bar_range = [-1 1];
    colmap = colormap(jet);
    close all
elseif nargin < 3
    Bar_range = [-1 1];
    colmap = colormap(jet);
    close all
elseif nargin < 4
    colmap = colormap(jet);
    close all
end
%% draw figure
figure
% load FC matrix
FC_draw = FC_input;
[X_range Y_range] = size(FC_draw);

% figure size
set(gcf,'Position',[50 50 950 900]);
set(gca,'Position',[.06 .06 .90 .90]);

% set XY axis range
imagesc(FC_draw)
colorbar;
set(gca,'Clim', Bar_range);
set(gca, 'YDir','reverse');
set(gca,'XLim',[0.5 X_range+0.5]);%X????????
set(gca,'XTick',[1:1:X_range]);%?????????
set(gca,'YLim',[0.5 Y_range+0.5]);%X????????
set(gca,'YTick',[1:1:Y_range]);%?????????
% set(gca,'XTickLabel',{'a','b','c'})

% save and load axis
SaveAxis = axes('Position',get(gca,'Position'),'XAxisLocation','top','YAxisLocation','right');
set(SaveAxis,'XTick',[],'YTick',[]);

% draw figure
imagesc(FC_draw)
set(gca,'Clim', Bar_range);
colormap(colmap);
set(gca,'YTick',[])
set(gca,'XTick',[])

% diagnol entries
axis_ratio = (1:1:X_range);
for i = 1:X_range
    tx = num2str(ICN_label(i));
    text(axis_ratio(i),axis_ratio(i),tx,'Fontsize',8,'HorizontalAlignment','center');
end

% draw lines
for i = 1:X_range
    hold on
    Linerange = ((1-0.5):1:(X_range+0.5));
    axisposi  = (i+0.5).*ones(length(Linerange));
    plot(Linerange, axisposi);
    plot(axisposi, Linerange);
end

