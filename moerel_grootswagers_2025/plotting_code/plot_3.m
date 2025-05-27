%% FIGURE 3: We're using real RDMs here, compare reconstruction methods
% Load the reconstruction accuracy  
load('analysis_3.mat');

%% Plot
p_col = flipud(winter(3));
cm = parula(199);
cm = [1,1,1;cm];
font_size = 13;
x_idx = 1:80;

fh = figure(2);clf
fh.Position = [100,100,600,800];

%%% PLOT THE CORRELATION WITH THE ORIGINAL RDM
% Use the full RDMs or the reconstructed part only
x_pos = [0.08,0.58];
y_pos = fliplr(linspace(0.055,0.78,4));
x_pos_title = [0.008,0.508];
titles = [{'A) Missing values'},{'B) Full matrix'}];
dataset_title = [{'Bracci 2019 (27 × 27)'},{'Robinson 2025 (36 × 36)'},{'Mur 2013 (92 × 92)'},{'Grootswagers 2024 (256 × 256)'}];
for dataset_idx = 1:4
    for plot_type = 1:2 % Correlation with full RDM or reconstructed part only
    
        % Title
        ax_t = axes('Position',[x_pos_title(plot_type),0.94,0.39,0.1]);
        text(0,0.5,titles{plot_type},'FontWeight','bold','FontSize',font_size+2);
        ax_t.Visible = 'off';

        % Select the data to plot
        if plot_type==1
            toplot = squeeze(corr_orig_reconstructed(:,:,dataset_idx,:));
            chance = zeros(1,numel(x_idx));
        else
            toplot = squeeze(corr_orig_full(:,:,dataset_idx,:));
            chance = 1-x_idx/100;
        end
    
        % Main plot
        ax = axes('Position',[x_pos(plot_type),y_pos(dataset_idx),0.39,0.165],'LineWidth',1.5);
        hold on;
        ph = gobjects(1,4);
        % Plot chance: substitute with random values
        ph(4) = plot(x_idx,chance,'--k','LineWidth',2);
        for i = 1:3 % Loop over methods
            ci = prctile(toplot(:,:,i),[2.5 97.5]);
            h = fill([x_idx,fliplr(x_idx)],[ci(1,:),fliplr(ci(2,:))],p_col(i,:),'FaceAlpha',0.1,'EdgeColor','none');
            ph(i) = plot(x_idx,nanmean(toplot(:,:,i)),'Color',p_col(i,:),'LineWidth',2); % There are some nans for large % missing for the graph method
        end    
        ax.FontSize = font_size;
        xlim([1,80]);
        ax.XTick = 10:10:90;
        ylim([-0.05,1]);
        xlabel('Percentage missing');
        ylabel('Reconstruction acc. (\rho)');
        lh = legend(ph,'Geometric','Graph','MDS','Chance','Location','SW');
        lh.Position(2) = lh.Position(2)+0.01;
        lh.BackgroundAlpha = 0.7;
        lh.EdgeColor = 'none';
        title(dataset_title{dataset_idx});
    end
end

% Save
saveas(gcf,'Figure3.png');




