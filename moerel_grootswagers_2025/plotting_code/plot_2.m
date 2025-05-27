%% FIGURE 2: We're using simulated RDMs here
% Load the reconstruction accuracy  
load('analysis_2.mat');

%% Plot
% Set plotting parameters
p_col = cool(4);
cm = parula(199);
cm = [1,1,1;cm];
percentage_missing = [25,50,75];
font_size = 16;
x_idx = 1:80;

fh = figure(1);clf
fh.Position = [100,100,800,800];

%%% PLOT THE CORRELATION WITH THE FULL RDM
% Use the full RDMs or the reconstructed part only
x_pos = [0.07,0.57];
x_pos_title = [0.008,0.508];
titles = [{'A) Missing values'},{'B) Full matrix'}];
for plot_type = 1:2 % Correlation with full RDM or reconstructed part only

    % Title
    ax_t = axes('Position',[x_pos_title(plot_type),0.9,0.39,0.1]);
    text(0,0.5,titles{plot_type},'FontWeight','bold','FontSize',font_size+2);
    ax_t.Visible = 'off';

    % Select the data to plot
    if plot_type==1
        toplot = corr_orig_reconstructed;
        chance = zeros(1,numel(x_idx));
    else
        toplot = corr_orig_full;
        chance = 1-x_idx/100;
    end

    % Main plot
    ax = axes('Position',[x_pos(plot_type),0.555,0.39,0.37],'LineWidth',1.5);
    hold on;
    ph = gobjects(1,5);
    for i = 1:4 % Loop over datasets
        ci = prctile(toplot(:,:,i),[2.5 97.5]);
        h = fill([x_idx,fliplr(x_idx)],[ci(1,:),fliplr(ci(2,:))],p_col(i,:),'FaceAlpha',0.2,'EdgeColor','none');
        ph(i) = plot(x_idx,mean(toplot(:,:,i)),'Color',p_col(i,:),'LineWidth',2);
    end

    % Plot chance: substitute with random values
    ph(5) = plot(x_idx,chance,'--k','LineWidth',2);

    ax.FontSize = font_size;
    xlim([1,80]);
    ax.XTick = 10:10:90;
    ylim([-0.05,1]);
    xlabel('Percentage missing');
    ylabel('Reconstruction accuracy (\rho)');
    lh = legend(ph,'32 × 32','64 × 64','128 × 128','256 × 256','Chance');
    lh.Position(1:2) = [x_pos(plot_type)+0.01,0.6];
    legend boxoff;
end

%%% VISUALISE THE RECONSTRUCTION
rng(111);
rdm_size = 32;
orig_rdm = squareform(pdist(randn(rdm_size,2)/4)); % Rescale for plotting
orig_rdm(eye(size(orig_rdm))==1) = 0; % Set the diagonal to 0

% Make the RDM with missing values
x_pos= [0.03,linspace(0.38,0.74,3)-0.002];
y_pos = [0.185,0.01]+0.085;

% Show the original full & zoomed in RDM
ax = axes('Position',[x_pos(1),y_pos(2),0.32,0.32]);
imagesc(orig_rdm,[-0.01,1]);
hold on
colormap(cm)
th = title('Original (32 × 32)');
th.Units = 'normalized';
th.Position(2) = th.Position(2)+0.01;
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';
ax.YDir = 'normal';

% Show the missing & reconstructed values
orig_rdm = squareform(orig_rdm);
for percentage_missing_idx = 1:3
    % Plot the background
    ax_b = axes('Position',[x_pos(percentage_missing_idx+1)-0.01,y_pos(plot_type)-0.01,0.145+0.02,0.355]);
    fill([0,1,1,0],[0,0,1,1],[0.9,0.9,0.9],'EdgeColor','none');
    ax_b.Visible = 'off';

    % Preallocate the RDM
    X = orig_rdm;
    % Delete a precentage of the values
    X(randsample(1:numel(orig_rdm),ceil(percentage_missing(percentage_missing_idx)*numel(orig_rdm)/100))) = NaN;
    % Make it 2D
    X = squareform(X);
    % Complete the missing values
    Y = complete_rdm(X,false);
    % Combine for plotting
    rdm_example = cat(3,X,Y);

    for plot_type = 1:2
        % Main plot
        ax = axes('Position',[x_pos(percentage_missing_idx+1),y_pos(plot_type),0.145,0.145]);
        imagesc(rdm_example(:,:,plot_type),[-0.01,1]);
        if plot_type==1
            hold on
            annotation('arrow',[x_pos(percentage_missing_idx+1)+0.004,x_pos(percentage_missing_idx+1)+0.004],[y_pos(plot_type)-0.002,y_pos(plot_type)-0.026],'Color','k','LineWidth',3,'HeadStyle','cback1');
            annotation('arrow',[x_pos(percentage_missing_idx+1)+0.145-0.004,x_pos(percentage_missing_idx+1)+0.145-0.004],[y_pos(plot_type)-0.002,y_pos(plot_type)-0.026],'Color','k','LineWidth',3,'HeadStyle','cback1');
            sel_title = sprintf('%02d%% missing',percentage_missing(percentage_missing_idx));
        else
            sel_title = sprintf('Accuracy: %.02f',corr(orig_rdm',squareform(Y)',"type","Pearson"));
        end
        colormap(cm)
        th = title(sel_title);
        th.Units = 'normalized';
        th.Position(2) = th.Position(2)+0.01;
        ax.YAxis.Visible = 'off';
        ax.XAxis.Visible = 'off';
        ax.YDir = 'normal';
    end

end

% Add colorbar
cb = colorbar;
cb.Position = [0.91,y_pos(2),0.02,0.32];
cb.Label.String = 'Dissimilarity';
cb.Label.FontSize = font_size-4;
cb.Limits = [0,1];

% Title
ax_t = axes('Position',[x_pos_title(1),0.42,0.39,0.1]);
text(0,0.5,'C) Example reconstructions (simulated 32 × 32)','FontWeight','bold','FontSize',font_size+2);
ax_t.Visible = 'off';

% Save
saveas(gcf,'Figure2.png');



