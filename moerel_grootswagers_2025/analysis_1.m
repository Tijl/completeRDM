%% We're using real RDMs here
clearvars; close all; clc

% Set the random seed (so we get the same random RDMs each time)
rng(1);
num_iter = 1000;            % Number of iterations (i.e. random cells masked)
perc_missing = 1:80;        % percentage of missing values we want to test

% We are using RDMs from 4 different datasets
% Get the data in the same format
num_datasets = 4;
behav_rdms = cell(1,num_datasets);
% Bracci, S., Ritchie, J. B., Kalfas, I. & Op de Beeck, H. The ventral visual pathway represents animal appearance over animacy, unlike human behavior and deep neural networks. J. Neurosci. 1714–18 (2019) doi:10.1523/JNEUROSCI.1714-18.2019.
% Data can be found here: https://osf.io/k6qne
x = load('behav_rdms/bracci2019.mat');
behav_rdms(1) = {x.OSF.behaviData};
% Robinson, A. K., Grootswagers, T., Shatek, S. M., Behrmann, M., & Carlson, T. A. Dynamics of visual object coding within and across the hemispheres: Objects in the periphery. Science Advances, 11(1), eadq0889. (2025). 
% Data can be found here: https://osf.io/9us3y/files/4ebef5bf-1963-4768-b510-111afe9c0c1c
% We're using the 'concept task' RDM (stats_concept_similarity.mat)
x = load('behav_rdms/robinson2025concept.mat');
behav_rdms(2) = {x.RDMmean};
behav_rdms{2}(eye(size(behav_rdms{2}))==1) = 0; % Set the diagonal to 0
% Mur, M. et al. Human Object-Similarity Judgments Reflect and Transcend the Primate-IT Object Representation. Front Psychol 4, (2013).
% Data can be found here: https://github.com/rsagroup/rsatoolbox_matlab/tree/develop/Demos/92imageData
x = load('behav_rdms/mur2013.mat');
behav_rdms(3) = {mean(cat(3,x.rdms_behav_92(:).RDM),3)};
% Grootswagers, T., Robinson, A. K., Shatek, S. M. & Carlson, T. A. Mapping the dynamics of visual feature coding: Insights into perception and integration. PLOS Computational Biology 20, e1011760 (2024).
% Data can be found here: https://github.com/Tijl/features-eeg
x = load('behav_rdms/grootswagers2024.mat');
behav_rdms(4) = {x.RDMmean};
behav_rdms{4}(eye(size(behav_rdms{4}))==1) = 0; % Set the diagonal to 0
clear x;

% Pre-allocate the output
corr_orig_full = zeros(num_iter,numel(perc_missing),num_datasets);
corr_orig_reconstructed = zeros(num_iter,numel(perc_missing),num_datasets);

% Loop over datasets
for dataset_idx = 1:num_datasets

    % Show an update
    fprintf('rdm dataset %.0f of %.0f\n',dataset_idx,num_datasets);

    % Select the behavioural RDM
    Xt = squareform(behav_rdms{dataset_idx});
   
    % Loop over iterations
    for iter = 1:num_iter
        fprintf('   iteration %.0f of %.0f\n',iter,num_iter);
    
        % Loop over different percentages of missing values
        parfor perc_missing_idx = 1:numel(perc_missing) % percentage of missing values

            % Make the RDM with missing values
            X = nan;
            while max(conncomp(graph(max(squareform(X),0),'upper')))>1
                
                % Preallocate the RDM
                X = Xt;
    
                % Delete a precentage of the values
                X(randsample(1:numel(Xt),ceil(perc_missing(perc_missing_idx)*numel(Xt)/100))) = NaN;
            end

            %%% RECONSTRUCTION %%%
            % Make it 2D
            X = squareform(X);
            % Complete the missing values
            Y = complete_rdm(X,false);
    
            % Check the correlation between the full original and reconstructed RDMs
            corr_orig_full(iter,perc_missing_idx,dataset_idx) = corr(squareform(Y)',Xt',"type","Pearson");
    
            % Check the correlation between reconstructed values only for the original and completed RDMs
            X_2d = Xt'; % Full RDM
            Y_2d = squareform(Y)'; % Reconstructed RDM
            nan_idx = isnan(squareform(X)'); % Get the missing/reconstructed values
            corr_orig_reconstructed(iter,perc_missing_idx,dataset_idx) = corr(X_2d(nan_idx),Y_2d(nan_idx),"type","Pearson");
    
        end
    
    end
end

% Save 
save('analysis_1.mat','corr_orig_full','corr_orig_reconstructed');
