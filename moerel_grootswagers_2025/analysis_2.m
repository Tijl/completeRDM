%% We're using simulated RDMs here
clearvars; close all; clc

% Set the random seed (so we get the same random RDMs each time)
rng(1);
num_iter = 1000;            % Number of iterations (i.e. random cells masked)
perc_missing = 1:80;        % percentage of missing values we want to test

% We are using RDMs of different sizes
all_rdm_sizes = [32,64,128,256];    % The different sizes of the RDMs we want to test
num_rdm_sizes = numel(all_rdm_sizes);

% Pre-allocate the output
corr_orig_full = zeros(num_iter,numel(perc_missing),num_rdm_sizes);
corr_orig_reconstructed = zeros(num_iter,numel(perc_missing),num_rdm_sizes);

% Loop over the different RDMs sizes
for rdm_size_idx = 1:num_rdm_sizes

    % Set the rdm size
    fprintf('rdm size %.0f of %.0f\n',rdm_size_idx,num_rdm_sizes);

    % Simulate the RDM
    Xt = pdist(randn(all_rdm_sizes(rdm_size_idx),2));

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
            corr_orig_full(iter,perc_missing_idx,rdm_size_idx) = corr(squareform(Y)',Xt',"type","Pearson");
    
            % Check the correlation between reconstructed values only for the original and completed RDMs
            X_2d = Xt'; % Full RDM
            Y_2d = squareform(Y)'; % Reconstructed RDM
            nan_idx = isnan(squareform(X)'); % Get the missing/reconstructed values
            corr_orig_reconstructed(iter,perc_missing_idx,rdm_size_idx) = corr(X_2d(nan_idx),Y_2d(nan_idx),"type","Pearson");
    
        end
  
    end
end

% Save 
save('analysis_2.mat','corr_orig_full','corr_orig_reconstructed');
