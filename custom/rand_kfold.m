%% Ramdomised stratification for 5-fold validation analysis

% -----------------------------------------------------------------------
% This script calls the subject lists for the Discovery and Replication
% samples, which are stored in the respective .mat files, perform
% randomised stratification, and generate the design matrix for DCM
% analysis.
% -----------------------------------------------------------------------

clear

num_fold = 5; % specify number of folds

% Generate random list of participants for each stratified sample
D_list = load('list_Discovery.mat').list;
R_list = load('list_Replication.mat').list;
full_list = cell2mat([D_list; R_list]);

ind = crossvalind('Kfold',full_list,num_fold);

% Save subject ID in new .mat for DCM estimation & create design matrices
for G = 1:num_fold

    list = num2cell(full_list(ind==G));
    save(['./list_k',num2str(G),'.mat'], 'list');

    dm.X = ones(length(list),1);
    dm.labels = {'Mean'};
    save(['../dm/M_validation_k',num2str(G),'.mat'], 'dm');

end

