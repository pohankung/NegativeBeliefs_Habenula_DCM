%% Parametric Empirical Bayes (PEB), k-fold validation models (informed)

% This script reproduces the habenula effective connectivity results
% reported in the Kung et al. manuscript (Supplementary Figure 4).

% -----------------------------------------------------------------------
% 1) Please ensure that your SPM12 folder (r7771) is listed in your MATLAB set
% path. These results were obtained using Matlab R2023a. Values may
% slightly differ from the manuscript depending on OS and Matlab version.
%
% 2) Hierarchical PEB model inversion needs to be completed prior to this
% step (i.e., Run_PEB_Discovery.m).
% -----------------------------------------------------------------------

clear
close all

% Loop through the k-fold validation samples to perform PEB inversion
% -----------------------------------------------------------------------
num_fold = 5; % 5-fold validation

for k_id = 1:num_fold
    
    % Load respective GCM & design matrix
    load(['../data/GCM_validation_k',num2str(k_id),'.mat']);
    load(['../dm/M_validation_k',num2str(k_id),'.mat']);

    X = dm.X;
    K = width(X);
    X(:,2:K)=X(:,2:K)-mean(X(:,2:K)); % all covariates are mean-centered
    X_labels = dm.labels;
        
    M = struct();
    M.Q = 'all';
    M.X = X;
    M.Xnames = X_labels;

    % Update 3rd-level prior expectation (and covariance) with BMA posterior
    % from the Discovery model
    % -----------------------------------------------------------------------
    % [PEB, RCM] = spm_dcm_peb(DCM, M, {'A','B'});
    % Main Inputs:
    % - DCM:    Structure array of inverted DCMs that describe the same network
    %            model explored in the Discovery dataset. 
    % - M:      2nd-level design matrix containing group effect and other related parameters.
    %   1) M.bE:    3rd-level prior expectation. Here, we input the posterior
    %               expectations obtained via BMR (i.e., BMA.Ep)
    %   2) M.bC:    3rd-level prior covariance. Here, we input the posterior
    %               covariances obtained via BMR (i.e., BMA.Cp)
    % -----------------------------------------------------------------------
    
    % Load Discovery PEB model (Bayesian-averaged) to extract posteriors
    load('../analysis/BMA_search_AB_Discovery.mat');
        
        % Replace generic priors with Discovery model posteriors
        % This step effectively limits PEB model inversion to only consider
        % connections with posterior probability >.95 in the Discovery model
    
        Ep_base = zeros(size(spm_vec(DCM{1,1}.M.pE)));
        for i = 1:size(BMA.Ep)
            ind = BMA.Pind(i); % grab index of connections with Pp >.95
            Ep_base(ind) = BMA.Ep(i); % plug posterior expectation into base array
        end
        bBMA_Ep = Ep_base;
        M.bE = bBMA_Ep;
        
        Cp_base = zeros(size(DCM{1,1}.M.pC));
        for i = 1:size(BMA.Cp)
            ind = BMA.Pind(i); % grab index of connections with Pp >.95
            Cp_base(ind,ind) = BMA.Cp(i,i); % plug posterior covariance into base array
        end
        bBMA_Cp = Cp_base;
        M.bC = bBMA_Cp;


    % Invert and summarise informed PEB model
    % -----------------------------------------------------------------------
    % Hierarchical (PEB) inversion of DCMs using BMR and Variational Laplace
    [PEB, RCM] = spm_dcm_peb(DCM, M, {'A','B'});
    
    % Hierarchical (PEB) model comparison and averaging
    BMA = spm_dcm_peb_bmc(PEB);
    save(['./BMA_search_AB_validation_k',num2str(k_id),'.mat'], 'BMA'); % output used for visualisation

end

close all

% Loop through the k-fold validation samples to perform PEB inversion to
% produce comparison plots of connectivity parameters of the Discovery
% model and 5-fold validation samples
% -----------------------------------------------------------------------
% spm_plot_ci(E,C,x,j,s)
%
% Plot mean and conditional confidence intervals
% FORMAT spm_plot_ci(E,C,x,j,s)
% E - expectation (structure or array)
% C - variance or covariance (structure or array)
% x - domain
% j - rows of E to plot
% s - string to specify plot type:e.g. '--r' or 'exp'
%
% If E is a row vector with two elements, confidence regions will be
% plotted; otherwise, bar charts with confidence intervals are provided
%
% Copyright (C) 2008-2015 Wellcome Trust Centre for Neuroimaging, adapted
% to generate 95% CI
% -----------------------------------------------------------------------
    
addpath ../custom/
num_fold = 5; % 5-fold validation

for k_id = 1:num_fold
    
    % Extract model parameters
    reference = 'Discovery';
    load(['./BMA_search_AB_',reference,'.mat']);
        rBMA = BMA;    
        rE = rBMA.Ep;
        rC = rBMA.Cp;
        rj = find(rBMA.Pp); % pull index of connectivity with non-zero Pp
        x = '';
        s = '--r';
        ref_line_col = [101 99 96]/256;
        comp_bar_col = [064 176 166]/256;
        
    comparison = ['k',num2str(k_id)];
    load(['./BMA_search_AB_validation_',comparison,'.mat']);
        cBMA = BMA;    
        cE = cBMA.Ep;
        cC = cBMA.Cp;
        cj = 1:length(cC);
        ep_names = {'Hb-Hb';'Hb-PCC';'PCC-Hb';'PCC-PCC';'HC-Hb';'HC-HC';...
                    'pOFC-Hb';'pOFC-pOFC';'Hb-PCC';'Hb-pOFC';'Hb-PCC'};
        
        all_Pp = [cBMA.Pp; rBMA.Pp(rj)];
        all_Pp_flipped = flip(all_Pp,1);
    
    % Invoke custom spm_plot_ci.m to plot reference model parameters
    figure('Position', [300 300 1350 500]),
    spm_plot_ci(rE,rC,x,rj,s);
    hold on
    
    % Overlay comparison model parameters
    spm_plot_ci(cE,cC,x,cj,s);
      
        % Find and adjust covariance line
        covar_lines = findobj(gca, 'Type', 'line');

            % Reference model
            for i = length(rj)+1:length(rj)*2
                set(covar_lines(i), 'LineWidth', 2, 'Color', ref_line_col, 'Marker', '_');
            end
            % Comparison model (hidden)
            for i = 1:length(cj)
                set(covar_lines(i), 'Visible', 'off');
            end
        
        % Find and adjust posterior expectation bar
        ep_bars = findobj(gca, 'Type', 'bar');
        
            % Reference model
            set(ep_bars(2), 'BarWidth', 0.4, 'EdgeColor', 'none', 'FaceAlpha', 1)
            % Comparison model
            set(ep_bars(1), 'BarWidth', 0.8, 'FaceColor', comp_bar_col,...
                'EdgeColor', 'none', 'FaceAlpha', 0.6)
    
    % Add legends and labels
    legend([ep_bars(2), ep_bars(1)], {reference, comparison}, FontSize=16);
    
    xlim([0 (length(rj)+1)]);
    xlabel('Connections', FontSize=24, FontWeight='bold');
    xlabelHandle = get(gca, 'XLabel'); 
    xlabelPosition = get(xlabelHandle, 'Position'); 
    xlabelPosition(2) = xlabelPosition(2) - 0.09; 
    set(xlabelHandle, 'Position', xlabelPosition); 
    set(gca,'XTickLabel',ep_names,'fontsize',16,'XTickLabelRotation', 0);
    
    ylabel('Posterior expectation', FontSize=24, FontWeight='bold');
    set(gca,'YTickLabel',get(gca, 'YTickLabel'),'fontsize',16);
    set(findall(gca, '-property', 'FontName'), 'FontName', 'Arial');
    
    % Add asterisk to mark connectivity showing consistent estimates across
    % models (i.e., Validation model connectivity with Pp > .95 while also
    % falling within the 95% CI of the Discovery model connectivity
    % estimate)
    for i = length(rj)+1:length(rj)*2
        ci_l = covar_lines(i).YData(1);
        ci_u = covar_lines(i).YData(2);
        if all_Pp_flipped(i)>0.95
            if ci_u < 0
                a_loc = ci_l;
                a_align = 'cap';
            else
                a_loc = ci_u;
                a_align = 'baseline';
            end
            text(covar_lines(i).XData(1), a_loc, '*', 'FontSize', 24,...
                'HorizontalAlignment', 'center', 'VerticalAlignment',...
                a_align);
        end
    end
    
    % Secondary axes label
    text(9.3, -0.66, 'Modulatory Connectivity', 'Color', [0.5 0.5 0.5], 'FontSize', 16, 'VerticalAlignment', 'top');
    text(3.4, -0.66, 'Intrinsic Connectivity', 'Color', [0.5 0.5 0.5], 'FontSize', 16, 'VerticalAlignment', 'top');
    
    % Additional visualisation features
    line([8.5 8.5], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'HandleVisibility', 'off');
    legend_position = get(legend, 'Position');
    legend_position(1) = 0.15; 
    legend_position(2) = legend_position(2);
    set(legend, 'Position', legend_position);
    
    hold off
    
    saveas(figure(k_id), ['./',reference,'_',comparison,'_EpCp.png']);

end

close all