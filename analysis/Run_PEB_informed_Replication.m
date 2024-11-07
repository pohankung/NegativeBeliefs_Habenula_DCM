%% Parametric Empirical Bayes (PEB), Replication model (informed)

% This script reproduces the habenula effective connectivity results
% reported in the Kung et al. manuscript (Table 4, Supplementary Table 6).

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

% Load GCM & design matrix
load('../data/GCM_Replication.mat');
load('../dm/M_Replication.mat');

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
save('./BMA_search_AB_Replication.mat', 'BMA'); % output used for visualisation


% Review BMA results
% -----------------------------------------------------------------------
% Second-level effect - Mean (Table 4a, Supplementary Table 6)
%   Threshold: Free energy, Strong evidence (Pp>.95)
%   Display as matrix: 
%     1) A-matrix (endogenous connectivity)
%     2) B-matrix (modulatory connectivity; Input Chal OR Rep)
% -----------------------------------------------------------------------
spm_dcm_peb_review(BMA, DCM)


% Produce comparison plots of modulatory connectivity parameters of the
% Discovery and Replication models (Figure 4c)
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

% Extract model parameters
reference = 'Discovery';
load(['./BMA_search_AB_',reference,'.mat']);
    rBMA = BMA;
    rE = rBMA.Ep(11:22); % only consider modulatory connectivity
    rBMA.Cp = diag(rBMA.Cp);
    rC = rBMA.Cp(11:22);
    rj = find(rBMA.Pp(11:22)); %pull index of connectivity with non-zero Pp
    x = '';
    s = '--r';
    ref_line_col = [101 99 96]/256;
    comp_bar_col = [255 190 106]/256;

comparison = 'Replication';
load(['./BMA_search_AB_',comparison,'.mat']);
    cBMA = BMA;
    cE = cBMA.Ep(9:11); % only consider modulatory connectivity
    cBMA.Cp = diag(cBMA.Cp);
    cC = cBMA.Cp(9:11);
    cj = 1:length(cC);
    ep_names = {'Hb-PCC';'Hb-pOFC';'Hb-PCC'};

all_Pp = [cBMA.Pp(9:11); rBMA.Pp(rj)];
all_Pp_flipped = flip(all_Pp,1);

% Invoke custom spm_plot_ci.m to plot reference model parameters
figure('Position', [300 300 1200 600]),
spm_plot_ci(rE,rC,x,rj,s);
hold on

% Overlay comparison model parameters
spm_plot_ci(cE,cC,x,cj,s);

    % Find and adjust covariance line
    covar_lines = findobj(gca, 'Type', 'line');
    
        % Reference model
        for i = length(rj)+1:length(rj)*2
            set(covar_lines(i), 'LineWidth', 2, 'Color', ref_line_col,...
                'Marker', '_');
        end
        % Comparison model (hidden)
        for i = 1:length(cj)
            set(covar_lines(i), 'Visible', 'off');
        end
    
    % Find and adjust posterior expectation bar
    ep_bars = findobj(gca, 'Type', 'bar');
    
        % Reference model
        set(ep_bars(2), 'BarWidth', 0.4, 'EdgeColor', 'none',...
            'FaceAlpha', 1)
        % Comparison model
        set(ep_bars(1), 'BarWidth', 0.8, 'FaceColor', comp_bar_col,...
            'EdgeColor', 'none', 'FaceAlpha', 0.6)

% Add legends and labels
legend([ep_bars(2), ep_bars(1)], {reference, comparison}, FontSize=26);

xlim([0 (length(rj)+1)]);
xlabel('Modulatory connectivity', FontSize=30, FontWeight='bold');
xlabelHandle = get(gca, 'XLabel');
xlabelPosition = get(xlabelHandle, 'Position');
xlabelPosition(2) = xlabelPosition(2) - 0.045;
set(xlabelHandle, 'Position', xlabelPosition);
set(gca,'XTickLabel',ep_names,'fontsize',24,'XTickLabelRotation', 0);

ylabel('Posterior expectation', FontSize=30, FontWeight='bold');
set(gca,'YTickLabel',get(gca, 'YTickLabel'),'fontsize',24);
set(findall(gca, '-property', 'FontName'), 'FontName', 'Arial');

% Add asterisk to mark connectivity showing consistent estimates across
% models (i.e., Replication model connectivity with Pp > .95 while also
% falling within the 95% CI of the Discovery model connectivity estimate)
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
        text(covar_lines(i).XData(1), a_loc, '*', 'FontSize', 32,...
            'HorizontalAlignment', 'center', 'VerticalAlignment',...
            a_align);
    end
end

% Additional visualisation features
line([2.5 2.5], [0 1.2], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'HandleVisibility', 'off');

hold off
saveas(figure(4), ['./',reference,'_',comparison,'_EpCp.png']);