%% Parametric Empirical Bayes (PEB), Effects of PTQ on habenula effective connectivity

% This script reproduces the habenula effective connectivity results
% reported in the Kung et al. manuscript (Supplementary Table 4).

% -----------------------------------------------------------------------
% Please ensure that your SPM12 folder (r7771) is listed in your MATLAB set
% path. These results were obtained using Matlab R2023a. Values may
% slightly differ from the manuscript depending on OS and Matlab version.
% -----------------------------------------------------------------------

% This section runs a PEB model assessing the association between
% perseverative thinking and habenula connectivity. The design matrix
% included PTQ total score as a covariate regressor in addition to the
% intercept term (single column of ones), which represent the average
% connectivity estimate across the sample.

clear
close all

% Load GCM & design matrix
load('../data/GCM_Discovery.mat');
load('../dm/M_Discovery_PTQ.mat');

X = dm.X;
K = width(X);
X(:,2:K)=X(:,2:K)-mean(X(:,2:K)); % all covariates are mean-centered
X_labels = dm.labels;
    
M = struct();
M.Q = 'all';
M.X = X;
M.Xnames = X_labels;
    
% Hierarchical (PEB) inversion of DCMs using BMR and Variational Laplace
[PEB, RCM] = spm_dcm_peb(DCM, M, {'A','B'});
    
% Hierarchical (PEB) model comparison and averaging
BMA = spm_dcm_peb_bmc(PEB);


% Review BMA results
% -----------------------------------------------------------------------
% Second-level effect - PTQ (Supplementary Table 4)
%   Threshold: Free energy, Strong evidence (Pp>.95)
%   Display as matrix: 
%     1) A-matrix (endogenous connectivity)
%     2) B-matrix (modulatory connectivity; Input Chal OR Rep)
% -----------------------------------------------------------------------
spm_dcm_peb_review(BMA, DCM)