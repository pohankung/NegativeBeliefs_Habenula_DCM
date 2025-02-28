# Habenula neural circuitry drives negative self-cognitions (in prep)
The repository contains [de-identified effective connectivity data](/data/) and the codes used to support the main findings in the [Kung et al. manuscript](DOI). 

The present study leverages ultra-high field (UFH) MRI and dynamic causal model (DCM) to characterise habenula effective connectivity during the processing of negative self-cognitions commonly reported in psychopathology. Out-of-sample reliability of the current findings is assessed in an independent replication sample and a randomised stratified 5-fold validation analysis.

## Discovery model
The [Run_PEB_Discovery.m](analysis/Run_PEB_Discovery.m) script reproduces the DCM results of habenula effective connectivity in the discovery sample. 

The [Run_PEB_Discovery_CNBTQ.m](analysis/Run_PEB_Discovery_CNBTQ.m) and [Run_PEB_Discovery_PTQ.m](analysis/Run_PEB_Discovery_PTQ.m) scripts reproduces the DCM exploring the association between habenula effective connectivity and participant's endorsement of negative self-cognitions and tendency of perseverative thinking, respectively.

The [Run_PEB_Discovery_LeftHemisphere.m](analysis/Run_PEB_Discovery_LeftHemisphere.m) script reproduces the exploratory DCM results of the left-lateralised habenula connectivity model. 

## Replication and 5-fold validation models
The **Run_PEB_Replication.m** scripts reproduce habenula connectivity results in the replication sample. The prior distribution of the replication model is furnished with the posterior distribution of the discovery model, taking into account information on habenula connectivity that was obtained through the discovery model when inverting the replication PEB model[<sup>1</sup>](analysis/Run_PEB_informed_Replication.m). Key modulatory connectivity estimates from this "informed" replication model are compared to that of the discovery model and visualised using the same script.

A "non-informed" replication model with the generic prior distribution is also inverted[<sup>2</sup>](analysis/Run_PEB_noninformed_Replication.m), and the results are summarised in the Supplementary Information.

The script [Run_PEB_informed_kfold.m](analysis/Run_PEB_informed_kfold.m) reproduces the effective connectivity results across the 5-fold validation models and the comparison plot included in the Supplementary Information.

## Decisional bias to restructure or repeat self-cognitions

## Notes
1. A custom version of `spm_plot_ci` is used to calculate the 95% confidence interval of the model parameter estimates. This is included [here](custom/spm_plot_ci.m).
2. Randomised stratified samples for the 5-fold validation analysis are generated with `crossvalind` in the custom code [rand_kfold.m](custom/rand_kfold.m).

