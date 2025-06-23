This file gives an overview of technical considerations and scripts used in generating the results for my master thesis.

In this thesis, R version 4.5.0 and Python version 3.13.3 was used.

All script used in this thesis contain a description at the start, describing the inteded purpose of the script.

**Code folder summary**
Data folder: Contains various datafiles used in thesis. Note that large datafiles, such as genomic data, is not included here due to file size constraints.
GEMMA: Contains files necessary for running GEMMA (GWAS) programme.
GetGemma: Contains GEMMA results.
Illustration_viz: Contains scripts used for visualizing all results, such as generating figures used in thesis.
Models folder: Contains the XGBoost models with tuned hyperparameters used for prediction and SHAP analysis.
Markov folder: Contains script used for NTNU markov-server. Scripts here are implemented to run parallel processes.
Results folder: Contains correlation results from genomic prediction for all models.
Various scripts in Code folder: The Code folder contains all scripts used for generating results. This includes generation of CV, implementations of prediction models, implementation of inference methods and more.
