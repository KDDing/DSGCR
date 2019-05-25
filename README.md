# DSGCR

Data.xlsx is the benchmark dataset.
BMA.R is for calculating the drug biological similarity based on best match average method.
CrossValidation.R is for five-fold cross validation, involving in constructing training set and testing set, as well as calculating the FPR and TPR based on the synergy score.
FFCV.R is for implementing Five-Fold Cross Validation (FFCV) experiment.
Lapla.R is for normalizing the similarity matrix and calculating the graph Laplacian matrix.
LOOCV.R is for implementing Leave-One-Out Cross Validation (LOOCV) experiment.
Novel.R is for scoring drug synergy based on the all known synergistic drug combinations.
Para.R is for conducting the grid search method to determine all parameter combinations.
SDCLR.R is for implementing SDCLR algorithm via LOOCV experiment.
Sim_Rank.R is for implementing SimRankLN method to calculate protein similarity.
