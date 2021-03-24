# MD-scoring
Script and associated data for MD-based scoring of HADDOCK clusters

This repository is part of a protocol based on molecular dynamics (MD) simulations would allow to distinguish native from non-native models to complement scoring functions used in docking. First, models for 25 protein-protein complexes were generated using HADDOCK. Next, MD simulations complemented with machine learning were used to discriminate between native and non-native complexes based on a combination of metrics reporting on the stability of the initial models. Here one can access used scripts, jupyter notebooks of machine learning models and extracted features which served as input for machine learning.


## MD trajectories
MD trajectories of 100 ns for 25 protein-protein complexes (8 HADDOCK models x 2 replicas, reference crystal structure in 4 replicas per model) are available at [zenodo](http://doi.org/10.5281/zenodo.4629895)

## Obtaining complex properties -> features as input for ML classifiers
All trajectories were analyzed by scripts in the `gmx_scripts` directory 

## Features
All properties were extracted into features .csv files. Feature files for different training and validation sets are in the `features` directory. 

## RF classifiers
* First, the most fitting(accurate) ML classifier was selected. This was done comparing 9 different classifiers from the Scikit-learn library. The jupyter notebook used for this purpose is `Classifier selection`.
* To optimize the most accurate classifier (Random Forest in our case) two runs of optimization were run (Random and Grid Search). Both of these procedures are shown in the `Parameters_search` notebook.
* To estimate the accuracy of the Random Forest classifier we first tested it the training set using K-fold cross-validation (100x) in `Cross_validaiton`.
* To test the classifier on an external validations set(s) the`External_validation` notebook was used.