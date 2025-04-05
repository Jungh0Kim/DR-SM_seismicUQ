DRSM_seismicUQ v1.0
========

Matlab code for Dimensionality reduction as a surrogate model (DRSM) for seismic response UQ, from:

Kim, J., & Wang, Z. (2025). Uncertainty quantification for seismic response using dimensionality reduction‐based stochastic simulator. Earthquake Engineering & Structural Dynamics, 54(2), 471-490.
https://doi.org/10.1002/eqe.4265.

Included toolboxes:
- drtoolbox - https://lvdmaaten.github.io/drtoolbox/
- Netlab - http://www1.aston.ac.uk/ncrg/ (Minor: https://uk.mathworks.com/matlabcentral/fileexchange/2654-netlab)

Notes:
 - Some functions in the Netlab toolbox have been modified to support the conditional distribution model used in DRSM.
 - This script implements stochastic simulator for predicting responses of the 3 story SAC building subjected to recorded ground motions.
 
 How to Run:
 - Download the dataset: https://www.dropbox.com/scl/fi/gimufs0r2yplz7d1rzg0i/Train_Test_dataset_seismic_3SAC.mat?rlkey=0fc8mma2h3aft76bqnfj0tzmm&st=77sfdpx5&dl=0
 - Run main_DRSM_3SAC.m to reproduce the results
