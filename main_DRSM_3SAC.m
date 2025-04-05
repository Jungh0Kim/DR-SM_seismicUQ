%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jungho Kim, UC Berkeley
% junghokim@berkeley.edu
%
% This code is part of the following publication:
% Kim, J., & Wang, Z. (2025). Uncertainty quantification for seismic response using
% dimensionality reduction‐based stochastic simulator. Earthquake Engineering & Structural Dynamics, 54(2), 471-490.
% https://doi.org/10.1002/eqe.4265
%
% This script presents stochastic simulator predictions for 3SAC 
% building responses subjected to recorded ground motions.
% Inputs = [Hazard parameters, Structural parameters, Ground motion time histories]
% Outputs = [IDRs]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; tic;
set(0,'DefaultFigureColor',[1 1 1])
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',14);

%% Add paths

addpath('drtoolbox\corefunc')
addpath('drtoolbox\gui')
addpath('drtoolbox\techniques')
addpath('netlab_GMR')

%% Load training and testing datasets

% When generating training data, using latin hypercube sampling
% with sample decorrelation is highly recommended
% e.g., lhs_train = lhsdesign(n_train, n_rvs, 'criterion','correlation');
% x_train = norminv(lhs_train);

load('Train_Test_dataset_seismic_3SAC.mat')

%%% Training set (700 data)
x_H_train = x_H_train(:,1:2);    % Hazard parameters: [M, R]
x_S_train;                       % Structural parameters: [xi, E, fy_b, fy_c, Ha_b, Ha_c]
x_GM_train_cell;                 % Ground motion records at discretized time steps
y_train;                         % Building EDPs: [IDR1, IDR2, IDR3]

%%% Testing set (1000 data)
x_H_test = x_H_test(:,1:2);
x_S_test;
x_GM_test_cell;
y_test;                          % True responses to be predicted

[n_train, Y_size] = size(y_train);
n_test = size(y_test,1);

% Here we saved uncorrelated x_S in an uncorrelated standard Gaussian space.
% You can transform data using probabilistic transformation (see Transform.m)

%% Set parameters

rng(14)

iter_n = 1300;         % Length of random sequence
cutoff_n = 300;        % Burn-in
iM = 7;                % Number of components in GMM for condi. distr. model

% Physics-based mapping (H_p in Eq. (5))
PDR_func_name = 'PDR_GM';
PDR_func = str2func(PDR_func_name);

% 'Conventional' DR mapping (H_DR in Eq. (5))
DR_method = 'PCA';

%% Augmented input-output set

% Training inputs
PDR_train_vars = PDR_func(x_GM_train_cell);
x_train_prime = [x_H_train x_S_train PDR_train_vars];

% Augmented training dataset
Tot_train = [x_train_prime, y_train];

%% Construct dimension reduction mapping and conditional distribution model

n_dim_redc = 26;     % This should be selected based on Algorithm 1
% From our experience, the selection of n_dim_redc is not critical
% since the core of the proposed method is to "extract" a stochastic
% surrogate model from the low-dimensional representation of the
% input-output space.
% You may choose a suitable "dimensionality reduction algorithm" and
% "conditional distribution model", and tune their parameters and the
% reduced dimensionality accordingly.

% Input-output mapping
[Psi_z, Psi_z_mapping] = compute_mapping(Tot_train, DR_method, n_dim_redc);

% Conditional distribution modeling
sMixR = GMR(Psi_z, y_train, iM, 'full');

toc;

fprintf('\n Training finished. \n')

%% Extract stochastic surrogate modeling for testing dataset

% Testing inputs
PDR_test_vars = PDR_func(x_GM_test_cell);
x_test_prime = [x_H_test x_S_test PDR_test_vars];

% Simulate random sequence using transition kernel in Eq.(4).
% This fixed-point iteration is performed independently for each test input
% to extract a stationary distribution from the transition kernel.
% Parallel computing (e.g., parfor) can be used to accelerate this process.
mu_y_train = mean(y_train);
y_random_sequence = zeros(iter_n+1, Y_size, n_test);
for i = 1:n_test
    x_test_i = x_test_prime(i,:);
    y_start = mu_y_train; % starting point
    
    y_sequence_i = DRSM_sequence(y_start, x_test_i, Psi_z_mapping, sMixR, iter_n);
    y_random_sequence(:,:,i) = y_sequence_i;

    disp(i);
end

% Remove burn-in
y_random_sequence = y_random_sequence(cutoff_n+1:end, :, :);

% Final prediction mean and std.
y_pred_mean = mean(y_random_sequence, 1);
y_pred_mean = reshape(y_pred_mean, Y_size, n_test);
y_pred_mean = y_pred_mean';
y_pred_std = std(y_random_sequence);
y_pred_std = reshape(y_pred_std, Y_size, n_test);
y_pred_std = y_pred_std';

toc;

fprintf('\n Prediction finished. \n')

%% Prediction error

% Relative mean squared error (Eq. (14))
for kk=1:Y_size
    rMSE(:,kk) = mean((y_test(:,kk) - y_pred_mean(:,kk)).^2)/std(y_test(:,kk)).^2;
end

rMSE

%% Plot prediction

kp = 2;
shaded = [0.7 0.7 0.7];
str_set = {'IDR_1','IDR_2','IDR_3'};

% Plot prediction mean/+-2std based on sorted testing data
figure()
for kk=1:Y_size
    [~,y_dix] = sort(y_pred_mean(:,kk));
    y_test_p = y_test(y_dix,kk);
    y_pred_p = y_pred_mean(y_dix,kk);
    y_pred_std_p = y_pred_std(y_dix,kk);
    y_pred_CI_p = [y_pred_p + kp*y_pred_std_p; flip(y_pred_p - kp*y_pred_std_p,1)];

    subplot(1, 3, kk)
    fill([(1:length(y_pred_p))'; flip((1:length(y_pred_p))',1)], y_pred_CI_p,...
        shaded, 'EdgeColor', shaded, 'FaceAlpha', 0.9);  hold on
    plot(1:length(y_pred_p), y_pred_p,'linewidth',2,'color','k');
    plot(1:length(y_test_p), y_test_p,'b.','markersize',5)
    xlabel('Test data sorted in ascending order')
    ylabel(strcat('$Y_{',num2str(kk),'} \,(=',str_set{1,kk},')$'))
    title("Error:" + num2str(rMSE(kk)))
    if kk==1
        lgnd = legend('Mean 2std',' Mean',' True y','Location','northwest');
        set(lgnd,'FontSize',12.2,'NumColumns',1);
    end
    hold off
end
set(gcf,'unit','centimeters','position',[0 0 35 10]);
