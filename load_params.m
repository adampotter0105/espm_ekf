% Time to Simulate in EKF
DT = 5; % Time step for the EKS, model will integrate over each of these steps

addpath('Cell_model')

%% Specify Cell Cycling
% SPECIFY: Finite Difference / Volume Discretization Parameters
param.Nr = 10;      % Number of radial discretization grids in ESPM
param.Nx_n = 10;    % Number of cartesian discretization grids in ESPM
param.Nx_s = 10;    % Number of cartesian discretization grids in ESPM
param.Nx_p = 10;    % Number of cartesian discretization grids in ESPM
param.Nsei = 10;    % SEI Layer Discretization

% Initialize all model parameters
run ModelParameters.m

% For a custom input current profile 
% Current convention: Positive for discharge, negative for charge
param.t_duration = DT;  % simulation time in seconds
param.dt =  DT;               % sampling time            
param.t_data = [0:param.dt:param.t_duration]'; % time vector

% Specify number of additional cycles beyond initial charge / discharge
% 'Run_module' script will concatenate alternating charge / discharge
% current profiles as needed to meet the input # of cycles
param.cycles = 0;

param.T_amb = 298; % [K]
% SPECIFY: degree of variation in model parameters and initialize

param = CellParametersFunction_MC_modstudy(param);
param.I_data = ones(length(param.t_data),1);

%% INITIAL CONDITIONS
% Specify initial SOC (SOC_ref = 1 means it as at 100% SOC)
param.SOC_ref = 1;                              % Average SOC


% Initial lithium concentration in solid phase for all cells corresponding
% to the cell SOC
[cs_initial, csn0, csp0] = conc_initial_sd(param.SOC_ref, param);

% Initial lithium concentration in electrolyte phase for all cells
param.ce_states = param.Nx_n + param.Nx_s + param.Nx_p;
ce_initial = param.ce0*ones(param.Nc*param.ce_states,1);


% Group all initial state variables
x_initial = [cs_initial; ce_initial];
n_states = max(size(x_initial));

%% Generate matrices for solid phase, electrolyte phase, thermal, and aging

% Solid phase discretization matrices
[param.A_sd, param.B_sd_n, param.B_sd_p] = matrices_solidphase(param);

% Electrolyte phase discretization matrices
[param.A_en_d, param.A_es_d, param.A_ep_d, param.B_e_d] = matrices_elecphase(param);

%% Separate compiled coefficient matrices (i.e. Define A1_en, A2_en, A3_en, etc.) 
run separate_elec_matrices.m
