%% Specify Cell Cycling
% SPECIFY: Finite Difference / Volume Discretization Parameters
param.Nr = 4;      % Number of radial discretization grids in ESPM
param.Nx_n = 4;    % Number of cartesian discretization grids in ESPM
param.Nx_s = 4;    % Number of cartesian discretization grids in ESPM
param.Nx_p = 4;    % Number of cartesian discretization grids in ESPM
param.Nsei = 4;    % SEI Layer Discretization


% For a custom input current profile 
% Current convention: Positive for discharge, negative for charge
param.t_duration = DT;  % simulation time in seconds
param.dt =  DT;               % sampling time            
param.t_data = [0:param.dt:param.t_duration]'; % time vector

% Specify number of additional cycles beyond initial charge / discharge
% 'Run_module' script will concatenate alternating charge / discharge
% current profiles as needed to meet the input # of cycles
param.cycles = 0;

if usePrevParams == false
    param.T_amb = 298; % [K]
    % Model Parameters for 2Ah Cell
    run ModelParametersReduced.m
    param = CellParametersFunction_MC_modstudy(param);
    param.I_data = ones(length(param.t_data),1);

    
    %% INITIAL CONDITIONS
    % Specify initial SOC (SOC_ref = 1 means it as at 100% SOC)
    param.SOC_ref = 1;                              % Average SOC

end

%% Additional parameters needed for this ESPM model
%FOR SENSITIVITY ANALYSIS:
param.De_scaling = 1; %Electrolyte Diffusivity Scaling Factor
param.K_el_scaling = 1; %Electrolyte Conductivity Scaling Factor

%Kinetic & Transport Parameter Updates
param.De_scaling = param.De_scaling*ones([param.Nc,1]); 
param.K_el_scaling = param.K_el_scaling*ones([param.Nc,1]);
param.kn_ref = param.kn_ref*ones([param.Nc,1]);          
param.kp_ref = param.kp_ref*ones([param.Nc,1]);
param.Dsn_ref = param.Dsn_ref*ones([param.Nc,1]); 
param.Dsp_ref = param.Dsp_ref*ones([param.Nc,1]);

param.Tref = 25+273; 

% Activation energy for temperature Arrhenius dependency [J/mol]
param.Ea_Dsp = 25000;
param.Ea_Dsn = 50000;
param.Ea_kp = 30000;
param.Ea_kn = 30000;
%% End of additional parameters

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