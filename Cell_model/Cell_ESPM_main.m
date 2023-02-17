%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Simulate a lithium-ion cells  %%%%%%%%

% % Model description: Coupled electrochemical
% Electrochemical model: Enhanced Single Particle Model (ESPM)

% % Specify the following inputs:
% 1) simulation time: t_duration
% 2) current input: I_data
% 3) Number of radial discretization grids in ESPM: Nr
% 4) Number of cartesian discretization grids in ESPM: Nx_n , Nx_s, Nx_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

addpath('Cell_model')

%% SPECIFY: Input current profile
Crate = [1]; 
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
param.t_duration = 3600;  % simulation time in seconds
param.dt =  1;               % sampling time            
param.t_data = [0:param.dt:param.t_duration]'; % time vector


% applied current = Nc*nominal current value
  I_multiply = param.capacity*Crate;
  param.I_data = I_multiply*ones(length(param.t_data),1); % current vector


% Specify number of additional cycles beyond initial charge / discharge
% 'Run_module' script will concatenate alternating charge / discharge
% current profiles as needed to meet the input # of cycles
param.cycles = 0;

% SPECIFY: degree of variation in model parameters and initialize

% Specify initial SOC (SOC_ref = 1 means it as at 100% SOC)
param.SOC_ref = 1;                              % Average SOC


% SPECIFY: initial conditions of all cells
% Specify Ambient Temperature
param.T_amb = 298; %[K]


% Initial lithium concentration in solid phase for all cells corresponding
% to the cell SOC
[cs_initial, csn0, csp0] = conc_initial_sd(param.SOC_ref, param);

% Initial lithium concentration in electrolyte phase for all cells
param.ce_states = param.Nx_n + param.Nx_s + param.Nx_p;
ce_initial = param.ce0*ones(param.Nc*param.ce_states,1);


% Group all initial state variables
x_initial = [cs_initial; ce_initial];


[param_mod] = CellParametersFunction_MC_modstudy(param);

% NOTE: reduced number of outputs for memory management
[V_cell,T_c, T_s, soc_n, soc_p, I_c,param] = Cell_Run(x_initial,param_mod);

%% Plot Current for both positve and negative perturb.

% figure(1); hold on; grid on; set(gcf,'color','w');
% tlt = tiledlayout(2,4);
% tlt.TileSpacing = 'compact'; tlt.Padding = 'normal';  
% 
% nexttile([1 2])
% for i=1:length(Crate)
%     plot((1:length(V_cell))/3600, V_cell, 'LineWidth',1.5); hold on; grid on;
% end
% ylim([2.5 4.2])
% ylabel('Voltage [V]', 'interpreter','Latex')
% xlabel('Time [h]', 'interpreter','Latex')
% 
% yyaxis right
%     plot((1:length(soc_n))/3600, soc_n, '--', 'LineWidth',1.5); hold on; grid on;
% 
% ylabel('SOC [-]', 'interpreter','Latex')
% xlabel('Time [h]', 'interpreter','Latex')


% nexttile
% for i=1:length(Crate)
%     plot((1:t_cutoff(i))/3600, soc_bulk_p(1,1:t_cutoff(i),i), 'LineWidth',1.5); hold on; grid on;
% end
% ylabel('SOC [-]', 'interpreter','Latex')
% xlabel('Time [h]', 'interpreter','Latex')
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
