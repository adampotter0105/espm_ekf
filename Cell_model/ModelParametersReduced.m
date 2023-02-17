%% Initialize Model Parameters (SI Units)
% Defining a structure named "params" that contains all parameter values
% BATTERY DETAILS (i.e. geometry, chemistry)
% FOR ELECTROCHEMICAL & THERMAL MODELS: 
% Cell Name: US18650VTC4
% Manufacturer: Sony
% Rated Capacity = 2 A-h
% Nominal Voltage = 3.7V / Maximum Voltage = 4.2V / Minimum Voltage = 2.5V
% Cathode Chemistry = NMC / Anode Chemistry = Graphite

% FOR AGING MODELS - SEI LAYER GROWTH & SEMI-EMPIRICAL: 
% Cell Name: ANR26650
% Manufacturer: A123
% Rated Capacity = 2.3 A-h
% Nominal Voltage = 3.3V / Maximum Voltage = 3.6V / Minimum Voltage = 2V
% Cathode Chemistry = LFP / Anode Chemistry = Graphite
param.Nc = 1; 
param.Rg = 8.314;                 % Gas constant
param.F = 96487;                  % Faraday's constant
param.alpha_cell = 0.5;           % Anode/Cathode transfer coefficient
param.ce0 = 1200;                 % Average electrolyte concentration [mol/m^3]
param.capacity = 2;               % Nominal cell capacity (Ah)   

param.c_n_max = 31080;            % Maximum anode concentration
param.c_p_max = 51830;            % Maximum cathode concentration

% param.Dsn_ref = 2.87e-14;         % Anode diffusion coefficient
% param.Dsp_ref = 4.85e-14;         % Cathode diffusion coefficient

%TEST - Anirudh's values
param.Dsn_ref = 1.4493e-14;         % Anode diffusion coefficient
param.Dsp_ref = 5.8811e-14;         % Cathode diffusion coefficient

param.Rs_n = 4.2642e-06;                % Anode particle radius
param.Rs_p = 5.1592e-07;                % Cathode particle radius
param.A = 0.1048;                  % Area
param.Ln = 4.0000e-05;                 % Anode thickness
param.Lp = 3.6550e-05;                 % Cathode thickness
param.epsilon_n = 0.6187;          % Anode solid phase volume fraction
param.epsilon_p = 0.5800;          % Cathode solid phase volume fraction
param.theta100_n = 0.928;         % Anode stoichiometry at 100% SOC
param.theta100_p = 0.3504;        % Cathode stoichiometry at 100% SOC

% param.kn_ref = 3.48e-10;          % Anode reaction rate constant
% param.kp_ref = 4.164e-10;         % Cathode reaction rate constant

%Test - Anirudh's values
param.kn_ref = 1.0e-10;          % Anode reaction rate constant
param.kp_ref = 1.0e-10;         % Cathode reaction rate constant

param.theta0_n = 0.002;           % Anode stoichiometry at 0% SOC
param.theta0_p = 0.9986;          % Cathode stoichiometry at 0% SOC

%Introduce stochastic variation to Lumped Contact Resistance
 param.R_l_ref = 0.032;            % Lumped resistance
% param.R_l_ref = 1; 
% param.R_l_deviation = param.R_l_ref/20; % Lumped resistance Stochastic Deviation
% param.R_l = param.R_l_ref*ones(param.Nc,1) + (param.R_l_deviation)*randn(param.Nc,1);
param.R_l = param.R_l_ref;

% ESPM paramters - TRY USING THIS - TEST
% param.Ln = 40e-6;                   % Anode Thickness [m]
param.Ls = 2.5000e-05;                   % Separator Thickness [m]
% param.Lp = 36.55e-6;                % Cathode Thickness [m]
% param.Rs_n = 4.2642e-6;
% param.Rs_p = 5.1592e-7;
% param.A = 0.1048;
% param.epsilon_n = 0.6187;
% param.epsilon_p = 0.5800;
% param.c_n_max = 31080;            % Maximum anode concentration
% param.c_p_max = 51830;            % Maximum cathode concentration
% ESPM parameters - TRY USING THIS - TEST


% Specific interfacial (electroactive) surface area for anode and cathode
param.a_sn = 3*param.epsilon_n/param.Rs_n;   
param.a_sp = 3*param.epsilon_p/param.Rs_p; 

% Electrolyte porosity
param.eps_filler_n = 0.038;
param.eps_filler_p = 0.12; %Value from ESPM code
% param.eps_filler_p = 0.012; %Original value in module code
% param.eps_filler_s = 0.4;
param.eps_el_n_ref = 1-param.epsilon_n-param.eps_filler_n;
param.eps_el_p = 1-param.epsilon_p-param.eps_filler_p;

% param.eps_filler_s = 0.4;
% param.eps_el_s = 1-param.eps_filler_s; %NOTE: VERIFY THIS

param.eps_el_s = 0.4; %NOTE: Switched this with value above

%Additional Electrolyte Parameters
%NOTE: NEED TO VERIFY IF THIS SEPARATOR THICKNESS WAS ALSO IDENTIFIED W/
%ELECTRODE THICKNESSES SPECIFIED ABOVE
param.t0 = 0.369;                   % Transference Number
param.brugg = 1.5;                  % Bruggeman coefficient
param.Ln_el = 40e-6;                   % Anode Thickness [m]
param.Ls = 25e-6;                   % Separator Thickness [m]
param.Lp_el = 36.55e-6;                % Cathode Thickness [m]

%TEST
% param.Ln = 40e-6;                   % Anode Thickness [m]
% param.Ls = 25e-6;                   % Separator Thickness [m]
% param.Lp = 36.55e-6;                % Cathode Thickness [m]

param.brugg_n = param.brugg;
param.brugg_s = param.brugg;
param.brugg_p = param.brugg;

% Derived parameters for electrolyte phase discretization
% Grid interval size for anode, cathode and separator regions

%For finite difference method:
% param.delta_n = param.Ln/(param.Nx_n -1);
% param.delta_s = param.Ls/(param.Nx_s -1);
% param.delta_p = param.Lp/(param.Nx_p -1);

%For finite volume method:
param.delta_n = param.Ln/(param.Nx_n);
param.delta_s = param.Ls/(param.Nx_s);
param.delta_p = param.Lp/(param.Nx_p);


%FOR SENSITIVITY ANALYSIS:
param.De_scaling = 1; %Electrolyte Diffusivity Scaling Factor
param.K_el_scaling = 1; %Electrolyte Conductivity Scaling Factor
