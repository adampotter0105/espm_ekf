function [param] = ModelParametersFunction(param)

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
param.capacity = 2;  % Nominal cell capacity (Ah) - Reference Value   

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
%param.R_l_ref = 1;
% param.R_l_deviation = param.R_l_ref/20; % Lumped resistance Stochastic Deviation
% param.R_l = param.R_l_ref*ones(param.Nc,1) + (param.R_l_deviation)*randn(param.Nc,1);

% ESPM paramters - TRY USING THIS - TEST
% param.Ln = 40e-6;                   % Anode Thickness [m]
param.Ls = 2.5000e-05;                   % Separator Thickness [m]



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


% Derived parameters for solid phase discretization
% Grid interval size for anode and cathode particle
param.delta_xn = param.Rs_n/(param.Nr-1);    
param.delta_xp = param.Rs_p/(param.Nr-1);


%% Restructure all parameters to be perturbed within pack
%Primary Parameters
param.epsilon_n = param.epsilon_n*ones([param.Nc,1]);
param.epsilon_p = param.epsilon_p*ones([param.Nc,1]);
param.Ln = param.Ln*ones([param.Nc,1]);
param.Lp = param.Lp*ones([param.Nc,1]);
param.Ls = param.Ls*ones([param.Nc,1]);
param.A = param.A*ones([param.Nc,1]);
%param.A = [0.1048; 0.1058; 0.1068; 0.1078; 0.1088; 0.1098];
param.R_l = param.R_l*ones([param.Nc,1]);
%param.R_l = [0.03; 0.0305; 0.031; 0.0315;0.032; 0.0325];
param.eps_filler_n = param.eps_filler_n*ones([param.Nc,1]);
param.eps_filler_p = param.eps_filler_p*ones([param.Nc,1]);
param.eps_el_s = param.eps_el_s*ones([param.Nc,1]);
param.ce0 = param.ce0*ones([param.Nc,1]);

%Downstream Parameters
param.a_sn = param.a_sn*ones([param.Nc,1]);
param.a_sp = param.a_sp*ones([param.Nc,1]);
param.delta_n = param.delta_n*ones([param.Nc,1]);
param.delta_p = param.delta_p*ones([param.Nc,1]);
param.delta_s = param.delta_s*ones([param.Nc,1]);
param.eps_el_n_ref = param.eps_el_n_ref*ones([param.Nc,1]);
param.eps_el_p = param.eps_el_p*ones([param.Nc,1]);
param.capacity = param.capacity*ones([param.Nc,1]);  

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

end