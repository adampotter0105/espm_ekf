clear
% Define Parameters for EKF Simulation
global param
rng(2022);    % For repeatable results
DT = 1; % Seconds for EKF to run
simTime = 90; % seconds
tspan = 0:DT:simTime;
numSteps = length(tspan);
param.t_data = tspan;

% Whether to use Params from the casadi model
usePrevParams = true;

% Battery Initial Conditions
T_amb = 25; % [oC]
SOC_IC = 0.8; % Initial SOC
Q_IC = 4.8769; % Nominal capacity [Ah]
Lsei_IC = 5e-9; % Intial width of sei [m]

% Set Current Profile  
target_C_rate = 1;
I_data = target_C_rate*param.capacity*ones(numSteps,1) + 1*sin(tspan'/5); % current vector
param.I_data = I_data;

% Calculate columb-counted SOC at each timestep
SOC_cc = zeros(1,length(param.t_data));
SOC_cc(1) = SOC_IC;
for i = 2:length(numSteps)
    SOC_cc(i) = SOC(i-1) - param.I_data(i-1)/(3600*Q_IC);
end
param.SOC_cc = SOC_cc;

% Simulate Full Model For Real Data
cd('ESPM_casadi\cell_model\')
fprintf("Simulating advanced real battery \n")
[V_real, ~, ~, ~, SOC_n_real, SOC_p_real, cs_n_real, cs_p_real,...
          ~, ~,~, ~,~, ~, ~, ~, ~, ~,...
          ~, ~, ~, ~, ~, ~, param] = ESPM_sim(0, 0, param.t_data, param.I_data, ...
          SOC_cc, SOC_IC,Q_IC, Lsei_IC, T_amb) ;
% run ModelParameters.m
cd('../../Cell_model/')
% Generate Parameters for Reduced Model
run load_params.m % Generate and System Params and Associated Matrices

% Noise and Covariance Matrices
processNoise = diag(100*ones(1,n_states)); % Process noise matrix  TODO: Find these numbers for real
measureNoise = 1.0006e-08; % Measurement noise matrix (measurement is 1x1).

% Simulate same model next to it
[cs_initial, ~, ~] = conc_initial_sd(param.SOC_ref, param);
x_initial = [cs_initial; ce_initial];
trueStates = NaN(n_states,numSteps);
trueStates(:,1) = x_initial;
measurements = NaN(1,numSteps);
measurements(1) = li2voltage(trueStates(:,1)) + sqrt(measureNoise)*randn(1,1);
fprintf("Simulating same-model battery \n")
for i = 2:length(tspan)
    param.t_data = tspan(i-1:i);
    param.I_data = I_data(i-1:i);
    trueStates(:,i) = ESPMmodel(trueStates(:,i-1),DT) + sqrt(processNoise)*randn(n_states,1);  
    measurements(:,i) = li2voltage(trueStates(:,i)) + sqrt(measureNoise)*randn(1,1);
end

% Making Starting State Intentially Wrong    
x0_error_factor = -0.05; % percent error below real SOC
[cs_initial_noisy, ~, ~] = conc_initial_sd(param.SOC_ref*(1+x0_error_factor), param);
x_initial_noisy = [cs_initial_noisy; ce_initial];
n_states = size(x_initial_noisy,1);

%% Create the EKF
fprintf("Creating EKF \n")
ekf = trackingEKF(State=x_initial_noisy, StateTransitionFcn=@ESPMmodel,ProcessNoise=processNoise, ...
    MeasurementFcn=@li2voltage,MeasurementNoise=measureNoise);

%% Run the EKF
estimateStates_ekf = NaN(n_states,numSteps); % preallocate state estimates
estimateStates_ekf(:,1) = ekf.State;
runtimes1 = zeros(1,length(tspan)-1);

for i=2:length(tspan)
    param.t_data = tspan(i-1:i);
    param.I_data = I_data(i-1:i);
    tic
    predict(ekf,DT);
    estimateStates_ekf(:,i) = correct(ekf,measurements(i));
    runtimes1(i) = toc;
    fprintf("EKF step: %d \n",i-1)
end
fprintf("EKF Average Runtime for EKF Step: %d seconds \n", mean(runtimes1))

cd('../')

%% Reduced Model Functions
function x_out = ESPMmodel(x_in, dt)  % param
    global param
    tspan = param.t_data;
    
    %% Solve ODEs 
    reltol=5.0e-09; abstol=5.0e-07;
    event_formatted = @(t,x) physical_event_function(t,x,param); 
    options=odeset('RelTol',reltol,'AbsTol',abstol,'Events',event_formatted);
    [t_out, x_out,te,xe,ie] = ode15s(@(t_out, x_out) Cell_ode(t_out, x_out, param), tspan, x_in, options);
    
    %Transpose state matrix into row form to match established data structure 
    x_out = real(x_out'); %Keep the result real
    x_out = x_out(:,end);
end

% Lithium Concentration to Voltage Estimator
function V_cell =li2voltage(x_out)
    global param
    %% Separate electrochemical, thermal & aging state variables from x_out matrix
    %Define solid concentrations
    cs = x_out(1:(param.Nr-1)*param.Nc*2,:);            %All solid concentrations
    cs_n = cs(1:(param.Nr-1)*param.Nc,:);               %Anode Concentrations
    cs_p = cs((param.Nr-1)*param.Nc+1:end,:);           %Cathode Concentrations
    
    index_cs = (param.Nr-1)*param.Nc*2;                 %Index for final solid concentration
    index_ce = index_cs + param.Nc*param.ce_states;     %Index for final electrolyte concentration
    
    %Define electrolyte concentrations
    ce = x_out(index_cs+1:index_ce,:);                  %All electrolyte concentrations
    ce_n = ce(1:param.Nx_n*param.Nc,:);                 %Negative Electrolyte Region
    ce_s = ce(param.Nx_n*param.Nc+1:(param.Nx_n + param.Nx_s)*param.Nc,:); %Separator Region
    ce_p = ce((param.Nx_n + param.Nx_s)*param.Nc+1:param.ce_states*param.Nc,:); %Positive Region
    
    %% Temperature dependent transport and kinetics for each cell
    
    for j = 1:length(cs(1,:))
        I_cell(1,j) = param.I_data(j);  
        T_core(1,j) = param.T_amb;
        T_surf(1,j) = param.T_amb;
        % Surface concentration -> surface stoichiometry
        for i = 1:param.Nc
            theta_surf_n(i,j) = cs_n(i*(param.Nr-1),j)/param.c_n_max;
            theta_surf_p(i,j) = cs_p(i*(param.Nr-1),j)/param.c_p_max;
        end
        
        % Temperature dependence
        [Dsn(:,j), Dsp(:,j), kn(:,j), kp(:,j)] = arrhenius_temp(param, T_core(:,j));
    
        
        % Define electrolyte boundary concentrations for each cell
        for i = 1:param.Nc
            index1n = (i-1)*param.Nx_n+1;
            index2n = i*param.Nx_n;
            index1s = (i-1)*param.Nx_s+1;
            index2s = i*param.Nx_s;
            index1p = (i-1)*param.Nx_p+1;
            index2p = i*param.Nx_p;
          
            %Average Electrolyte Concentrations for conc. dependent parameters
            ce_n_avg(i,j) = mean(ce_n(index1n:index2n,j));
            ce_s_avg(i,j) = mean(ce_s(index1s:index2s,j));
            ce_p_avg(i,j) = mean(ce_p(index1p:index2p,j));
            %ce_all(Row = 'Grid point', Column = 'Time point', 3rd Dimension = 'Cell #')
            ce_all(:,j,i) = [ce_n(index1n:index2n,j); ce_s(index1s:index2s,j); ce_p(index1p:index2p,j)];
            
            %ce_all_avg(Row = 'Cell #', Column = 'Time point')
            ce_all_avg(i,j) = mean(ce_all(:,j,i));
            
        end
        eps_el_n = param.eps_el_n_ref;
        % Electrolyte Conductivity
        [K_el_eff_n(:,j),K_el_eff_s(:,j),K_el_eff_p(:,j)] = conductivity_update(ce_n_avg(:,j)...
            ,ce_s_avg(:,j),ce_p_avg(:,j),eps_el_n,T_core(:,j),param);
        
        % Open circuit potential and overpotential
        ocp_p(:,j) = U_p(theta_surf_p(:,j));
        ocp_n(:,j) = U_n(theta_surf_n(:,j));
        
       
        % NOTE TO SELF: replace all param.I_data with I_cell
        % Using I_cell, update current-dependent electrochemical parameters
        eta_p(:,j) = eta_cathode(theta_surf_p(:,j), ce_p_avg(:,j), T_core(:,j), I_cell(:,j), param, kp(:,j));
        eta_n(:,j) = eta_anode(theta_surf_n(:,j), ce_n_avg(:,j), T_core(:,j), I_cell(:,j), param, kn(:,j));
    
    
        % Cell voltages
        [phi_e(:,j),R_el(:,j)] = electrolyte_potential(ce_all_avg(:,j), ce_all(:,j,:),K_el_eff_n(:,j),...
            K_el_eff_s(:,j), K_el_eff_p(:,j),T_core(:,j), I_cell(:,j), param);
        V_cell(:,j) = V_calculation(ocp_p(:,j),ocp_n(:,j),eta_p(:,j),eta_n(:,j),...
            phi_e(:,j),I_cell(:,j),param);
        V_oc = ocp_p - ocp_n;
end


%% Calculate SOC
%[soc_bulk_n,soc_bulk_p] = soc_calculation(cs_n,cs_p,param);
end
