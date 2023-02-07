clear
rng(2022);    % For repeatable results
run load_params.m % Run the ECM Model and System Params Script
simTime = 100; % seconds
tspan = 0:DT:simTime;
numSteps = length(tspan);

trueInitialState = x_initial; 
global param

% Control Param:  applied current = Nc*nominal current value
target_C_rate = 5;
param.I_data = target_C_rate*param.capacity*ones(length(param.t_data),1); % current vector

% Noise and Covariance Matrices
% TODO: Find these numbers for real
processNoise = diag(100*ones(1,n_states)); % Process noise matrix
measureNoise = 1.0006e-08; % Measurement noise matrix (measurement is 1x1).

% Using EPSM Model as Truth
trueStates = NaN(n_states,numSteps);
trueStates(:,1) = x_initial;
measurements = NaN(1,numSteps);

fprintf("Simulating real battery \n")
for i = 2:length(tspan)
    if i ~= 1
        trueStates(:,i) = ESPMmodel(trueStates(:,i-1),DT) + sqrt(processNoise)*randn(n_states,1);  
    end
    measurements(:,i) = li2voltage(trueStates(:,i)) + sqrt(measureNoise)*randn(1,1);
end

% Making Starting State wrong
x0_error_factor = 0.05;
x_initial_noisy = x_initial;
x0_err = x_initial_noisy(1:param.Nr-1)*x0_error_factor; % Move a percentage of anode lithium to cathode 
x_initial_noisy(1:param.Nr-1) = x_initial_noisy(1:param.Nr-1) - x0_err;
x_initial_noisy(param.Nr:2*(param.Nr-1)) = x_initial_noisy(param.Nr:2*(param.Nr-1)) + x0_err*param.Rs_n/param.Rs_p; 

% Create the EKF
fprintf("Creating EKF \n")
filter = trackingEKF(State=x_initial_noisy, StateTransitionFcn=@ESPMmodel,ProcessNoise=processNoise, ...
    MeasurementFcn=@li2voltage,MeasurementNoise=measureNoise);
% filter = trackingEKF(State=x_initial_noisy,StateTransitionFcn=@ESPMmodel,MeasurementFcn=@li2voltage);

% Run the EKF
estimateStates = NaN(size(trueStates)); % preallocate state estimates
estimateStates(:,1) = filter.State;

runtimes = zeros(1,length(tspan)-1);
for i=2:length(tspan)
    tic
    predict(filter,DT);
    estimateStates(:,i) = correct(filter,measurements(:,i));
    runtimes(i) = toc;
    fprintf("EKF step: %d \n",i-1)
end
fprintf("Average Runtime for EKF Step: %d seconds", mean(runtimes))

%% Helper Functions

% ESPM Model State Integrator
function x_out = ESPMmodel(x_in, ~)  % param
global param
tspan = param.t_data;

%% Solve ODEs 
reltol=5.0e-06; abstol=5.0e-06;
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
    %V_oc = ocp_p - ocp_n;
end


%% Calculate SOC
%[soc_bulk_n,soc_bulk_p] = soc_calculation(cs_n,cs_p,param);
end