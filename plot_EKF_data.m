%% Plot EKF Results
dV0 = li2voltage(estimateStates(:,1)) - li2voltage(trueStates(:,1))

% Estimate Voltage using state estimates
V_estimate(1,numSteps) = 0;
V_real(1,numSteps) = 0;
SOC_n_estimate(1,numSteps) = 0;
SOC_p_estimate(1,numSteps) = 0;
SOC_n_real(1,numSteps) = 0;
SOC_p_real(1,numSteps) = 0;
for i = 1:numSteps
    [V_estimate(i), SOC_n_estimate(i), SOC_p_estimate(i)] = li2voltage(estimateStates(:,i));
    [V_real(i), SOC_n_real(i), SOC_p_real(i)] = li2voltage(trueStates(:,i));
end

figure(1)
plot(tspan,V_estimate)
hold on
plot(tspan,V_real)
plot(tspan,measurements)
legend(["EKF Voltage", "Real Voltage","Measured Voltage"])
ylim([3.5 4.5])
xlabel("Time (s)")
ylabel("Voltage (V)")
hold off
improvePlot

figure(2)
plot(tspan,SOC_n_estimate*1e2)
hold on
plot(tspan,SOC_p_estimate*1e2)
plot(tspan,SOC_n_real*1e2)
plot(tspan,SOC_p_real*1e2)
legend(["EKF SOC Anode", "EKF SOC Cathode", "Real SOC Anode", "Real SOC Cathode"])
xlabel("Time (s)")
ylabel("SOC (%)")
hold off
improvePlot

figure(3)
yyaxis left
plot(tspan,sum(estimateStates(1:param.Nr-1,:),1))
hold on
plot(tspan,sum(trueStates(1:param.Nr-1,:),1))
yyaxis right
plot(tspan,sum(estimateStates(param.Nr:2*(param.Nr-1),:),1))
plot(tspan,sum(trueStates(param.Nr:2*(param.Nr-1),:),1))
legend(["EKF Li Anode", "Real Li Anode", "EKF Li Cathode", "Real Li Cathode"])
xlabel("Time (s)")
ylabel("Li Concentration")
hold off
improvePlot

function [V_cell, soc_n, soc_p] =li2voltage(x_out)
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
[soc_bulk_n,soc_bulk_p] = soc_calculation(cs_n,cs_p,param);
soc_n = mean(soc_bulk_n);
soc_p = mean(soc_bulk_p);
end