%ESPM Model State Integrator
function dx_out = ESPMmodel(x_initial, ~)  % param
global param
%% Solve ODEs 
tspan = param.t_data;
param.TIMEOLD = datetime('now','Format','HHmmss');

import casadi.* %call new solver
% %NB: Set to 5e-12 for HPPC; 5e-5 for capacity and UDDS
% reltol=5.0e-8; abstol=reltol*0.001;
% 
% event_formatted = @(t,x) physical_event_function_tk(t,x,param); 
% options=odeset('RelTol',reltol,'AbsTol',abstol,'Events',event_formatted);
% [t_out, x_out,te,xe,ie] = ode15s(@(t_out, x_out) Module_ode(t_out, x_out, param), tspan, x_initial, options);
% 
% % Transpose state matrix into row form to match established data structure 
% x_out = x_out';

%================Define ODE using SUNDIALS
if var(param.I_data)>1e-4
    totN=size(x_initial,1);
    sundials_u = SX.sym('sundials_u',totN);
    input_crt=SX.sym('input_crt');
    fun_ode=Module_ode_sundials(sundials_u,input_crt, param);
    dae = struct('x',sundials_u,'p',input_crt,'ode',[fun_ode],'alg',[]);
    delta_time_ref=mode(param.t_data(end)-param.t_data(end-1));
    
    opts= struct('tf',delta_time_ref,'abstol',5.0e-8*0.001,'reltol',5.0e-8,'linear_solver','csparse',...
    'newton_scheme','gmres','show_eval_warnings',false);  %
    % opts= struct('tf',delta_time_ref,'abstol',5.0e-8*0.001,'reltol',5.0e-8,...
    % 'newton_scheme','gmres','show_eval_warnings',false);  %
    F1 = integrator('F', 'idas', dae,opts);
    u0=x_initial';
    
    
    cc=1;
    t_out=0;
    for kk=0:1:size(param.t_data,2)-2
        %     kk
        delta_time=param.t_data(kk+2)-param.t_data(kk+1);
        t_out=[t_out;t_out(end)+delta_time];
        if delta_time~=delta_time_ref
            opts= struct('tf',delta_time,'abstol',5.0e-8*0.001,'reltol',5.0e-8,'linear_solver','csparse',...
            'newton_scheme','gmres','show_eval_warnings',false);  %
            % opts= struct('tf',delta_time,'abstol',5.0e-8*0.001,'reltol',5.0e-8,'linear_solver','csparse',...
            % 'show_eval_warnings',false);  %
            F1 = integrator('F', 'idas', dae,opts);
        end
        iapp=param.I_data(kk+1);
        sol = F1('x0',u0,'p',iapp);
        u0=full(sol.xf);
        
        rex(:,cc)=full(sol.xf);
        cc=cc+1;
        
        
        %==========Monitor concentraion constraion not used for identificaiton
        if cc>10    
            monitorxnow=rex(:,cc-2);
            monitorxpre=rex(:,cc-3);
            
            value1now = min(monitorxnow(1:param.Nc*(param.Nr-1),:)) - (param.theta0_n*param.c_n_max+10);
            value1pre = min(monitorxpre(1:param.Nc*(param.Nr-1),:)) - (param.theta0_n*param.c_n_max+10);
            
            value2now = max(monitorxnow(1:param.Nc*(param.Nr-1))) - (param.theta100_n*param.c_n_max-10);
            value2pre = max(monitorxpre(1:param.Nc*(param.Nr-1))) - (param.theta100_n*param.c_n_max-10);
            
            value3now = min(monitorxnow(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - (param.theta100_p*param.c_p_max+10);
            value3pre = min(monitorxpre(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1))) - (param.theta100_p*param.c_p_max+10);
            
            value4now = max(monitorxnow(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1),:)) - (param.theta0_p*param.c_p_max-10);
            value4pre = max(monitorxpre(param.Nc*(param.Nr-1)+1:2*param.Nc*(param.Nr-1),:)) - (param.theta0_p*param.c_p_max-10);
            
        end     
    end
    
    
    x_out=rex;
    x_out=[x_initial,x_out];
    
elseif var(param.I_data)<1e-4
    totN=size(x_initial,1);
    sundials_u = SX.sym('sundials_u',totN);
    input_crt=SX.sym('input_crt');
    fun_ode=Module_ode_sundials(sundials_u,input_crt, param);
    dae = struct('x',sundials_u,'p',input_crt,'ode',[fun_ode],'alg',[]);
    %================aat
    additional_flag=0;
    div_num=11;
    if length(tspan)>div_num
        feasible_num=floor(length(tspan)/div_num);
        tspanuse=reshape(tspan(1:feasible_num*div_num),[],div_num);
        tspanuse=tspanuse';
        if tspanuse(end,end)~=tspan(end)
            additional_flag=1;
            additional_tspan=tspan(feasible_num*div_num+1:end);
        end
    else
        tspanuse=tspan;
    end
    
    rex=[];
    for kk=1:1:size(tspanuse,1)
        if kk==1
           tgrid=tspanuse(kk,:);
        elseif kk~=size(tspanuse,1)
           tgrid=tspanuse(kk,:)-tspanuse(kk-1,end);    
        else
            if additional_flag==1
                 tgrid=[tspanuse(kk,:),additional_tspan]-tspanuse(kk-1,end);    
            else
                tgrid=tspanuse(kk,:)-tspanuse(kk-1,end);   
            end
        end
        
        if kk==1
            opts = struct('tf',tgrid(end),'grid',tgrid,'output_t0',true,'abstol',5.0e-8,'reltol',5e-4,'linear_solver','csparse',...
            'newton_scheme','gmres','show_eval_warnings',false);  %设定仿真结束时间
            F1 = integrator('F', 'idas', dae,opts);
            u0=x_initial';
            
            sol = F1('x0',u0,'p',mean(param.I_data));
            
            
            rex=[rex,full(sol.xf)];
    
        else
            opts = struct('tf',tgrid(end),'grid',tgrid,'output_t0',true,'abstol',5.0e-8,'reltol',5e-4,'linear_solver','csparse',...
            'newton_scheme','gmres','show_eval_warnings',false);  %设定仿真结束时间
            F1 = integrator('F', 'idas', dae,opts);
            u0=rex(:,end);
            u0=u0';
            
            sol = F1('x0',u0,'p',mean(param.I_data));
            
            
            rex=[rex,full(sol.xf)];
        end
    end
    
    x_out=rex; % Only return final state
    t_out=tspan';
end
%====================
I_dummy = param.I_data(1:size(x_out,2)); %Store exact current profiles used in I_dummy

% For no additional cycles, output the relevant variables
%param.I_data = I_dummy;
%param.t_data = t_out;
dx_out = x_out(:,end);

end

% Lithium Concentration to Voltage Estimator
function V_cell = li2voltage(x_out)
global param
%% Separate electrochemical, thermal & aging state variables from x_out matrix
%Define solid concentrations
cs = x_out(1:(param.Nr-1)*param.Nc*2,:);            %All solid concentrations
cs_n = cs(1:(param.Nr-1)*param.Nc,:);               %Anode Concentrations
cs_p = cs((param.Nr-1)*param.Nc+1:end,:);           %Cathode Concentrations
if ~isreal(cs_p)                                                    
    cs_p = abs(cs_p);
end

index_cs = (param.Nr-1)*param.Nc*2;                 %Index for final solid concentration
index_ce = index_cs + param.Nc*param.ce_states;     %Index for final electrolyte concentration

%Define electrolyte concentrations
ce = x_out(index_cs+1:index_ce,:);                  %All electrolyte concentrations
ce_n = ce(1:param.Nx_n*param.Nc,:);                 %Negative Electrolyte Region
ce_s = ce(param.Nx_n*param.Nc+1:(param.Nx_n + param.Nx_s)*param.Nc,:); %Separator Region
ce_p = ce((param.Nx_n + param.Nx_s)*param.Nc+1:param.ce_states*param.Nc,:); %Positive Region

%Define cell temperatures
T_cell = x_out(index_ce+1:index_ce+param.Nc*2,:);
index_thermal = index_ce+param.Nc*2;                %Index for final temperature state

%Define aging states
index_aging = index_thermal+7*param.Nc;                                %Index for final aging state
L_sei = x_out(index_thermal+1:index_thermal+param.Nc,:);               %SEI Layer Thickness
Q = x_out(index_thermal+param.Nc+1:index_thermal+2*param.Nc,:);        %ODE-updated Capacity 
aina_n = x_out(index_thermal+2*param.Nc+1:index_thermal+3*param.Nc,:); %Inactive area evolution, anode
aina_p = x_out(index_thermal+3*param.Nc+1:index_thermal+4*param.Nc,:); %Inactive area evolution, cahode
c_sei = x_out(index_thermal+4*param.Nc+1:index_thermal+5*param.Nc,:);  %SEI concentration
c_li = x_out(index_thermal+5*param.Nc+1:index_thermal+6*param.Nc,:);   %Plated li concentration
L_film = x_out(index_thermal+6*param.Nc+1:index_thermal+7*param.Nc,:); %Film Layer Thickness

%Define SEI Solvent Concentration States
Csolv = x_out(index_aging+1:end,:); %Solvent Concentration in SEI 

%% Temperature dependent transport and kinetics for each cell
% for j = 1:length(cs)
for j = 1:size(cs,2)
    % Surface concentration -> surface stoichiometry
    for i = 1:param.Nc
        theta_surf_n(i,j) = cs_n(i*(param.Nr-1),j)/param.c_n_max;
        theta_surf_p(i,j) = cs_p(i*(param.Nr-1),j)/param.c_p_max;
        T_core(i,j) = T_cell(2*i-1,j);
        T_surf(i,j) = T_cell(2*i,j);
    end
    
    % Fracture area
    af_n(:,j) = param.a_sn*param.k_lam_n*(param.af_tinit+param.t_data(j))*param.flag_aina_n;
    af_p(:,j) = param.a_sp*param.k_lam_p*(param.af_tinit+param.t_data(j))*param.flag_aina_p;
    
    % Temperature dependence
    [Dsn(:,j), Dsp(:,j), kn(:,j), kp(:,j)] = arrhenius_temp(param,T_core(:,j));
    
    % SEI layer dependence
    [eps_el_n(:,j),R_sei(:,j)] = sei_lam_effect(L_film(:,j),param,aina_n(:,j),af_n(:,j));
    
    % LAM dependence
    eps_el_p(:,j) = lam_effect(param,aina_p(:,j),af_p(:,j));
    
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
        ce_all(:,j,i) = [ce_n(index1n:index2n,j); ce_s(index1s:index2s,j); ce_p(index1p:index2p,j)];
        
        %ce_all_avg(Row = 'Cell #', Column = 'Time point')
        ce_all_avg(i,j) = mean(ce_all(:,j,i));
        
    end
    
    % Electrolyte Conductivity
    [K_el_eff_n(:,j),K_el_eff_s(:,j),K_el_eff_p(:,j)] = conductivity_update(ce_n_avg(:,j)...
        ,ce_s_avg(:,j),ce_p_avg(:,j),eps_el_n(:,j),eps_el_p(:,j),T_core(:,j),param);
    
    % Open circuit potential and overpotential
    ocp_p(:,j) = U_p(theta_surf_p(:,j));
    ocp_n(:,j) = U_n(theta_surf_n(:,j));
    eta_p(:,j) = eta_cathode(theta_surf_p(:,j), ce_p_avg(:,j), T_core(:,j), param.I_data(j), param, kp(:,j), aina_p(:,j), af_p(:,j));
    eta_n(:,j) = eta_anode(theta_surf_n(:,j), ce_n_avg(:,j), T_core(:,j), param.I_data(j), param, kn(:,j), aina_n(:,j), af_n(:,j));  

    % Cell voltages
    [phi_e(:,j),R_el(:,j)] = electrolyte_potential(ce_all_avg(:,j), ce_all(:,j,:),K_el_eff_n(:,j),...
        K_el_eff_s(:,j), K_el_eff_p(:,j),T_core(:,j), param.I_data(j), param);
    [V_cell(:,j), R_l(:,j)] = V_calculation(ocp_p(:,j),ocp_n(:,j),eta_p(:,j),eta_n(:,j),...
        phi_e(:,j),R_sei(:,j),param.I_data(j),param.SOC_cc(j),param);
    V_oc = ocp_p - ocp_n;
    
    % Side reaction current density
    phi_sn(:,j) = ocp_n(:,j) + eta_n(:,j) + R_sei(:,j)*param.I_data(j);
    i_s(:,j) = side_current_solvent(R_sei(:,j), T_core(:,j), ...
        phi_sn(:,j), param, param.I_data(j), theta_surf_n(:,j),Csolv(:,j));
    i_lpl(:,j) = side_current_plating(R_sei(:,j),T_core(:,j),phi_sn(:,j),phi_e(:,j),param,param.I_data(j));
end

% figure; hold on
% plot(eta_p); plot(-eta_n)

%% Calculate SOC
[soc_bulk_n,soc_bulk_p] = soc_calculation(cs_n,cs_p,param);
end