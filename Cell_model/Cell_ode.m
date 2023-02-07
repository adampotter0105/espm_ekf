function [dx_dt] = Cell_ode(t_in, x_in, param)
% Solving solid phase PDE using ODE solver
% Interpolate current
input_crt = interp1(param.t_data, param.I_data, t_in);   


%% Separate electrochemical, thermal & aging state variables from x_in vector
%Define solid concentrations
cs = x_in(1:(param.Nr-1)*2);  %Solid Concentration states
cs_n = cs(1:(param.Nr-1));             %Anode Solid Conc.
cs_p = cs((param.Nr-1)+1:end);         %Cathode Solid Conc. 
index_cs = (param.Nr-1)*2;             %Index for final solid concentration

%Define electrolyte concentrations
index_ce = index_cs + param.Nc*param.ce_states; %Index for final electrolyte concentration
ce = x_in(index_cs+1:index_ce);                 %All electrolyte concentrations
ce_n = ce(1:param.Nx_n*param.Nc);               %Negative Electrolyte Region
ce_s = ce(param.Nx_n+1:(param.Nx_n + param.Nx_s)); %Separator Region
ce_p = ce((param.Nx_n + param.Nx_s)+1:param.ce_states); %Positive Region

%% Temperature dependent transport and kinetics for each cellparam.kappa_sei
%Allocate memory
theta_surf_n(param.Nc,1) = cs_n(param.Nr-1)/param.c_n_max;
theta_surf_p(param.Nc,1) = cs_p(param.Nr-1)/param.c_p_max;

%% Update parameters based on aging & temperature
T_core = param.T_amb; 
T_surf = param.T_amb; 

eps_el_n = param.eps_el_n_ref;
% Temperature dependence
[Dsn, Dsp, kn, kp] = arrhenius_temp(param, T_core);

%% Calculate average electrolyte concentrations
%Allocate Memory
ce_n_avg(param.Nc,1) = 0;
ce_s_avg(param.Nc,1) = 0;
ce_p_avg(param.Nc,1) = 0;
ce_all((param.Nx_n+param.Nx_s+param.Nx_p),param.Nc) = 0;
ce_all_avg(param.Nc,1) = 0;

for i = 1:param.Nc

    index1n = (i-1)*param.Nx_n+1; %1st negative electrolyte grid point for cell 'i'
    index2n = i*param.Nx_n;       %Last negative electrolyte grid point for cell 'i'
    index1s = (i-1)*param.Nx_s+1; %1st separator grid point for cell 'i'
    index2s = i*param.Nx_s;       %Last separator grid point for cell 'i'
    index1p = (i-1)*param.Nx_p+1; %1st positive electrolyte grid point for cell 'i'
    index2p = i*param.Nx_p;       %Last positive electrolyte grid point for cell 'i'
    
   %Average Electrolyte Concentrations
    ce_n_avg(i,1) = mean(ce_n(index1n:index2n));
    ce_s_avg(i,1) = mean(ce_s(index1s:index2s));
    ce_p_avg(i,1) = mean(ce_p(index1p:index2p));
    ce_all(:,i) = [ce_n(index1n:index2n); ce_s(index1s:index2s); ce_p(index1p:index2p)];
    ce_all_avg(i,1) = mean(ce_all(:,i));
    
end
%% Update Electrolyte Concentration- / Temperature- / Aging - dependent parameters
% Electrolyte Diffusivity

[De_n,De_s,De_p] = diffusivity_update(ce_n,ce_s,ce_p,T_core,param);

    De_n_minus0_5 = zeros(param.Nc*param.Nx_n,1);
    De_n_plus0_5 = zeros(param.Nc*param.Nx_n,1);
    
    De_s_minus0_5 = zeros(param.Nc*param.Nx_s,1);
    De_s_plus0_5 = zeros(param.Nc*param.Nx_s,1);
    
    De_p_minus0_5 = zeros(param.Nc*param.Nx_p,1);
    De_p_plus0_5 = zeros(param.Nc*param.Nx_p,1);

for i = 1:param.Nc
    
    index1n = (i-1)*param.Nx_n+1;
    index2n = i*param.Nx_n;       
    index1s = (i-1)*param.Nx_s+1; 
    index2s = i*param.Nx_s;       
    index1p = (i-1)*param.Nx_p+1; 
    index2p = i*param.Nx_p;       
    
    %Create vectors of diffusivity at grid points i - 1 or i + 1
    De_n_minus1 = [0; De_n(index1n:index2n-1,1)];
    De_n_plus1 = [De_n(index1n+1:index2n,1); 0];
    De_n_dummy = De_n(index1n:index2n,1);
    
    %Calculate diffusivities at grid points i + 1/2 and i - 1/2
    De_n_minus0_5(index1n:index2n) = 2*De_n_dummy.*De_n_minus1./(De_n_minus1 + De_n_dummy);
    De_n_plus0_5(index1n:index2n) = 2*De_n_dummy.*De_n_plus1./(De_n_plus1 + De_n_dummy);
    
    
    De_s_minus1 = [0; De_s(index1s:index2s-1,1)];
    De_s_plus1 = [De_s(index1s+1:index2s,1); 0];
    De_s_dummy = De_s(index1s:index2s,1);
    
    De_s_minus0_5(index1s:index2s) = 2*De_s_dummy.*De_s_minus1./(De_s_minus1 + De_s_dummy);
    De_s_plus0_5(index1s:index2s) = 2*De_s_dummy.*De_s_plus1./(De_s_plus1 + De_s_dummy);
    
    De_p_minus1 = [0; De_p(index1p:index2p-1,1)];
    De_p_plus1 = [De_p(index1p+1:index2p,1); 0];
    De_p_dummy = De_p(index1p:index2p,1);
    
    De_p_minus0_5(index1p:index2p) = 2*De_p_dummy.*De_p_minus1./(De_p_minus1 + De_p_dummy);
    De_p_plus0_5(index1p:index2p) = 2*De_p_dummy.*De_p_plus1./(De_p_plus1 + De_p_dummy);
    
    %Define Electrolyte Interface Effective Diffusivities
    %Note: these formulas assume same geometry for each cell
    beta_n_s = param.delta_n./(param.delta_n + param.delta_s);
    De_eff_n_s(i) = (De_n(index2n)*eps_el_n(i)^param.brugg_n)*(De_s(index1s)*param.eps_el_s(i)^param.brugg_s)/(beta_n_s(i)*...
        (De_s(index1s)*param.eps_el_s(i)^param.brugg_s) + (1-beta_n_s(i))*(De_n(index2n)*eps_el_n(i)^param.brugg_n));
    param.delta_n_s = 1/2*(param.delta_n + param.delta_s);
    
    beta_s_p = param.delta_s./(param.delta_s + param.delta_p);
    De_eff_s_p(i) = (De_s(index2s)*param.eps_el_s(i)^param.brugg_s)*(De_p(index1p)*param.eps_el_p(i)^param.brugg_p)/(beta_s_p(i)*...
        (De_p(index1p)*param.eps_el_p(i)^param.brugg_p) + (1-beta_s_p(i))*(De_s(index2s)*param.eps_el_s(i)^param.brugg_s));
    param.delta_s_p = 1/2*(param.delta_s + param.delta_p);
    
end

% Electrolyte Conductivity
[K_el_eff_n,K_el_eff_s,K_el_eff_p] = conductivity_update(ce_n_avg,ce_s_avg,ce_p_avg,eps_el_n,T_core,param);

% Open circuit potential and overpotential
ocp_p = U_p(theta_surf_p);
ocp_n = U_n(theta_surf_n);

I_cell = input_crt; 

eta_p = eta_cathode(theta_surf_p, ce_p_avg, T_core, I_cell, param, kp);
eta_n = eta_anode(theta_surf_n, ce_n_avg, T_core, I_cell, param, kn);

% Cell voltages
[phi_e, R_el] = electrolyte_potential(ce_all_avg, ce_all,K_el_eff_n,K_el_eff_s,K_el_eff_p,T_core, I_cell, param);
V_cell = ocp_p - ocp_n + eta_p - eta_n + phi_e - param.R_l.*I_cell;


%% Solve solid phase ODEs

% Coefficients of discretized ODEs
alpha_n = Dsn/(param.delta_xn^2);
alpha_p = Dsp/(param.delta_xp^2);

W_mat_n = zeros(param.Nc*(param.Nr-1)); %Coeff. matrix for all anode solid conc. states for Nc cells
W_mat_p = zeros(param.Nc*(param.Nr-1)); %Coeff. matrix for all cathode solid conc. states for Nc cells

for i = 1:param.Nc; 
    %Define W_mat for all solid conc. states for each cell
    %Repmat creates W_mat by tiling Alpha(i) in (Nr - 1)x(Nr -1) matrix
    W_mat_n((i-1)*(param.Nr - 1)+1:i*(param.Nr - 1),(i-1)*(param.Nr - 1)+1:i*(param.Nr - 1)) = ...
        repmat(alpha_n(i,:),(param.Nr - 1),(param.Nr - 1));
    W_mat_p((i-1)*(param.Nr - 1)+1:i*(param.Nr - 1),(i-1)*(param.Nr - 1)+1:i*(param.Nr - 1)) = ...
        repmat(alpha_p(i,:),(param.Nr - 1),(param.Nr - 1));    
end

A_csn = param.A_sd.*W_mat_n; %Full Anode Solid Conc. State A-matrix
A_csp = param.A_sd.*W_mat_p; %%Full Cathode Solid Conc. State A-matrix



%Use modified intercalation current for anode
%NOTE: Sign convention is such that '+ i_s' is correct sign in
%inter_crt expression
dcsn_dt = A_csn*cs_n + param.B_sd_n.*(I_cell);
dcsp_dt = A_csp*cs_p + param.B_sd_p.*I_cell;

dcs_dt = [dcsn_dt; dcsp_dt];

%% Solve Electrolyte Phase ODE
% Coefficients of discretized electrolyte ODE's
alpha_el_n = eps_el_n.^(param.brugg_n-1)./(param.delta_n.^2);

alpha_el_s = param.eps_el_s.^(param.brugg_s-1)./(param.delta_s.^2);

alpha_el_p = param.eps_el_p.^(param.brugg_p-1)./(param.delta_p.^2);

beta_el_n = (1-param.t0)./(eps_el_n.*param.F.*param.A.*param.Ln);
beta_el_s = 0; 
beta_el_p = (1-param.t0)./(param.eps_el_p.*param.F.*param.A.*param.Lp);
%NOTE: +/- for anode / cathode is incorporated later

%DEFINE W_mat matrices
W_e_n = zeros(param.Nc*param.Nx_n);
W_e_s = zeros(param.Nc*param.Nx_s);
W_e_p = zeros(param.Nc*param.Nx_p);

for i = 1:param.Nc; 
    %Define W_mat for all electrolyte conc. states for each cell
    %Repmat creates W_mat by tiling Alpha(i) in (N_x)x(N_x) matrix
    W_e_n((i-1)*param.Nx_n+1:i*param.Nx_n,(i-1)*param.Nx_n+1:i*param.Nx_n) = ...
        repmat(alpha_el_n(i),param.Nx_n,param.Nx_n);

    W_e_s((i-1)*param.Nx_s+1:i*param.Nx_s,(i-1)*param.Nx_s+1:i*param.Nx_s) = ...
        repmat(alpha_el_s(i),param.Nx_s,param.Nx_s);
    
    W_e_p((i-1)*param.Nx_p+1:i*param.Nx_p,(i-1)*param.Nx_p+1:i*param.Nx_p) = ...
        repmat(alpha_el_p(i),param.Nx_p,param.Nx_p);
    
end

 Y_e_n = [];
 Y_e_s = [];
 Y_e_p = [];
 
    for i = 1:param.Nc
    n_dummy = beta_el_n(i)*ones(param.Nx_n,1);
    Y_e_n = [Y_e_n; n_dummy];
    
    s_dummy = beta_el_s*ones(param.Nx_s,1);
    Y_e_s = [Y_e_s; s_dummy];
    
    p_dummy = beta_el_p(i)*ones(param.Nx_p,1);
    Y_e_p = [Y_e_p; p_dummy];
    end

    A_tot_n = De_n_minus0_5.*param.A1_en_d + De_n_plus0_5.*param.A2_en_d;
    A_tot_s = De_s_minus0_5.*param.A1_es_d + De_s_plus0_5.*param.A2_es_d;
    A_tot_p = De_p_minus0_5.*param.A1_ep_d + De_p_plus0_5.*param.A2_ep_d;
    
A_e_n = A_tot_n.*W_e_n;

B_e_n = Y_e_n.*param.B_en_d;

A_e_s = A_tot_s.*W_e_s;
B_e_s = Y_e_s.*param.B_es_d;

A_e_p = A_tot_p.*W_e_p;

B_e_p = Y_e_p.*param.B_ep_d;
    
%% Define Electrolyte Interface ODE's
%NOTE: interfaces between electrolyte regions are handled separately to simplify indexing
%Can try streamlining to handle interfaces with same algorithm as bulk in future versions
for i = 1:param.Nc
    index2n = i*param.Nx_n;
    index1s = (i-1)*param.Nx_s+1;
    index2s = i*param.Nx_s;
    index1p = (i-1)*param.Nx_p+1;
    
dCen_interf_ns(i) = 1/eps_el_n(i)*De_eff_n_s(i)/(param.delta_n(i)*param.delta_n_s(i))*(ce_s(index1s) - ce_n(index2n))...
    -De_n_minus0_5(index2n)*eps_el_n(i)^(param.brugg_n-1)/(param.delta_n(i)^2)*(ce_n(index2n) - ce_n(index2n-1)) + beta_el_n(i)*I_cell(i);

dCes_interf_ns(i) = De_s_plus0_5(index1s)*param.eps_el_s(i)^(param.brugg_s-1)/(param.delta_s(i)^2)*(ce_s(index1s+1) - ce_s(index1s))...
    -1/param.eps_el_s(i)*De_eff_n_s(i)/(param.delta_s(i)*param.delta_n_s(i))*(ce_s(index1s) - ce_n(index2n));


dCes_interf_sp(i) = 1/param.eps_el_s(i)*De_eff_s_p(i)/(param.delta_s(i)*param.delta_s_p(i))*(ce_p(index1p) - ce_s(index2s))...
    -De_s_minus0_5(index2s)*param.eps_el_s(i)^(param.brugg_s-1)/(param.delta_s(i)^2)*(ce_s(index2s) - ce_s(index2s-1));

dCep_interf_sp(i) = De_p_plus0_5(index1p)*param.eps_el_p(i)^(param.brugg_s-1)/(param.delta_p(i)^2)*(ce_p(index1p+1) - ce_p(index1p))...
    -1/param.eps_el_p(i)*De_eff_s_p(i)/(param.delta_p(i)*param.delta_s_p(i))*(ce_p(index1p) - ce_s(index2s)) - beta_el_p(i)*I_cell(i);
end

%%   
% Define Electrolyte phase ODE's

% Define dummy vectors for I_cell of appropriate size
I_cell_dummy_en = ones(param.Nc*(param.Nx_n),1);
I_cell_dummy_ep = ones(param.Nc*(param.Nx_p),1);

for i = 1:param.Nc

I_cell_dummy_en((i-1)*(param.Nx_n)+1:i*(param.Nx_n),1) = I_cell_dummy_en((i-1)*...
    (param.Nx_n)+1:i*(param.Nx_n),1)*I_cell(i);
I_cell_dummy_ep((i-1)*(param.Nx_p)+1:i*(param.Nx_p),1) = I_cell_dummy_ep((i-1)*...
    (param.Nx_p)+1:i*(param.Nx_p),1)*I_cell(i);

end

%NOTE: 
%Input current term is added in negative electrolyte region
%Input current term is subtracted in positive electrolyte region
dce_n = A_e_n*ce_n + B_e_n.*I_cell_dummy_en;

% dce_s = A_e_s*ce_s + B_e_s*I_cell;
dce_s = A_e_s*ce_s; %NO intercalation/deintercalation in separator - don't complicate it w/ I_cell term


dce_p = A_e_p*ce_p - B_e_p.*I_cell_dummy_ep;

for i = 1:param.Nc
    
    index2n = i*param.Nx_n;
    index1s = (i-1)*param.Nx_s+1;
    index2s = i*param.Nx_s;
    index1p = (i-1)*param.Nx_p+1;
    
    dce_n(index2n) = dCen_interf_ns(i);
    dce_s(index1s) = dCes_interf_ns(i);
    dce_s(index2s) = dCes_interf_sp(i);
    dce_p(index1p) = dCep_interf_sp(i);
end

dce_dt = [dce_n;dce_s;dce_p];



%% Final formulation
dx_dt = [dcs_dt; dce_dt];
end

