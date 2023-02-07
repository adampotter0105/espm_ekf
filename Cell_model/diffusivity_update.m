function [De_n,De_s,De_p] = diffusivity_update(ce_n,ce_s,ce_p,T_core,param)

%Update Diffusivity in each region based on changing concentration
%Define molarities of each electrolyte region [mol/L] (convert from mol/m3)

M_n = ce_n/1000; 
M_s = ce_s/1000;
M_p = ce_p/1000;

%Create vectors of core temperatures to match dimensions of Concentration vectors
T_c_n = repelem(T_core,param.Nx_n)';
T_c_s = repelem(T_core,param.Nx_s)';
T_c_p = repelem(T_core,param.Nx_p)';


%NOTE: De_scaling factor included for sensitivity analysis
De_n = repelem(param.De_scaling,param.Nx_n).*10.^(-(4.43 + 54./(T_c_n - (229+5*M_n)) + 0.22*M_n))*(1/1e4); %[m2/s]

De_s = repelem(param.De_scaling,param.Nx_s).*10.^(-(4.43 + 54./(T_c_s - (229+5*M_s)) + 0.22*M_s))*(1/1e4); %[m2/s]

De_p = repelem(param.De_scaling,param.Nx_p).*10.^(-(4.43 + 54./(T_c_p - (229+5*M_p)) + 0.22*M_p))*(1/1e4); %[m2/s]


end