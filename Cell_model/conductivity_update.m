function [K_el_eff_n,K_el_eff_s,K_el_eff_p] = conductivity_update(ce_n_avg,ce_s_avg,ce_p_avg,eps_el_n,T_core,param)

%Define molarities of each electrolyte region [mol/L] (convert from mol/m3)
% T_c_celsius = T_core - 273; %For some reason, Valoen chose to use fit conductivity function in Celsisus 
T_c_celsius = T_core; %CORRECTED: It needs to be in [K], not Celsius 


M_n = ce_n_avg/1000; 
M_s = ce_s_avg/1000;
M_p = ce_p_avg/1000;

K_el_temp_n = (M_n.*((-10.5 + 0.074*T_c_celsius - 6.96e-5*T_c_celsius.^2) + ...
    M_n.*(0.668 - 0.0178*T_c_celsius + 2.8e-5*T_c_celsius.^2) + ...
    M_n.^2.*(0.494 - 8.86e-4*T_c_celsius)).^2).*1e-3/1e-2;
K_el_temp_p = (M_p.*((-10.5 + 0.074*T_c_celsius - 6.96e-5*T_c_celsius.^2) + ...
    M_p.*(0.668 - 0.0178*T_c_celsius + 2.8e-5*T_c_celsius.^2)+ ...
    M_p.^2.*(0.494 - 8.86e-4*T_c_celsius)).^2).*1e-3/1e-2;
K_el_temp_s = (M_s.*((-10.5 + 0.074*T_c_celsius - 6.96e-5*T_c_celsius.^2) + ...
    M_s.*(0.668 - 0.0178*T_c_celsius + 2.8e-5*T_c_celsius.^2)+ ...
    M_s.^2.*(0.494 - 8.86e-4*T_c_celsius)).^2).*1e-3/1e-2;


%NOTE: K_el_scaling factor included for sensitivity analyses
K_el_eff_n = param.K_el_scaling.*K_el_temp_n.*eps_el_n.^(param.brugg_n);
K_el_eff_s = param.K_el_scaling.*K_el_temp_s.*param.eps_el_s.^(param.brugg_s);
K_el_eff_p = param.K_el_scaling.*K_el_temp_p.*param.eps_el_p.^(param.brugg_p);

end
