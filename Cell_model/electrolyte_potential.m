function [phi_e, R_el] = electrolyte_potential(ce_all_avg, ce_all,K_el_eff_n,K_el_eff_s,K_el_eff_p,T_core, input_crt, param)
% Electrolyte potential calculation
%Convert from mol/m3 to mol/L
M_all = ce_all_avg/1000;

%Def'n: nu(C,T) = (1 - t0)*(1 + dlnf / dlnCe)
nu_T = 0.602 - 0.24*M_all.^0.5 + 0.982*(1-0.0052*(T_core-293)).*M_all.^1.5;

%NOTE: changed structure of ce_all to incorporate time-steps (in main
%script)
R_el = (param.Ln.*K_el_eff_n.^(-1) + 2*param.Ls.*K_el_eff_s.^(-1) + ...
    param.Lp.*K_el_eff_p.^(-1)).*(2*param.A).^(-1);

% phi_e = -R_el*input_crt + 2*param.Rg*T_core.*...
%     (1-param.t0).*nu_T.*((K_el_eff_n+K_el_eff_s+K_el_eff_p)/3).*...
%     (log(ce_all(end,:)')-log(ce_all(1,:)'))./param.F;

% Removed term "(1-param.t0)" because nu(T) already includes (1-param.t0)
phi_e = -R_el.*input_crt + 2*param.Rg*T_core.*...
    nu_T.*((K_el_eff_n+K_el_eff_s+K_el_eff_p)/3).*...
    (log(ce_all(end,:)')-log(ce_all(1,:)'))./param.F;

end