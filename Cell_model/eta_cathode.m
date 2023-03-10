function eta_p = eta_cathode(x, ce_p_avg, T_core, input_crt, param, kp)

% Calculates overpotential for the specified surface concentration

i0_p= kp*param.F.*(ce_p_avg.^0.5).*...
    ((param.c_p_max.*x).^0.5).*...
    ((param.c_p_max-param.c_p_max.*x).^0.5);   

eta_p = asinh(input_crt./(2*param.A.*param.a_sp.*...
    param.Lp.*i0_p)).*(param.Rg.*T_core)./...
    (param.F*param.alpha_cell);


end