function [Dsn, Dsp, kn, kp] = arrhenius_temp(param, T_core)

Dsn = param.Dsn_ref.*exp((1/param.Tref - 1./T_core).*param.Ea_Dsn/param.Rg);
Dsp = param.Dsp_ref.*exp((1/param.Tref - 1./T_core).*param.Ea_Dsp/param.Rg);
kn = param.kn_ref.*exp((1/param.Tref - 1./T_core).*param.Ea_kn/param.Rg);
kp = param.kp_ref.*exp((1/param.Tref - 1./T_core).*param.Ea_kp/param.Rg);
%Dsolv = param.Dsolv_ref.*exp((1/param.Tref - 1./T_core).*param.Ea_Dsolv/param.Rg);

end