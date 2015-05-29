function result = gsw_enthalpy(SA,t,p)

%% 
% result = gsw_enthalpy(SA,t,p)
%
% specific enthalpy of seawater 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
%
% result              : specific enthalpy                  [J/kg]


%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_enthalpy do not agree    ****')
end

n0 = 0; n1 = 1;

result = gsw_g(n0,n0,n0,SA,t,p) - (t+273.15).*gsw_g(n0,n1,n0,SA,t,p);

end