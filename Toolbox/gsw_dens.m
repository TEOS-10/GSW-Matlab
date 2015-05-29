function result = gsw_dens(SA,t,p)

%% 
% result = gsw_dens(SA,t,p)
%
% density of seawater 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
%
% result              : density                            [kg/m^3]

%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_dens do not agree    ****')
end

n0 = 0; n1 = 1;

result = 1d0./gsw_g(n0,n0,n1,SA,t,p);

return