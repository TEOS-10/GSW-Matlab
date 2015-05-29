function result = gsw_pden(SA,t,p,pr)

%% 
% result = gsw_pden(SA,t,p,pr)
%
% potential density of seawater 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
% pr                  : reference (gauge) pressure         [dbar]
%
% result              : potential density                  [kg/m^3]
%

%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_pden do not agree    ****')
end

theta = gsw_ptmp(SA,t,p,pr);

result = gsw_dens(SA,theta,pr);

end