function result = gsw_cp(SA,t,p)

%% 
% result = gsw_cp(SA,t,p)
%
% isobaric heat capacity of seawater
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
%
% result              : heat capacity                      [J/(kg*K)]

%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_cp do not agree    ****')
end

n0 = 0; n2 = 2;

result = -(t+273.15).*gsw_g(n0,n2,n0,SA,t,p);

end