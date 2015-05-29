function result = gsw_alpha_t(SA,t,p)

%%
% result = gsw_alpha_t(SA,t,p)
%
% thermal expansion coefficient of seawater 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
%
% result              : thermal expansion coefficient      [1/K]
%                       wrt (in situ) temperature

%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_alpha_t do not agree    ****')
end

n0 = 0; n1 = 1;

result = gsw_g(n0,n1,n1,SA,t,p)./gsw_g(n0,n0,n1,SA,t,p);

end