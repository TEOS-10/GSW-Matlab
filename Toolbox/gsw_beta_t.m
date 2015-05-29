function result = gsw_beta_t(SA,t,p)

%% 
% result = gsw_beta_t(SA,t,p)
%
% haline contraction coefficient of seawater 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
%
% result              : haline contraction coefficient     [kg/g]
%                       at constant (in situ) temperature

%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_beta_t do not agree    ****')
end

n0 = 0; n1 = 1; 

result = -gsw_g(n1,n0,n1,SA,t,p)./gsw_g(n0,n0,n1,SA,t,p);

end