function result = gsw_kappa_t(SA,t,p)

%% 
% result = gsw_kappa_t(SA,t,p)
%
% isothermal compressibility of seawater 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
%
% result              : isothermal compressibility         [1/Pa]


%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_kappa_t do not agree    ****')
end

n0 = 0; n1 = 1; n2 = 2;

result = -1.d4 * gsw_g(n0,n0,n2,SA,t,p)./gsw_g(n0,n0,n1,SA,t,p);

end