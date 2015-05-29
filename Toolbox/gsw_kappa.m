function result = gsw_kappa(SA,t,p)

%% 
% result = gsw_kappa(SA,t,p)
%
% isentropic compressibility of seawater 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
%
% result              : isentropic compressibility         [1/dbar]

%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_kappa do not agree    ****')
end

n0 = 0; n1 = 1; n2 = 2;

g_tt = gsw_g(n0,n2,n0,SA,t,p); g_tp = gsw_g(n0,n1,n1,SA,t,p);

result = 1.d4 * (g_tp.*g_tp - g_tt.*gsw_g(n0,n0,n2,SA,t,p)) ./ (gsw_g(n0,n0,n1,SA,t,p).*g_tt);

end