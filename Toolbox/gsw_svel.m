function result = gsw_svel(SA,t,p)

%% 
% result = gsw_svel(SA,t,p)
%
% sound speed of seawater 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
%
% result              : sound speed                        [m/sec]

%%

if gsw_check_arrays(SA,t,p)
    error('****    input array dimensions in gsw_svel do not agree    ****')
end

n0 = 0; n1 = 1; n2 = 2;

g_tt = gsw_g(n0,n2,n0,SA,t,p); g_tp = gsw_g(n0,n1,n1,SA,t,p);

result = gsw_g(n0,n0,n1,SA,t,p) .* ...
             sqrt(g_tt./(g_tp.*g_tp - g_tt.*gsw_g(n0,n0,n2,SA,t,p)));
         
end