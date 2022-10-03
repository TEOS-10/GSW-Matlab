function beta_const_pt_t_exact = gsw_beta_const_pt_t_exact(SA,t,p)

% gsw_beta_const_pt_t_exact                  saline contraction coefficient
%                                         at constant potential temperature
%==========================================================================
%
% USAGE:  
%  beta_const_pt_t_exact = gsw_beta_const_pt_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the saline (i.e. haline) contraction coefficient of seawater   
%  at constant potential temperature with a reference pressure of 0 dbar.
%   
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  beta_const_pt_t_exact  =  saline contraction coefficient        [ kg/g ]
%                            at constant potential temperature
%                            and with a reference pressure of 0 dbar.
%  
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.19.2) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_beta_const_pt_t_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_beta_const_pt_t_exact: SA and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_beta_const_pt_t_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

db2Pa = 1e-4;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).

pt0 = gsw_pt0_from_t(SA,t,p);

x2 = sfac.*SA;
x = sqrt(x2);
y = 0.025*t;
y_pt = 0.025*pt0;
z = db2Pa*p; %Note.The input pressure (p) is sea pressure in units of dbar.

g_SA_T_mod = 1187.3715515697959 + z.*(1458.233059470092 + ...
        z.*(-687.913805923122 + z.*(249.375342232496 + z.*(-63.313928772146 + 14.09317606630898.*z)))) + ...
        x.*(-1480.222530425046 + x.*(2175.341332000392 + x.*(-980.14153344888 + 220.542973797483.*x) + ...
        y.*(-548.4580073635929 + y.*(592.4012338275047 + y.*(-274.2361238716608 + 49.9394019139016.*y))) - ...
        90.6734234051316.*z) + z.*(-525.876123559641 + (249.57717834054571 - 88.449193048287.*z).*z) + ...
        y.*(-258.3988055868252 + z.*(2298.348396014856 + z.*(-325.1503575102672 + 153.8390924339484.*z)) + ...
        y.*(-90.2046337756875 - 4142.8793862113125.*z + y.*(10.50720794170734 + 2814.78225133626.*z)))) + ...
        y.*(3520.125411988816 + y.*(-1351.605895580406 + ...
        y.*(731.4083582010072 + y.*(-216.60324087531103 + 25.56203650166196.*y) + ...
        z.*(-2381.829935897496 + (597.809129110048 - 291.8983352012704.*z).*z)) + ...
        z.*(4165.4688847996085 + z.*(-1229.337851789418 + (681.370187043564 - 66.7696405958478.*z).*z))) + ...
        z.*(-3443.057215135908 + z.*(1349.638121077468 + ...
        z.*(-713.258224830552 + (176.8161433232 - 31.68006188846728.*z).*z))));  
g_SA_T_mod = 0.5*sfac*0.025*g_SA_T_mod;

g_SA_pt_mod = 1187.3715515697959 + ...
        x.*(-1480.222530425046 + x.*(2175.341332000392 + x.*(-980.14153344888 + 220.542973797483.*x) + ...
        y_pt.*(-548.4580073635929 + y_pt.*(592.4012338275047 + y_pt.*(-274.2361238716608 + 49.9394019139016.*y_pt)))) + ...
        y_pt.*(-258.3988055868252 + ...
        y_pt.*(-90.2046337756875 + y_pt.*(10.50720794170734)))) + ...
        y_pt.*(3520.125411988816 + y_pt.*(-1351.605895580406 + ...
        y_pt.*(731.4083582010072 + y_pt.*(-216.60324087531103 + 25.56203650166196.*y_pt))));
g_SA_pt_mod = 0.5*sfac*0.025*g_SA_pt_mod;

gp = gsw_gibbs(0,0,1,SA,t,p);

factor = (g_SA_T_mod - g_SA_pt_mod)./(gp.*gsw_gibbs(0,2,0,SA,t,p));

beta_const_pt_t_exact = gsw_gibbs(0,1,1,SA,t,p).*factor ...
                         - gsw_gibbs(1,0,1,SA,t,p)./gp;
                     
if transposed
    beta_const_pt_t_exact = beta_const_pt_t_exact.';
end

end
