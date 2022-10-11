function [rho_SA, rho_CT, rho_P] = gsw_rho_first_derivatives_CT_exact(SA,CT,p)

% gsw_rho_first_derivatives_CT_exact       SA, CT and p partial derivatives
%                                                                of density
%==========================================================================
% 
% USAGE:  
% [rho_SA, rho_CT, rho_P] = gsw_rho_first_derivatives_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the three (3) partial derivatives of in-situ density with 
%  respect to Absolute Salinity, Conservative Temperature and pressure.  
%  Note that the pressure derivative is done with respect to pressure in 
%  Pa, not dbar.  
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_rho_first_derivatives(SA,CT,p), which uses the computationally
%  efficient 75-term expression for specific volume in terms of SA, CT and 
%  p (Roquet et al., 2015)
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_SA  =  partial derivative of density           [ (kg/m^3)(g/kg)^-1 ]
%                 with respect to Absolute Salinity
%  rho_CT  =  partial derivative of density                  [ kg/(m^3 K) ]
%                 with respect to Conservative Temperature
%  rho_P   =  partial derivative of density                 [ kg/(m^3 Pa) ]
%                 with respect to pressure in Pa
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual. 
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_rho_first_derivatives_CT_exact:  Requires three inputs')
end %if
if ~(nargout == 3 | nargout == 4)
   error('gsw_rho_first_derivatives_CT_exact:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_rho_first_derivatives_CT_exact: SA and CT must have same dimensions')
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
    error('gsw_rho_first_derivatives_CT_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

t = gsw_t_from_CT(SA,CT,p);
T0 = gsw_T0;

db2Pa = 1e-4;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).

pt0 = gsw_pt0_from_t(SA,t,p);

x2 = sfac.*SA;
x = sqrt(x2);
y = 0.025*t;
y_pt = 0.025*pt0;
z = db2Pa*p; %Note.The input pressure (p) is sea pressure in units of dbar.

g_SA_T_mod_part = 1187.3715515697959 + z.*(1458.233059470092 + ...
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
g_SA_T_mod = 0.5*sfac*0.025*g_SA_T_mod_part;
   
g_SA_mod_part = 8645.36753595126 + ...
        x.*(-7296.43987145382 + x.*(8103.20462414788 + ...
        y_pt.*(2175.341332000392 + y_pt.*(-274.2290036817964 + ...
        y_pt.*(197.4670779425016 + y_pt.*(-68.5590309679152 + 9.98788038278032.*y_pt)))) + ...
        x.*(-5458.34205214835 - 980.14153344888.*y_pt + ...
        x.*(2247.60742726704 - 340.1237483177863.*x + 220.542973797483.*y_pt))) + ...
        y_pt.*(-1480.222530425046 + ...
        y_pt.*(-129.1994027934126 + ...
        y_pt.*(-30.0682112585625 + y_pt.*(2.626801985426835 ))))) + ...
        y_pt.*(1187.3715515697959 + ...
        y_pt.*(1760.062705994408 + y_pt.*(-450.535298526802 + ...
        y_pt.*(182.8520895502518 + y_pt.*(-43.3206481750622 + 4.26033941694366.*y_pt)))));
g_SA_mod = 0.5*sfac*g_SA_mod_part;   

g_p = gsw_gibbs(0,0,1,SA,t,p);

g_psq_g_tt = g_p.*g_p.*gsw_gibbs(0,2,0,SA,t,p);

g_tp = gsw_gibbs(0,1,1,SA,t,p);

factora = g_SA_T_mod - g_SA_mod./(T0+pt0);

factor = factora./g_psq_g_tt;

rho_SA = g_tp.*factor - gsw_gibbs(1,0,1,SA,t,p)./(g_p.*g_p);

rho_CT = g_tp.*gsw_cp0./((T0 + pt0).*g_psq_g_tt);

rho_P = (g_tp.^2 - gsw_gibbs(0,2,0,SA,t,p).*gsw_gibbs(0,0,2,SA,t,p))./(g_psq_g_tt);

if transposed
    rho_SA = rho_SA.'; 
    rho_CT = rho_CT.';
    rho_P = rho_P.';
end

end
