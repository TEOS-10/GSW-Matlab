function [h_SA, h_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT,p)

% gsw_enthalpy_first_derivatives_CT_exact     first derivatives of enthalpy
%==========================================================================
%
% USAGE:
%  [h_SA, h_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following two derivatives of specific enthalpy, h,
%   (1) h_SA, the derivative with respect to Absolute Salinity at 
%       constant CT and p, and
%   (2) h_CT, derivative with respect to CT at constant SA and p. 
%  Note that h_P is specific volume, v, it can be calulated by calling
%  gsw_specvol_CT_exact(SA,CT,p).
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_enthalpy_first_derivatives(SA,CT,p) which uses the computationally
%  efficient 75-term expression for specific volume in terms of SA, CT and
%  p (Roquet et al., 2015).   
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
%  h_SA  =  The first derivative of specific enthalpy with respect to 
%           Absolute Salinity at constant CT and p.     
%                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
%  h_CT  =  The first derivative of specific enthalpy with respect to 
%           CT at constant SA and p.                           [ J/(kg K) ]
%
% AUTHOR: 
%  Trevor McDougall.                                   [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%    See Eqns. (A.11.18), (A.11.15) and (A.11.12) of this TEOS-10 Manual.   
%
%  McDougall, T.J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%    See Eqns. (18) and (22)
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_enthalpy_first_derivatives_CT_exact:  Requires three inputs')
end %if

if ~(nargout == 2)
   error('gsw_enthalpy_first_derivatives_CT_exact:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
    error('gsw_enthalpy_first_derivatives_CT_exact: SA and CT do not have the same dimensions')
end %if

if (mp == 1) & (np == 1)              % p is a scalar - fill to size of SA.
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_enthalpy_first_derivatives_CT_exact: The dimensions of p do not agree')
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
pt0 = gsw_pt_from_CT(SA,CT);  

temp_ratio = (gsw_T0 + t)./(gsw_T0 + pt0);

h_CT = gsw_cp0.*temp_ratio;         % from Eqn. (A.11.15) of IOC et al. (2010).

db2Pa = 1e-4;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).

x = sqrt(sfac.*SA);
y = 0.025*t;
z = db2Pa*p; %Note.The input pressure (p) is sea pressure in units of dbar.

g_SA_mod_t = 8645.36753595126 + z.*(-6620.98308089678 + ...
        z.*(769.588305957198 + z.*(-193.0648640214916 + (31.6816345533648 - 5.24960313181984.*z).*z))) + ...
        x.*(-7296.43987145382 + x.*(8103.20462414788 + ...
        y.*(2175.341332000392 + y.*(-274.2290036817964 + ...
        y.*(197.4670779425016 + y.*(-68.5590309679152 + 9.98788038278032.*y))) - 90.6734234051316.*z) + ...
        x.*(-5458.34205214835 - 980.14153344888.*y + ...
        x.*(2247.60742726704 - 340.1237483177863.*x + 220.542973797483.*y) + 180.142097805543.*z) + ...
        z.*(-219.1676534131548 + (-16.32775915649044 - 120.7020447884644.*z).*z)) + ...
        z.*(598.378809221703 + z.*(-156.8822727844005 + (204.1334828179377 - 10.23755797323846.*z).*z)) + ...
        y.*(-1480.222530425046 + z.*(-525.876123559641 + (249.57717834054571 - 88.449193048287.*z).*z) + ...
        y.*(-129.1994027934126 + z.*(1149.174198007428 + z.*(-162.5751787551336 + 76.9195462169742.*z)) + ...
        y.*(-30.0682112585625 - 1380.9597954037708.*z + y.*(2.626801985426835 + 703.695562834065.*z))))) + ...
        y.*(1187.3715515697959 + z.*(1458.233059470092 + ...
        z.*(-687.913805923122 + z.*(249.375342232496 + z.*(-63.313928772146 + 14.09317606630898.*z)))) + ...
        y.*(1760.062705994408 + y.*(-450.535298526802 + ...
        y.*(182.8520895502518 + y.*(-43.3206481750622 + 4.26033941694366.*y) + ...
        z.*(-595.457483974374 + (149.452282277512 - 72.9745838003176.*z).*z)) + ...
        z.*(1388.489628266536 + z.*(-409.779283929806 + (227.123395681188 - 22.2565468652826.*z).*z))) + ...
        z.*(-1721.528607567954 + z.*(674.819060538734 + ...
        z.*(-356.629112415276 + (88.4080716616 - 15.84003094423364.*z).*z)))));
  
g_SA_mod_t = 0.5.*sfac.*g_SA_mod_t;   
    
y_pt = 0.025*pt0;

g_SA_mod_pt = 8645.36753595126 + ...
        x.*(-7296.43987145382 + x.*(8103.20462414788 + ...
        y_pt.*(2175.341332000392 + y_pt.*(-274.2290036817964 + ...
        y_pt.*(197.4670779425016 + y_pt.*(-68.5590309679152 + 9.98788038278032.*y_pt)))) + ...
        x.*(-5458.34205214835 - 980.14153344888.*y_pt + ...
        x.*(2247.60742726704 - 340.1237483177863.*x + 220.542973797483.*y_pt))) + ...
        y_pt.*(-1480.222530425046 + y_pt.*(-129.1994027934126 + ...
        y_pt.*(-30.0682112585625 + y_pt.*2.626801985426835)))) + ...
        y_pt.*(1187.3715515697959 + y_pt.*(1760.062705994408 + y_pt.*(-450.535298526802 + ...
        y_pt.*(182.8520895502518 + y_pt.*(-43.3206481750622 + 4.26033941694366.*y_pt)))));
    
g_SA_mod_pt = 0.5*sfac*g_SA_mod_pt;   

h_SA = g_SA_mod_t - temp_ratio.*g_SA_mod_pt;

if transposed
    h_CT = h_CT.';
    h_SA = h_SA.';
end

end