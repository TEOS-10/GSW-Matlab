function [h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives(SA,CT,p)

% gsw_enthalpy_second_derivatives            second derivatives of enthalpy
%                                                        (48-term equation)
% =========================================================================
%
% USAGE:
%  [h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of specific
%  enthalpy (h),using the computationally-efficient 48-term expression for 
%  density in terms of SA, CT and p (IOC et al., 2010).
%   (1) h_SA_SA, second-order derivative with respect to Absolute Salinity 
%       at constant CT & p.
%   (2) h_SA_CT, second-order derivative with respect to SA & CT at 
%       constant p. 
%   (3) h_CT_CT, second-order derivative with respect to CT at constant SA 
%       and p. 
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2010).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
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
%  h_SA_SA  =  The second derivative of specific enthalpy with respect to 
%              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
%  h_SA_CT  =  The second derivative of specific enthalpy with respect to 
%              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
%  h_CT_CT  =  The second derivative of specific enthalpy with respect to 
%              CT at constant SA and p.                      [ J/(kg K^2) ]
%
% AUTHOR:   
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (26th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%   
%  McDougall, T. J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%    See Eqns. (18) and (22)
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_enthalpy_second_derivatives:  Requires three inputs')
end %if

if ~(nargout == 3)
   error('gsw_enthalpy_second_derivatives:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_enthalpy_second_derivatives: SA and CT do not have the same dimensions')
end %if

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
    error('gsw_enthalpy_second_derivatives: The dimensions of p do not agree')
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

SA(SA<0) = 0;

dSA = 5e-2;                % increment in Absolute Salinity is 5e-2 g kg^-1

SA_l = nan(size(SA));
SA_l(SA >= dSA) = SA(SA >= dSA) - dSA;
SA_l(SA<dSA) = 0;
SA_u = SA + dSA;  

dCT = 5e-3;          % increment in Conservative Temperature is 5e-3 deg C.
CT_l = CT - dCT;
CT_u = CT + dCT;

[hSA_l,dummy] = gsw_enthalpy_first_derivatives(SA_l,CT,p);
[hSA_u,dummy] = gsw_enthalpy_first_derivatives(SA_u,CT,p);

h_SA_SA = (hSA_u - hSA_l)./(SA_u - SA_l);

if (any(SA<1))
   [h_SA_SA_exact,dummy,dummy] = gsw_enthalpy_second_derivatives_CT_exact(SA,CT,p);
   
   h_SA_SA_48 = h_SA_SA;
   
   SA_taper = 2 - 2*SA(SA>0.5 & SA<1);

   h_SA_SA(SA>0.5 & SA<1) = SA_taper.*(h_SA_SA_exact(SA>0.5 & SA<1)) + (1-SA_taper).*h_SA_SA_48(SA>0.5 & SA<1);

   h_SA_SA(SA<=0.5) = h_SA_SA_exact(SA<=0.5);
end

[hSA_l,dummy] = gsw_enthalpy_first_derivatives(SA,CT_l,p);
[hSA_u,dummy] = gsw_enthalpy_first_derivatives(SA,CT_u,p);

h_SA_CT  = (hSA_u - hSA_l)./(CT_u - CT_l);

hCT_l = gsw_dynamic_enthalpy(SA,CT_l,p);
hCT = gsw_dynamic_enthalpy(SA,CT,p);
hCT_u = gsw_dynamic_enthalpy(SA,CT_u,p);

h_CT_CT = 4.*(hCT_u - 2.*hCT + hCT_l)./((CT_u - CT_l).^2);

if transposed
    h_SA_SA = h_SA_SA.';
    h_SA_CT = h_SA_CT.';
    h_CT_CT = h_CT_CT.';
end

end
