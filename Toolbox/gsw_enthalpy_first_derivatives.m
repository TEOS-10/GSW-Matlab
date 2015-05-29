function [h_SA, h_CT] = gsw_enthalpy_first_derivatives(SA,CT,p)

% gsw_enthalpy_first_derivatives              first derivatives of enthalpy
%==========================================================================
%
% USAGE:
%  [h_SA, h_CT] = gsw_enthalpy_first_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following two derivatives of specific enthalpy (h) of
%  seawater using the computationally-efficient 48-term expression for 
%  density in terms of SA, CT and p (IOC et al., 2010).  
%   (1) h_SA, the derivative with respect to Absolute Salinity at 
%       constant CT and p, and
%   (2) h_CT, derivative with respect to CT at constant SA and p. 
%  Note that h_P is specific volume (1/rho) it can be calulated by calling
%  gsw_specvol(SA,CT,p).
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
%  h_SA  =  The first derivative of specific enthalpy with respect to 
%           Absolute Salinity at constant CT and p.     
%                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
%  h_CT  =  The first derivative of specific enthalpy with respect to 
%           CT at constant SA and p.                           [ J/(kg K) ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.03 (20th May, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%    See Eqns. (A.11.18), (A.11.15) and (A.11.12) of this TEOS-10 Manual.   
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
   error('gsw_enthalpy_first_derivatives:  Requires three inputs')
end %if

if ~(nargout == 2)
   error('gsw_enthalpy_first_derivatives:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
    error('gsw_enthalpy_first_derivatives: SA and CT do not have the same dimensions')
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
    error('gsw_enthalpy_first_derivatives: The dimensions of p do not agree')
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
dSA = 0.01;                % increment in Absolute Salinity is 1e-2 g kg^-1
SA_l = nan(size(SA));
SA_l(SA >= dSA) = SA(SA >= dSA) - dSA;
SA_l(SA<dSA) = 0;
SA_u = SA + dSA;

hSA_l = gsw_dynamic_enthalpy(SA_l,CT,p);
hSA_u = gsw_dynamic_enthalpy(SA_u,CT,p);

h_SA = (hSA_u - hSA_l)./(SA_u - SA_l);

dCT = 0.005;          % increment in Conservative Temperature is 5e-3 deg C.
CT_l = CT - dCT;
CT_u = CT + dCT;

hCT_l = gsw_dynamic_enthalpy(SA,CT_l,p);
hCT_u = gsw_dynamic_enthalpy(SA,CT_u,p);

h_CT = gsw_cp0 + (hCT_u - hCT_l)./(CT_u - CT_l);

if transposed
    h_CT = h_CT.';
    h_SA = h_SA.';
end

end