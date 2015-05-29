function dynamic_enthalpy_CT = gsw_dynamic_enthalpy_CT(SA,CT,p)

% gsw_dynamic_enthalpy_CT                      dynamic enthalpy of seawater
%                                                        (48-term equation)
%==========================================================================
%
% USAGE:
%  dynamic_enthalpy_CT = gsw_dynamic_enthalpy_CT(SA,CT,p), or equivalently
%     dynamic_enthalpy = gsw_dynamic_enthalpy(SA,CT,p)
% 
%  Note that gsw_dynamic_enthalpy(SA,CT,p) is identical to 
%  gsw_dynamic_enthalpy_CT(SA,CT,p).  The extra "_CT" emphasises that the 
%  input temperature is Conservative Temperature, but the extra "_CT" part
%  of the function name is not needed. 
%
% DESCRIPTION:
%  Calculates dynamic enthalpy of seawater using the computationally-
%  efficient 48-term expression for density in terms of SA, CT and p
%  (IOC et al., 2010).  Dynamic enthalpy is defined as enthalpy minus
%  potential enthalpy (Young, 2010). 
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
%  dynamic_enthalpy_CT  =  dynamic enthalpy                        [ J/kg ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.2 of this TEOS-10 Manual. 
%
%  McDougall, T. J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%    See Eqns. (18) and (22)
%
%  Young, W.R., 2010: Dynamic enthalpy, Conservative Temperature, and the
%   seawater Boussinesq approximation. Journal of Physical Oceanography, 
%   40, 394-400.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_dynamic_enthalpy_CT: requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_dynamic_enthalpy_CT: SA and CT must have same dimensions')
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
    error('gsw_dynamic_enthalpy_CT: Inputs array dimensions arguments do not agree')
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

dynamic_enthalpy_CT = gsw_dynamic_enthalpy(SA,CT,p);

if transposed
    dynamic_enthalpy_CT = dynamic_enthalpy_CT.';
end

end
