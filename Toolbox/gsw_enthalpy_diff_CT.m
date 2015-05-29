function enthalpy_diff_CT = gsw_enthalpy_diff_CT(SA,CT,p_shallow,p_deep)

% gsw_enthalpy_diff_CT              difference of enthalpy at two pressures
%                                                        (48-term equation)
%==========================================================================
%
% USAGE:
%  enthalpy_diff_CT = gsw_enthalpy_diff_CT(SA,CT,p_shallow,p_deep), or 
%  equivalently enthalpy_diff = gsw_enthalpy_diff(SA,CT,p_shallow,p_deep)
% 
%  Note that gsw_enthalpy_diff(SA,CT,p_shallow,p_deep) is identical to 
%  gsw_enthalpy_diff_CT(SA,CT,p_shallow,p_deep).  The extra "_CT" 
%  emphasises that the input temperature is Conservative Temperature, but 
%  the extra "_CT" part of the function name is not needed. 
%
% DESCRIPTION:
%  Calculates the difference of the specific enthalpy of seawater between 
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  This 
%  function uses the computationally-efficient 48-term expression for
%  density in terms of SA, CT and p (IOC et al., 2010).  The output
%  (enthalpy_diff_CT) is the specific enthalpy evaluated at (SA,CT,p_deep)
%  minus the specific enthalpy at (SA,CT,p_shallow). 
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2010).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA         =  Absolute Salinity                                 [ g/kg ]
%  CT         =  Conservative Temperature (ITS-90)                [ deg C ]
%  p_shallow  =  upper sea pressure                                [ dbar ]
%                ( i.e. shallower absolute pressure - 10.1325 dbar ) 
%  p_deep     =  lower sea pressure                                [ dbar ]
%                ( i.e. deeper absolute pressure - 10.1325 dbar )
%
%  p_shallow and p_deep may have dimensions Mx1 or 1xN or MxN, 
%  where SA and CT are MxN.
%
% OUTPUT:
%  enthalpy_diff_CT  =  difference of specific enthalpy            [ J/kg ]
%                       (deep minus shallow)
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker.                     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) and (A.30.6) of this TEOS-10 Manual. 
%
%  McDougall, T. J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%    See Eqns. (18) and (22)
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
    error('gsw_enthalpy_diff_CT: requires four inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mpu,npu] = size(p_shallow);
[mpl,npl] = size(p_deep);

if (ms~=mt) | (ns~=nt)
    error('gsw_enthalpy_diff_CT: SA & CT need to have the same dimensions')
end

if (mpu == 1) & (npu == 1)                           % p_shallow is a scalar 
    p_shallow = p_shallow*ones(size(SA));
elseif (ns == npu) & (mpu == 1)                  % p_shallow is row vector,
    p_shallow = p_shallow(ones(1,ms), :);          % copy down each column.
elseif (ms == mpu) & (npu == 1)               % p_shallow is column vector,
    p_shallow = p_shallow(:,ones(1,ns));            % copy across each row.
elseif (ns == mpu) & (npu == 1)          % p_shallow is a transposed row vector,
    p_shallow = p_shallow.';                              % transposed then
    p_shallow = p_shallow(ones(1,ms), :);                % copy down each column.
elseif (ms == mpu) & (ns == npu)
    % ok
end

if (mpl == 1) & (npl == 1)                             % p_deep is a scalar  
    p_deep = p_deep*ones(size(SA));
elseif (ns == npl) & (mpl == 1)                     % p_deep is row vector,
    p_deep = p_deep(ones(1,ms), :);                % copy down each column.
elseif (ms == mpl) & (npl == 1)                  % p_deep is column vector,
    p_deep = p_deep(:,ones(1,ns));                  % copy across each row.
elseif (ns == mpl) & (npl == 1)          % p_deep is a transposed row vector,
    p_deep = p_deep.';                              % transposed then
    p_deep = p_deep(ones(1,ms), :);                % copy down each column.
elseif (ms == mpl) & (ns == npl)
    % ok
else
    error('gsw_enthalpy_diff_CT: Inputs array dimensions arguments do not agree')
end %if

if any(p_shallow - p_deep > 0);
    error('gsw_enthalpy_diff_CT: p_deep needs to be greater than or equal to p_shallow')
end
    
if ms == 1
    SA = SA.';
    CT = CT.';
    p_shallow = p_shallow.';
    p_deep = p_deep.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

enthalpy_diff_CT = gsw_enthalpy_diff(SA,CT,p_shallow,p_deep);

if transposed
    enthalpy_diff_CT = enthalpy_diff_CT.';
end

end
