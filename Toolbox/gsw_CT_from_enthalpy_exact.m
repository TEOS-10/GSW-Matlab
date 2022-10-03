function CT = gsw_CT_from_enthalpy_exact(SA,h,p)

% gsw_CT_from_enthalpy_exact         Conservative Temperature from specific
%                                                      enthalpy of seawater  
%==========================================================================
%
% USAGE:
%  CT = gsw_CT_from_enthalpy_exact(SA,h,p)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of seawater, given the Absolute 
%  Salinity, SA, specific enthalpy, h, and pressure p.  The specific 
%  enthalpy input is calculated from the full Gibbs function of seawater,
%  gsw_enthalpy_t_exact. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_CT_from_enthalpy(SA,h,p), which uses the computationally 
%  efficient 75-term expression for specific volume in terms of SA, CT  
%  and p (Roquet et al., 2015).   
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  h   =  specific enthalpy                                        [ J/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & h need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & h are MxN.
%
% OUTPUT:
%  CT  =  Conservative Temperature ( ITS-90)                      [ deg C ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%
%  McDougall T.J. and S.J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_CT_from_enthalpy_exact: requires three inputs')
end

[ms,ns] = size(SA);
[mh,nh] = size(h); 
[mp,np] = size(p);

if (mh ~= ms | nh ~= ns)
    error('gsw_CT_from_enthalpy_exact: SA and h must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(ms,ns);
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
    error('gsw_CT_from_enthalpy_exact: Inputs array dimensions do not agree; check p')
end %if

if ms == 1
    SA = SA.';
    h = h.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensures that SA is non-negative

CT_freezing = gsw_CT_freezing(SA,p,0); % This is the CT freezing temperature
h_below_freeze = gsw_cp0;  %This allows for water to be approx. 1C below the freezing temperature
h_freezing = gsw_enthalpy_CT_exact(SA,CT_freezing,p);
if any(h(:) < (h_freezing(:) - h_below_freeze)) %The input, seawater enthalpy h, is less than the enthalpy at the freezing temperature, i.e. the water is frozen
    [I] = find(h < h_freezing);
    SA(I) = NaN;
end 

h_40 = gsw_enthalpy_CT_exact(SA,40*ones(size(SA)),p);
if any(h(:) > h_40(:)) % The input, seawater enthalpy h, is greater than the enthalpy when CT is 40 C
    [I] = find(h > h_40);
    SA(I) = NaN;
end 

CT = CT_freezing + (40 - CT_freezing).*(h - h_freezing)./(h_40 - h_freezing); % First guess of CT
[dummy, h_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT,p);

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure 
%--------------------------------------------------------------------------

CT_old = CT;
f = gsw_enthalpy_CT_exact(SA,CT_old,p) - h;
CT = CT_old - f./h_CT ; % this is half way through the modified Newton's method (McDougall and Wotherspoon, 2013)
CT_mean = 0.5*(CT + CT_old);
[dummy, h_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT_mean,p);
CT = CT_old - f./h_CT ; % this is the end of one full iteration of the modified Newton's method

% After 1 iteration of this modified Newton-Raphson iteration, the error in
% CT is no larger than 5x10^-14 degrees C, which is machine precision for 
% this calculation. 

if transposed
    CT = CT.';
end

end
