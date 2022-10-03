function CT_maxdensity = gsw_CT_maxdensity(SA,p)

% gsw_CT_maxdensity                     Conservative Temperature of maximum 
%                                    density of seawater (75-term equation)
% =========================================================================
%
% USAGE:
%  CT_maxdensity = gsw_CT_maxdensity(SA,p)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of maximum density of seawater. 
%  This function returns the Conservative temperature at which the density
%  of seawater is a maximum, at given Absolute Salinity, SA, and sea 
%  pressure, p (in dbar).  This function uses the computationally-efficient
%  expression for specific volume in terms of SA, CT and p
%  (Roquet et al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA =  Absolute Salinity                                         [ g/kg ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  CT_maxdensity  =  Conservative Temperature at which            [ deg C ]
%                    the density of seawater is a maximum for
%                    given Absolute Salinity and pressure.
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.42 of this TEOS-10 Manual.  
%
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of 
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

if ~(nargin == 2)
   error('gsw_CT_maxdensity:  Requires two inputs')
end %if

[ms,ns] = size(SA);
[mp,np] = size(p);

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
    error('gsw_CT_maxdensity: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

dCT = 0.001;                      % the Conservative Temperature increment.

CT = 3.978 - 0.22072*SA;                         % the initial guess of CT.

dalpha_dCT = 1.1e-5;                  % the initial guess for d(alpha)_dCT.

for Number_of_iterations = 1:3
    CT_old = CT;
    alpha = gsw_alpha(SA,CT_old,p);
    CT = CT_old - alpha./dalpha_dCT; % this is half way through the modified method
    CT_mean = 0.5*(CT + CT_old);
    dalpha_dCT = (gsw_alpha(SA,CT_mean + dCT,p) ...
                  - gsw_alpha(SA,CT_mean - dCT,p))./(dCT + dCT);
    CT = CT_old - alpha./dalpha_dCT;
end

% After three iterations of this modified Newton-Raphson (McDougall and 
% Wotherspoon, 2013) iteration, the error in CT_maxdensity is typically no
% larger than 1x10^-15 degress C.  

CT_maxdensity = CT;

if transposed
    CT_maxdensity = CT_maxdensity.';
end

end
