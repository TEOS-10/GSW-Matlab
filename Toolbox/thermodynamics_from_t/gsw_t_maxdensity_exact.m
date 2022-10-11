function t_maxdensity_exact = gsw_t_maxdensity_exact(SA,p)

% gsw_t_maxdensity_exact                     in-situ temperature of maximum 
%                                                       density of seawater
% =========================================================================
%
% USAGE:
%  t_maxdensity_exact = gsw_t_maxdensity_exact(SA,p)
%
% DESCRIPTION:
%  Calculates the in-situ temperature of maximum density of seawater. 
%  This function returns the in-situ temperature at which the density
%  of seawater is a maximum, at given Absolute Salinity, SA, and sea 
%  pressure, p (in dbar).  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  t_maxdensity_exact  =  in-situ temperature at which            [ deg C ]
%                         the density of seawater is a maximum for
%                         given Absolute Salinity and pressure.
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
%  McDougall, T.J., and S.J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_t_maxdensity_exact:  Requires two inputs')
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
    error('gsw_t_maxdensity_exact: Inputs array dimensions arguments do not agree')
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

dt = 0.001; % the temperature increment for calculating the gibbs_PTT derivative.

t = 3.978 - 0.22072*SA;                    % the initial guess of t_maxden.

gibbs_PTT = 1.1e-8;                          % the initial guess for g_PTT.

for Number_of_iterations = 1:3
    t_old = t;
    gibbs_PT = gsw_gibbs(0,1,1,SA,t_old,p);
    t = t_old - gibbs_PT./gibbs_PTT ; % this is half way through the modified method (McDougall and Wotherspoon, 2013)
    t_mean = 0.5*(t + t_old);
    gibbs_PTT = (gsw_gibbs(0,1,1,SA,t_mean + dt,p) - ...
        gsw_gibbs(0,1,1,SA,t_mean - dt,p))./(dt + dt);
    t = t_old - gibbs_PT./gibbs_PTT;
end

% After three iterations of this modified Newton-Raphson iteration, the 
% error in t_maxdensity_exact is typically no larger than 1x10^-15 degress C.  

t_maxdensity_exact = t;

if transposed
    t_maxdensity_exact = t_maxdensity_exact.';
end

end
