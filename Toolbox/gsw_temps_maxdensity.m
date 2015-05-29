function [t_maxden, pt_maxden, CT_maxden] = gsw_temps_maxdensity(SA,p)

% gsw_temps_maxdensity          temperatures of maximum density of seawater
% =========================================================================
%
% USAGE:
%  [t_maxden, pt_maxden, CT_maxden] = gsw_temps_maxdensity(SA,p)
%
% DESCRIPTION:
%  Calculates the temperatures of maximum density of seawater. This
%  function returns the in-situ, potential, and Conservative temperatures
%  at which the density of seawater is a maximum, at given Absolute
%  Salinity, SA, and sea pressure, p (in dbar).  
%
% INPUT:
%  SA =  Absolute Salinity                                         [ g/kg ]
%  p  =  sea pressure                                              [ dbar ]
%        (ie. absolute pressure - 10.1325 dbar) 
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  t_maxden   =  In-situ temperature at which the density of seawater 
%                is a maximum for given Absolute Salinity and pressure,
%                measured on the ITS-90 temperature scale.        [ deg C ]
%  pt_maxden  =  potential temperature at which the density of seawater 
%                is a maximum for given Absolute Salinity and pressure. 
%                This is the potential temperature referenced to a sea
%                pressure of 0 dbar.                              [ deg C ]
%  CT_maxden  =  Conservative Temperature at which the density of 
%                seawater is a maximum for given Absolute Salinity 
%                and pressure.                                    [ deg C ]
%
% AUTHOR: 
%   Trevor McDougall & Paul Barker [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.42 of this TEOS-10 Manual.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_temps_maxdensity:  Requires two inputs')
end %if

if ~(nargout == 3)
   error('gsw_temps_maxdensity:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mp,np] = size(p);

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_temps_maxdensity: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA';
    p = p';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

n0 = 0; 
n1 = 1;

dt = 0.001; % the temperature increment for calculating the gibbs_PTT derivative.

t = 3.978 - 0.22072*SA;                    % the initial guess of t_maxden.

gibbs_PTT = 1.1e-8;                          % the initial guess for g_PTT.

for Number_of_iterations = 1:3
    t_old = t;
    gibbs_PT = gsw_gibbs(n0,n1,n1,SA,t_old,p);
    t = t_old - gibbs_PT./gibbs_PTT ; % this is half way through the modified method
    t_mean = 0.5*(t + t_old);
    gibbs_PTT = (gsw_gibbs(n0,n1,n1,SA,t_mean + dt,p) - ...
        gsw_gibbs(n0,n1,n1,SA,t_mean - dt,p))./(dt + dt);
    t = t_old - gibbs_PT./gibbs_PTT;
end

% After three iterations of this modified Newton-Raphson iteration, the 
% error in t_maxden is typically no larger than 1x10^-15 degress C.  

t_maxden  = t;

pt_maxden = gsw_pt0_from_t(SA,t_maxden,p);

CT_maxden = gsw_CT_from_pt(SA,pt_maxden);

if transposed
    t_maxden = t_maxden';
    pt_maxden = pt_maxden';
    CT_maxden = CT_maxden';
end

end
