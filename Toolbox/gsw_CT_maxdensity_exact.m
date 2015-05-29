function CT_maxdensity_exact = gsw_CT_maxdensity_exact(SA,p)

% gsw_CT_maxdensity_exact               Conservative Temperature of maximum 
%                                                       density of seawater
% =========================================================================
%
% USAGE:
%  CT_maxdensity_exact = gsw_CT_maxdensity_exact(SA,p)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of maximum density of seawater. 
%  This function returns the Conservative temperature at which the density
%  of seawater is a maximum, at given Absolute Salinity, SA, and sea 
%  pressure, p (in dbar).  
%
% INPUT:
%  SA =  Absolute Salinity                                         [ g/kg ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  CT_maxdensity_exact  =  Conservative Temperature at which      [ deg C ]
%                          the density of seawater is a maximum for
%                          given Absolute Salinity and pressure.
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
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
   error('gsw_CT_maxdensity_exact:  Requires two inputs')
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
    error('gsw_CT_maxdensity_exact: Inputs array dimensions arguments do not agree')
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

t_maxdensity_exact = gsw_t_maxdensity_exact(SA,p);
CT_maxdensity_exact = gsw_CT_from_t(SA,t_maxdensity_exact,p);

if transposed
    CT_maxdensity_exact = CT_maxdensity_exact.';
end

end
