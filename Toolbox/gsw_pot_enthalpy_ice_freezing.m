function pot_enthalpy_ice_freezing = gsw_pot_enthalpy_ice_freezing(SA,p,saturation_fraction)

% gsw_pot_enthalpy_ice_freezing_poly              potential enthalpy of ice
%                                                 at which seawater freezes
%==========================================================================
%
% USAGE:
%  pot_enthalpy_ice_freezing = ...
%                   gsw_pot_enthalpy_ice_freezing(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice at which seawater freezes.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in 
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default 
%    is 0, air free) 
%
%  p & saturation_fraction (if provided) may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  pot_enthalpy_ice_freezing = potential enthalpy of ice at freezing 
%                              of seawater                         [ J/kg ]
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.  
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3) 
   error('gsw_pot_enthalpy_ice_freezing:  Requires either two or three inputs') 
end

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_pot_enthalpy_ice_freezing: saturation fraction MUST be between zero and one.')
end

[ms,ns] = size(SA);
[mp,np] = size(p);
[msf,nsf] = size(saturation_fraction);

if (mp == 1) & (np == 1) 
    p = p*ones(ms,ns);
elseif (ns == np) & (mp == 1)
    p = p(ones(1,ms), :);
elseif (ms == mp) & (np == 1)
    p = p(:,ones(1,ns));
elseif (ns == mp) & (np == 1)
    p = p.';
    p = p(ones(1,ms), :);
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_pot_enthalpy_ice_freezing: Inputs array dimensions arguments do not agree')
end

if (msf == 1) & (nsf == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(SA));         % fill to size of SA
elseif (ns == nsf) & (msf == 1)                        % saturation_fraction is row vector,
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (nsf == 1)                     % saturation_fraction is column vector,
    saturation_fraction = saturation_fraction(:,ones(1,ns));        % copy across each row.
elseif (ns == msf) & (nsf == 1)           % saturation_fraction is a transposed row vector,
    saturation_fraction = saturation_fraction.';                           % transposed then
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (ns == nsf)
    % ok
else
    error('gsw_pot_enthalpy_ice_freezing: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

t_freezing = gsw_t_freezing(SA,p,saturation_fraction); 

pt0_ice = gsw_pt0_from_t_ice(t_freezing,p); 

pot_enthalpy_ice_freezing = gsw_pot_enthalpy_from_pt_ice(pt0_ice);

% set any values that are out of range to be NaN. 
pot_enthalpy_ice_freezing(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if transposed
    pot_enthalpy_ice_freezing = pot_enthalpy_ice_freezing.';
end

end