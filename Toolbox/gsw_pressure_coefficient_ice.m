function pressure_coefficient_ice = gsw_pressure_coefficient_ice(t,p)

% gsw_pressure_coefficient_ice                  pressure coefficient of ice
%==========================================================================
%
% USAGE:  
%  pressure_coefficient_ice = gsw_pressure_coefficient_ice(t,p)
%
% DESCRIPTION:
%  Calculates pressure coefficient of ice. 
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where t is MxN.
%
% OUTPUT:
%  pressure_coefficient_ice  =  pressure coefficient of ice          [Pa/K]
%   Note. The output units are Pa/K NOT dbar/K.
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%   See Eqn. (2.15.1) of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_pressure_coefficient_ice:  Requires two inputs')
end

[mt,nt] = size(t);
[mp,np] = size(p);

if (mp == 1) & (np == 1)              % p scalar - fill to size of t
    p = p*ones(size(t));
elseif (nt == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,mt), :);              % copy down each column.
elseif (mt == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,nt));               % copy across each row.
elseif (nt == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,mt), :);                % copy down each column.
elseif (mt == mp) & (nt == np)
    % ok
else
    error('gsw_pressure_coefficient_ice: Inputs array dimensions arguments do not agree')
end

if mt == 1
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

pressure_coefficient_ice = -gsw_gibbs_ice(1,1,t,p)./gsw_gibbs_ice(0,2,t,p);

if transposed
    pressure_coefficient_ice = pressure_coefficient_ice.';
end

end
