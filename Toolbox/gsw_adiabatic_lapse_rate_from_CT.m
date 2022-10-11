function adiabatic_lapse_rate = gsw_adiabatic_lapse_rate_from_CT(SA,CT,p)

% gsw_adiabatic_lapse_rate_from_CT                     adiabatic lapse rate
%==========================================================================
%
% USAGE:
%  adiabatic_lapse_rate = gsw_adiabatic_lapse_rate_from_CT(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the adiabatic lapse rate of sea water from Conservative
%  Temperature.
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
%  adiabatic_lapse_rate  =  adiabatic lapse rate                   [ K/Pa ]
%    Note.  The output is in unit of degress Celsius per Pa,
%      (or equivilently K/Pa) not in units of K/dbar. 
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.22.1) of this TEOS-10 Manual.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3) 
   error('gsw_adiabatic_lapse_rate_from_CT:  Requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_adiabatic_lapse_rate_from_CT: SA and CT must have same dimensions')
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
    error('gsw_adiabatic_lapse_rate_from_CT: Inputs array dimensions arguments do not agree')
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

t = gsw_t_from_CT(SA,CT,p);

adiabatic_lapse_rate = -gsw_gibbs(0,1,1,SA,t,p)./(gsw_gibbs(0,2,0,SA,t,p));

if transposed
    adiabatic_lapse_rate = adiabatic_lapse_rate.';
end

end
