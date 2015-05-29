function [Sstar, in_ocean] = gsw_Sstar_from_SA(SA,p,long,lat)

% gsw_Sstar_from_SA               Preformed Salinity from Absolute Salinity
%==========================================================================
%
% USAGE:
%  [Sstar, in_ocean] = gsw_Sstar_from_SA(SA,p,long,lat)
%
% DESCRIPTION:
%  Converts Preformed Salinity from Absolute Salinity. 
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  p, lat and long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where SA is MxN.
%
% OUTPUT:
%  Sstar     =  Preformed Salinity                                 [ g/kg ]
%  in_ocean  =  0, if long and lat are a long way from the ocean 
%            =  1, if long and lat are in the ocean
%  Note. This flag is only set when the observation is well and truly on
%    dry land; often the warning flag is not set until one is several 
%    hundred kilometres inland from the coast. 
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker    [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 3.0 (27th March, 2011)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_Sstar_from_SA:  Requires four inputs')
end %if

[ms,ns] = size(SA);
[mp,np] = size(p);

if (mp == 1) & (np == 1)                 % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)            % p is row vector,
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (np == 1)            % p is column vector,
    p = p(:,ones(1,ns));                 % copy across each row.
elseif (ns == mp) & (np == 1)            % p is a transposed row vector,
    p = p.';                              % transpose, then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_Sstar_from_SA: Inputs array dimensions arguments do not agree')
end %if

[mla,nla] = size(lat);

if (mla == 1) & (nla == 1)             % lat is a scalar - fill to size of SA
    lat = lat*ones(size(SA));
elseif (ns == nla) & (mla == 1)        % lat is a row vector,
    lat = lat(ones(1,ms), :);          % copy down each column.
elseif (ms == mla) & (nla == 1)        % lat is a column vector,
    lat = lat(:,ones(1,ns));           % copy across each row.
elseif (ns == mla) & (nla == 1)        % lat is a transposed row vector,
    lat = lat.';                        % transpose, then
    lat = lat(ones(1,ms), :);          % copy down each column.
elseif (ms == mla) & (ns == nla)
    % ok
else
    error('gsw_Sstar_from_SA: Inputs array dimensions arguments do not agree')
end %if

[mlo,nlo] = size(long);
[Iwest] =find(long < 0);
if ~isempty(Iwest)
    long(Iwest) = long(Iwest) + 360; 
end
if (mlo == 1) & (nlo == 1)            % long is a scalar - fill to size of SA
    long = long*ones(size(SA));
elseif (ns == nlo) & (mlo == 1)       % long is a row vector,
    long = long(ones(1,ms), :);       % copy down each column.
elseif (ms == mlo) & (nlo == 1)       % long is a column vector,
    long = long(:,ones(1,ns));        % copy across each row.
elseif (ns == mlo) & (nlo == 1)       % long is a transposed row vector,
    long = long.';                     % transpose, then
    long = long(ones(1,ms), :);       % copy down each column.
elseif (ms == mlo) & (ns == nlo)
    % ok
else
    error('gsw_Sstar_from_SA: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    p = p.';
    lat = lat.';
    long = long.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

r_1 = 0.35;

[Iocean] = find(~isnan(SA.*p.*lat.*long));

Sstar = nan(size(SA));
SAAR = nan(size(SA));
in_ocean = nan(size(SA));

[SAAR(Iocean), in_ocean(Iocean)] = gsw_SAAR(p(Iocean),long(Iocean),lat(Iocean));

Sstar(Iocean) = SA(Iocean).*(1 - r_1*SAAR(Iocean))./(1 + SAAR(Iocean));

% In the Baltic Sea, Sstar = SA, and note that gsw_delta_SA returns zero
% for dSA in the Baltic.

if transposed
    Sstar = Sstar.';
    in_ocean = in_ocean.';
end

end
