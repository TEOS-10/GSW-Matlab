function [Sstar, in_ocean] = gsw_Sstar_from_SP(SP,p,long,lat)

% gsw_Sstar_from_SP              Preformed Salinity from Practical Salinity
%==========================================================================
%
% USAGE:  
%  [Sstar, in_ocean] = gsw_Sstar_from_SP(SP,p,long,lat)
%
% DESCRIPTION:
%  Calculates Preformed Salinity from Absolute Salinity.
%  Since SP is non-negative by definition, this function changes any 
%  negative input values of SP to be zero.  
%
% INPUT:
%  SP   =  Practical  Salinity  (PSS-78)                       [ unitless ]
%  p    =  sea pressure                                            [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  long =  longitude in decimal degrees                      [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat  =  latitude in decimal degrees north                [ -90 ... +90 ] 
%
%  p, lat and long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where Sstar is MxN.
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
%  David Jackett, Trevor McDougall & Paul Barker       [ help@teos-10.org ]
%
% VERSION NUMBER: 3.0 (27th March, 2011)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.5 and appendices A.4 and A.5 of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm 
%   for estimating Absolute Salinity in the global ocean.  Submitted to 
%   Ocean Science. A preliminary version is available at Ocean Sci. Discuss.,
%   6, 215-242.  
%   http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==4)
   error('gsw_Sstar_from_SP:  Requires four inputs')
end %if

[ms,ns] = size(SP);
[mp,np] = size(p);

if (mp == 1) & (np == 1)               % p is a scalar - fill to size of SP
    p = p*ones(size(SP));
elseif (ns == np) & (mp == 1)          % p is row vector,
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (np == 1)          % p is column vector,
    p = p(:,ones(1,ns));                 % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_Sstar_from_SP: Inputs array dimensions arguments do not agree')
end %if

[mla,nla] = size(lat);

if (mla == 1) & (nla == 1)             % lat is a scalar - fill to size of SP
    lat = lat*ones(size(SP));
elseif (ns == nla) & (mla == 1)        % lat is a row vector,
    lat = lat(ones(1,ms), :);          % copy down each column.
elseif (ms == mla) & (nla == 1)        % lat is a column vector,
    lat = lat(:,ones(1,ns));           % copy across each row.
elseif (ns == mla) & (nla == 1)        % lat is a transposed row vector,
    lat = lat.';                         % transposed then
    lat = lat(ones(1,ms), :);           % copy down each column.
elseif (ms == mla) & (ns == nla)
    % ok
else
    error('gsw_Sstar_from_SP: Inputs array dimensions arguments do not agree')
end %if

[mlo,nlo] = size(long);
[Iwest] =find(long < 0);
if ~isempty(Iwest)
    long(Iwest) = long(Iwest) + 360; 
end
if (mlo == 1) & (nlo == 1)            % long is a scalar - fill to size of SP
    long = long*ones(size(SP));
elseif (ns == nlo) & (mlo == 1)       % long is a row vector,
    long = long(ones(1,ms), :);       % copy down each column.
elseif (ms == mlo) & (nlo == 1)       % long is a column vector,
    long = long(:,ones(1,ns));        % copy across each row.
elseif (ns == mlo) & (nlo == 1)       % long is a transposed row vector,
    long = long.';                      % transposed then
    long = long(ones(1,ms), :);        % copy down each column.
elseif (ms == mlo) & (ns == nlo)
    % ok
else
    error('gsw_Sstar_from_SP: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SP = SP.';
    p = p.';
    lat = lat.';
    long = long.';
    transposed = 1;
else
    transposed = 0;
end

[Iout_of_range] = find(p < 100 & SP > 120);
SP(Iout_of_range) = NaN;
[Iout_of_range] = find(p >= 100 & SP > 42);
SP(Iout_of_range) = NaN;

[Inan] = find(abs(SP) == 99999 | abs(SP) == 999999);
SP(Inan) = NaN;
[Inan] = find(abs(p) == 99999 | abs(p) == 999999);
p(Inan) = NaN;
[Inan] = find(abs(long) == 9999 | abs(long) == 99999);
long(Inan) = NaN;
[Inan] = find(abs(lat) == 9999 | abs(lat) == 99999);
lat(Inan) = NaN;

if ~isempty(find(p < -1.5 | p > 12000))
    error('gsw_Sstar_from_SP: pressure is out of range')
end
if ~isempty(find(long < 0 | long > 360))
    error('gsw_Sstar_from_SP: longitude is out of range')
end
if ~isempty(find(abs(lat) > 90))
    error('gsw_Sstar_from_SP: latitude is out of range')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
  
% These few lines ensure that SP is non-negative.
[I_neg_SP] = find(SP < 0);
if ~isempty(I_neg_SP)
    SP(I_neg_SP) = 0;
end

r_1 = 0.35;

[Iocean] = find(~isnan(SP.*p.*lat.*long));

Sstar = nan(size(SP));
SAAR = nan(size(SP));
in_ocean = nan(size(SP));

[SAAR(Iocean), in_ocean(Iocean)] = gsw_SAAR(p(Iocean),long(Iocean),lat(Iocean));

Sstar(Iocean) = (35.16504/35)*SP(Iocean).*(1 - r_1.*SAAR(Iocean));

%In the Baltic Sea, Sstar = SA.

Sstar_baltic(Iocean) = gsw_SA_from_SP_Baltic(SP(Iocean),long(Iocean),lat(Iocean));

Ibaltic = find(~isnan(Sstar_baltic(Iocean)));

Sstar(Iocean(Ibaltic)) = Sstar_baltic(Iocean(Ibaltic));

if transposed
    Sstar = Sstar.';
    in_ocean = in_ocean.';
end

end
