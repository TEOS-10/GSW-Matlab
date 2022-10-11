function Fdelta = gsw_Fdelta(p,long,lat)

% gsw_Fdelta                                                         Fdelta
%==========================================================================
%
% USAGE:  
%  Fdelta = gsw_Fdelta(p,long,lat)
%
% DESCRIPTION:
%  Calculates Fdelta from the Absolute Salinity Anomaly Ratio (SAAR).  It
%  finds SAAR by calling the function "gsw_SAAR(p,long,lat)" and then 
%  simply calculates Fdelta from 
% 
%     Fdelta = (1 + r1)SAAR/(1 - r1*SAAR)
%            = (SA/Sstar) - 1 
% 
%  with r1 being the constant 0.35 based on the work of Pawlowicz et al.
%  (2011). Note that since SAAR is everywhere less than 0.001 in the global
%  ocean, Fdelta is only slighty different to 1.35*SAAR. 
%
% INPUT:
%  p    =  sea pressure                                            [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  long =  longitude in decimal degrees                      [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat  =  latitude in decimal degrees north                [ -90 ... +90 ] 
%
%  lat & long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where p is MxN.
%
% OUTPUT:
%  Fdelta  =  ratio of SA to Sstar, minus 1                    [ unitless ]
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
%    See section 2.5 and appendices A.4 and A.5 of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and 
%   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
%   Ocean Science, 8, 1123-1134.  
%   http://www.ocean-sci.net/8/1123/2012/os-8-1123-2012.pdf 
%
% Pawlawicz, R., D.G. Wright and F.J. Millero, 2011; The effects of
%   biogeochemical processes on oceanic conductivity/salinty/density
%   relationships and the characterization of real seawater. Ocean Science,
%   7, 363-387.  http://www.ocean-sci.net/7/363/2011/os-7-363-2011.pdf
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
    error('gsw_Fdelta:  Requires three inputs')
end %if

[mp,np] = size(p);
[mla,nla] = size(lat);

if (mla == 1) & (nla == 1)             % lat is a scalar - fill to size of p
    lat = lat*ones(size(p));
elseif (np == nla) & (mla == 1)        % lat is a row vector,
    lat = lat(ones(1,mp), :);           % copy down each column.
elseif (mp == mla) & (nla == 1)        % lat is a column vector,
    lat = lat(:,ones(1,np));            % copy across each row.
elseif (np == mla) & (nla == 1)        % lat is a transposed row vector,
    lat = lat.';                         % transposed then
    lat = lat(ones(1,ms), :);           % copy down each column.
elseif (mp == mla) & (np == nla)
    % ok
else
    error('gsw_Fdelta: Inputs array dimensions arguments do not agree')
end %if

[mlo,nlo] = size(long);
long(long < 0) = long(long < 0) + 360; 

if (mlo == 1) & (nlo == 1)            % long is a scalar - fill to size of p
    long = long*ones(size(p));
elseif (np == nlo) & (mlo == 1)       % long is a row vector,
    long = long(ones(1,mp), :);        % copy down each column.
elseif (mp == mlo) & (nlo == 1)       % long is a column vector,
    long = long(:,ones(1,np));         % copy across each row. 
elseif (np == mlo) & (nlo == 1)       % long is a transposed row vector,
    long = long.';                      % transposed then
    long = long(ones(1,mp), :);        % copy down each column.
elseif (ms == nlo) & (mlo == 1)       % long is a transposed column vector,
    long = long.';                      % transposed then
    long = long(:,ones(1,np));        % copy down each column.
elseif (mp == mlo) & (np == nlo)
    % ok
else
    error('gsw_Fdelta: Inputs array dimensions arguments do not agree')
end %if

if mp == 1
    p = p.';
    lat = lat.';
    long = long.';
    transposed = 1;
else
    transposed = 0;
end

% change standard blank fill values to NaN's.
p(abs(p) == 99999 | abs(p) == 999999) = NaN;
long(abs(long) == 9999 | abs(long) == 99999) = NaN;
lat(abs(lat) == 9999 | abs(lat) == 99999) = NaN;

if any(find(p < -1.5 | p > 12000))
    error('gsw_Fdelta: pressure is out of range')
end
if any(find(long < 0 | long > 360))
    error('gsw_Fdelta: longitude is out of range')
end
if any(find(abs(lat) > 90))
    error('gsw_Fdelta: latitude is out of range')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
 
r_1 = 0.35;

SAAR = nan(size(p));

if any(~isnan(p + long + lat))
    [I] = find(~isnan(p + long + lat));
    SAAR(I) = gsw_SAAR(p(I),long(I),lat(I));
end

Fdelta = ((1 + r_1).*SAAR)./(1 - r_1*SAAR);

if transposed
    Fdelta = Fdelta.';
end

end
