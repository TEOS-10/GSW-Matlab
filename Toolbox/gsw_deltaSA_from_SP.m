function deltaSA  = gsw_deltaSA_from_SP(SP,p,long,lat)

% gsw_deltaSA_from_SP                             Absolute Salinity Anomaly
%                                                   from Practical Salinity
%==========================================================================
%
% USAGE:  
%  deltaSA = gsw_deltaSA_from_SP(SP,p,long,lat)
%
% DESCRIPTION:
%  Calculates Absolute Salinity Anomaly from Practical Salinity.  Since SP
%  is non-negative by definition, this function changes any negative input 
%  values of SP to be zero.  
%
% INPUT:
%  SP   =  Practical Salinity  (PSS-78)                        [ unitless ]
%  p    =  sea pressure                                            [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  long =  longitude in decimal degrees                      [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat  =  latitude in decimal degrees north                [ -90 ... +90 ] 
%
%  p, lat & long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where SP is MxN.
%
% OUTPUT:
%  deltaSA  =  Absolute Salinity Anomaly                           [ g/kg ]
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
%   Ocean Science, 8, 1117-1128.  
%   http://www.ocean-sci.net/8/1117/2012/os-8-1117-2012.pdf 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==4)
    error('gsw_deltaSA_from_SP:  Requires four inputs')
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
    error('gsw_deltaSA_from_SP: Inputs array dimensions arguments do not agree')
end %if

[mla,nla] = size(lat);

if (mla == 1) & (nla == 1)             % lat is a scalar - fill to size of SP
    lat = lat*ones(size(SP));
elseif (ns == nla) & (mla == 1)        % lat is a row vector,
    lat = lat(ones(1,ms), :);           % copy down each column.
elseif (ms == mla) & (nla == 1)        % lat is a column vector,
    lat = lat(:,ones(1,ns));            % copy across each row.
elseif (ns == mla) & (nla == 1)        % lat is a transposed row vector,
    lat = lat.';                         % transposed then
    lat = lat(ones(1,ms), :);           % copy down each column.
elseif (ms == mla) & (ns == nla)
    % ok
else
    error('gsw_deltaSA_from_SP: Inputs array dimensions arguments do not agree')
end %if

[mlo,nlo] = size(long);
long(long < 0) = long(long < 0) + 360; 

if (mlo == 1) & (nlo == 1)            % long is a scalar - fill to size of SP
    long = long*ones(size(SP));
elseif (ns == nlo) & (mlo == 1)       % long is a row vector,
    long = long(ones(1,ms), :);        % copy down each column.
elseif (ms == mlo) & (nlo == 1)       % long is a column vector,
    long = long(:,ones(1,ns));         % copy across each row. 
elseif (ns == mlo) & (nlo == 1)       % long is a transposed row vector,
    long = long.';                      % transposed then
    long = long(ones(1,ms), :);        % copy down each column.
elseif (ms == nlo) & (mlo == 1)       % long is a transposed column vector,
    long = long.';                      % transposed then
    long = long(:,ones(1,ns));        % copy down each column.
elseif (ms == mlo) & (ns == nlo)
    % ok
else
    error('gsw_deltaSA_from_SP: Inputs array dimensions arguments do not agree')
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

% remove out of range values.
SP(p < 100 & SP > 120) = NaN;
SP(p >= 100 & SP > 42) = NaN;

% change standard blank fill values to NaN's.
SP(abs(SP) == 99999 | abs(SP) == 999999) = NaN;
p(abs(p) == 99999 | abs(p) == 999999) = NaN;
long(abs(long) == 9999 | abs(long) == 99999) = NaN;
lat(abs(lat) == 9999 | abs(lat) == 99999) = NaN;

if any(p < -1.5 | p > 12000)
    error('gsw_deltaSA_from_SP: pressure is out of range')
end
if any(long < 0 | long > 360)
    error('gsw_deltaSA_from_SP: longitude is out of range')
end
if any(abs(lat) > 90)
    error('gsw_deltaSA_from_SP: latitude is out of range')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
 
% This ensures that SP is non-negative.
SP(SP < 0) = 0;

SA = gsw_SA_from_SP(SP,p,long,lat);
SR = gsw_SR_from_SP(SP);
deltaSA = SA - SR;

if transposed
    deltaSA = deltaSA.';
end

end
