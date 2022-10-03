function distance = gsw_distance(long,lat,p)

% gsw_distance                      spherical earth distance between points 
%                                in long, lat coordinates at given pressure 
%==========================================================================
%
% USAGE:  
%  distance = gsw_distance(long,lat,{p})
%
% DESCRIPTION:
%  Calculates the distance in metres between successive points in the 
%  vectors long and lat, computed using the Haversine formula on a 
%  spherical earth of radius 6,371 km, being the radius of a sphere having
%  the same volume as Earth.  For a sperical Earth of radius 6,371,000 m,
%  one nautical mile is 1,853.2488 m, thus one degree of latitude is
%  111,194.93 m.
%  Note. Distances are probably good to better than 1% of the "true" 
%   distance on the ellipsoidal earth.
%
% INPUT:
%  long  =  longitude in decimal degress                     [ 0 ... +360 ] 
%                                                       or [-180 ... +180 ]   
%  lat   =  latitude in decimal degress north               [ -90 ... +90 ]  
%
% OPTIONAL:
%  p     =  sea pressure ( default is 0 )                          [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%
%  lat and long need to have the same dimensions, Mx1 or 1xN or MxN.
%  p, if provided, may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where lat & long are Mx1 or 1xN or MxN.
%
% OUTPUT:
%  distance  =  Distance between points on a spherical                [ m ]
%               Earth at pressure (p)  
%  Note. The output is in m not km.
% 
% AUTHOR:  
%  6th November, 2000 by Rich Pawlowicz                [ help@teos-10.org ]
%  Note. This function was extracted from Rich Pawlowicz's m_map package,
%    which is available from http://www.eos.ubc.ca/~rich/map.html
%
% MODIFIED:
%  4th April, 2011 by Paul Barker and Trevor McDougall. 
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCE:
%  http://www.eos.ubc.ca/~rich/map.html
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3)
   error('gsw_distance:  Requires either two or three inputs')
end %if

if nargin == 2
  p = zeros(size(lat));
end %if
    
[mla,nla] = size(lat);
[mlo,nlo] = size(long);
[mp,np] = size(p);

if (mla ~= mlo | nla ~= nlo )
    error('*** Input array dimensions in gsw_distance do not agree  ***')
elseif mla == 1 & mlo == 1 & nla == 1 & nlo == 1
    error('*** No, you need more than one point to find a distance!  ***')  
end

transposed = 0;
if (mla == 1) & (nla == 1)            % lat is a scalar.  Fill to size of p
    lat = lat*ones(mp,np);
    long = long*ones(mp,np);
elseif (nla == np) & (mla == 1)                   % lat is row vector, 
    lat = lat(ones(1,mp),:);                %   copy down each column.
    long = long(ones(1,mp),:);
elseif (mla == mp) & (nla == 1)                   % lat is column vector.
    lat = lat.';
    long = long.';
    transposed = 1;
    lat = lat(ones(1,np),:);                
    long = long(ones(1,np),:);  
elseif (mla == np) & (nla == 1)                   % lat is column vector and p is a column vector
    lat = lat.';                                  
    long = long.';
    transposed = 1;
elseif (mla == mp) & (nla == np)               
    % ok
end %if

[mla,nla] = size(lat);

if (mp == 1) & (np == 1)              % p is a scalar.  Fill to size of lat
    p = p*ones(mla,nla);
elseif (np == nla) & (mp == 1)                   % p is row vector, 
    p = p(ones(1,mla),:);                %   copy down each column.
elseif (mp == nla) & (np == mla)                   % p is row vector, 
    p = p.';    
elseif (np == nla) & (np == 1)                   % p is column vector.
    p = p.';                                     % transpose.
    p = p(:,ones(1,nla));                 %   copy along each row.
elseif (mp == mla) & (np == nla)               
    % ok
else
    error('gsw_distance: p has wrong dimensions')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

pi180 = pi/180;
earth_radius = 6371000;                         % Earth's radius in metres.
dlong = pi180*(long(:,2:nla) - long(:,1:nla-1));
dlat = pi180*(lat(:,2:nla) - lat(:,1:nla-1));

a = (sin(0.5*dlat)).^2 + cos(lat(:,1:nla-1)*pi180).*cos(lat(:,2:nla)*pi180).*(sin(0.5*dlong)).^2;
angles = 2 * atan2(sqrt(a),sqrt(1-a));

p_mid = 0.5*(p(:,1:nla-1) + p(:,2:nla));
lat_mid = 0.5*(lat(:,1:nla-1) + lat(:,2:nla));
z = gsw_z_from_p(p_mid,lat_mid);        % Note. z is height and is negative
                                                            % in the ocean.
                                                            
distance = (earth_radius + z).*angles;   % Note. The output is in m not km.

if transposed
    distance = distance.';
end

end
