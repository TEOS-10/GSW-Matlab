function grav = gsw_grav(lat,p)

% gsw_grav                                       gravitational acceleration
%==========================================================================
%
% USAGE:  
%  grav = gsw_grav(lat,{p})
%
% DESCRIPTION:
%  Calculates acceleration due to gravity as a function of latitude and as
%  a function of pressure in the ocean.
%
% INPUT:
%  lat  =  latitude in decimal degress north                [ -90 ... +90 ] 
% 
% Optional:
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar ) 
%  (If pressure is not given then it is assumed that pressure = 0 dbar.)
%
%  p (if provided) may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where lat is MxN.
%
% OUTPUT:
%  grav  =  gravitational acceleration                           [ m s^-2 ]
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
%    See appendix D of this TEOS-10 Manual. 
%
%  Moritz, H., 2000: Geodetic reference system 1980. J. Geodesy, 74, 
%   pp. 128-133.
%
%  Saunders, P.M., and N.P. Fofonoff, 1976: Conversion of pressure to
%   depth in the ocean. Deep-Sea Res., pp. 109-111.
%   
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 1 | nargin == 2)
   error('gsw_grav:  Requires either one or two inputs, latitude and pressure')
end 

if nargin == 1
%    Assume no change in gravity with respect to pressure.
  p = zeros(size(lat));
end 

[ml,nl] = size(lat);
[mp,np] = size(p);

if (ml == 1) & (nl == 1)                        % lat is a scalar.  Fill to size of p
    lat = lat*ones(mp,np);
elseif (nl == np) & (ml == 1)                   % lat is row vector, 
    lat = lat(ones(1,mp),:);                %   copy down each column.
elseif (ml == mp) & (nl == 1)                   % lat is column vector,
    lat = lat(:,ones(1,np));                %   copy across each row.
elseif (nl == mp) & (nl == 1)          % p is a transposed row vector,
    lat = lat.';                              % transposed then
    lat = lat(ones(1,mp), :);                % copy down each column.
elseif (ml == mp) & (nl == np)               
    % ok
end

[ml,nl] = size(lat);

if (mp == 1) & (np == 1)               % p is a scalar.  Fill to size of lat
    p = p(1)*ones(ml,nl);
elseif (np == nl) & (mp == 1)                   % p is row vector,
    p = p(ones(1,ml),:);                %   copy down each column.
elseif (mp == ml) & (np == 1)                   % p is column vector,
    p = p(:,ones(1,nl));                %   copy across each row.
elseif (np == ml) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ml), :);                % copy down each column.
elseif (mp == ml) & (np == nl)
    % ok
else
    error('gsw_grav: p has wrong dimensions')
end

if ml == 1
    lat = lat.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

z = gsw_z_from_p(p,lat);
clear p

gamma = 2.26e-7;
deg2rad = pi/180;
sinlat = sin(lat*deg2rad);  % convert to radians
clear lat deg2rad

sin2 = sinlat.*sinlat;
clear sinlat
gs = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2)).*sin2) ;

grav = gs.*(1 - gamma*z);             % z is the height corresponding to p. 
                                        % Note. In the ocean z is negative.
if transposed
    grav = grav.';
end

end
