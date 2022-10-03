function Nsquared_lowerlimit = gsw_Nsquared_lowerlimit(p,long,lat)

% gsw_Nsquared_lowerlimit             specified profile of minimum buoyancy
%                                                         frequency squared
%==========================================================================
%
% USAGE:
%  Nsquared_lowerlimit = gsw_Nsquared_lowerlimit(p,long,lat)
%
% DESCRIPTION:
%  Calculates the minimum Nsquared such that a cast is stable.
%
% INPUT:
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  long =  longitude in decimal degrees                      [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat  =  latitude in decimal degrees north                [ -90 ... +90 ] 
%
%  lat & long may have dimensions 1x1 or Mx1 or 1xN or MxN, where p is MxN.
%
% OUTPUT:
%  Nsquared_lowerlimit = Minimum Brunt-Vaisala                [ rad^2/s^2 ]
%                         Frequency squared  
%
%  The units of N2 are radians^2 s^-2 however in may textbooks this is
%  abreviated to s^-2 as radians does not have a unit.  To convert the
%  frequency to hertz, cycles sec^-1, divide the frequency by pi, 
%  i.e. N/(2 x pi).
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (4th August, 2021)
%
% REFERENCES:
%  Barker, P.M., and T.J. McDougall, 2017: Stabilizing hydrographic 
%   profiles with minimal change to the water masses. J. Atmosph. Ocean. 
%   Tech., pp. 1935 - 1945, http://dx.doi.org/10.1175/JTECH-D-16-0111.1
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.18.3) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_Nsquared_lowerlimit:  Requires three inputs')
end 

[mp,np] = size(p);
[mla,nla] = size(lat);

if (mla == 1) & (nla == 1)  
    lat = lat*ones(mp,np);
elseif (np == nla) & (mla == 1)
    lat = lat(ones(1,mp), :);    
elseif (mp == mla) & (nla == 1)  
    lat = lat(:,ones(1,np));       
elseif (np == mla) & (nla == 1)     
    lat = lat.';                  
    lat = lat(ones(1,mp), :); 
elseif (mp == nla) & (mla == 1)
    lat = lat.';
    lat = lat(:,ones(1,np));
elseif (mp == mla) & (np == nla)
    % ok
else
    error('gsw_Nsquared_lowerlimit: Inputs array dimensions arguments do not agree')
end

[mlo,nlo] = size(long);
long(long < 0) = long(long < 0) + 360; 

if (mlo == 1) & (nlo == 1)          
    long = long*ones(mp,np);
elseif (np == nlo) & (mlo == 1)     
    long = long(ones(1,mp), :);
elseif (mp == mlo) & (nlo == 1)
    long = long(:,ones(1,np));
elseif (np == mlo) & (nlo == 1) 
    long = long.';
    long = long(ones(1,mp), :);
elseif (mp == nlo) & (mlo == 1)
    long = long.';
    long = long(:,ones(1,np));
elseif (mp == mlo) & (np == nlo)
    % ok
else
    error('gsw_Nsquared_lowerlimit: Inputs array dimensions arguments do not agree')
end

p(p < -1.5) = NaN;

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
    
Nsquared_lowerlimit = (0.25 + 0.75*(exp(-p./1000))) * 1e-7;
    
end

