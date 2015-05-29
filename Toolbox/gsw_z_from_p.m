function z = gsw_z_from_p(p,lat)

% gsw_z_from_p                                         height from pressure
%==========================================================================
%
% USAGE:  
%   z = gsw_z_from_p(p,lat)
%
% DESCRIPTION:
%  Calculates height from sea pressure using the computationally-efficient
%   25-term expression for density in terms of SA, CT and p.
%   (McDougall et al., 2010)  
%  Note. Height z is NEGATIVE in the ocean. ie. Depth is -z.  
%   Depth is not used in the GSW computer software library.  
%
% INPUT:
%   p   =  sea pressure                                            [ dbar ]
%           ( ie. absolute pressure - 10.1325 dbar )
%   lat =  latitude in decimal degrees north                [ -90 ... +90 ]
%   lat may have dimensions 1x1 or Mx1 or 1xN or MxN, where p is MxN.
%
% OUTPUT:
%   z   =  height                                                     [ m ]
%  Note. At sea level z = 0, and since z (HEIGHT) is defined to be
%    positive upwards, it follows that while z is positive in the 
%    atmosphere, it is NEGATIVE in the ocean.
%
% AUTHOR:  
% Trevor McDougall, Claire Roberts-Thomson & Paul Barker.
%                                                    [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
%   expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted 
%   to Ocean Science Discussions. 
%
%  Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.
%   
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_z_from_p: Requires two inputs, pressure and latitude')
end %if

[mp,np] = size(p);
[ml,nl] = size(lat);

if (ml == 1) & (nl == 1)              % lat scalar - fill to size of p
    lat = lat*ones(size(p));
elseif (nl == np) & (ml == 1)         % lat is row vector,
    lat = lat(ones(1,mp), :);              % copy down each column.
elseif (mp == ml) & (nl == 1)         % lat is column vector,
    lat = lat(:,ones(1,np));               % copy across each row.
elseif (mp == ml) & (np == nl)
    % ok
else
    error('gsw_z_from_p: Inputs array dimensions arguments do not agree')
end %if

if mp == 1
    p = p';
    lat = lat';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

gamma = 2.26e-7 ; 
DEG2RAD = pi/180;
X     = sin(lat*DEG2RAD);
sin2  = X.*X;
B     = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2)).*sin2); 
A     = -0.5*gamma*B ;
C     = gsw_enthalpy_SSO_0_CT25(p);
z     = -2*C./(B + sqrt(B.*B - 4.*A.*C));

if transposed
    z = z';
end

end
