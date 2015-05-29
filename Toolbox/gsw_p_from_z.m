function p = gsw_p_from_z(z,lat)

% gsw_p_from_z                                         pressure from height
%==========================================================================
%
% USAGE:
%   p = gsw_p_from_z(z,lat)
%
% DESCRIPTION:
%  Calculates sea pressure from height using computationally-efficient 
%  25-term expression for density, in terms of SA, CT and p 
%  (McDougall et al., 2010).
%  Note. Height (z) is NEGATIVE in the ocean.  Depth is -z.  
%    Depth is not used in the GSW computer software library. 
%
% INPUT:
%  z   =  height                                                      [ m ]
%   Note. At sea level z = 0, and since z (HEIGHT) is defined 
%     to be positive upwards, it follows that while z is 
%     positive in the atmosphere, it is NEGATIVE in the ocean.
%  lat =  latitude in decimal degrees north                 [ -90 ... +90 ]
%   
%  lat may have dimensions 1x1 or Mx1 or 1xN or MxN, where z is MxN.
%
% OUTPUT:
%   p  =  sea pressure                                             [ dbar ]
%         ( ie. absolute pressure - 10.1325 dbar )
%
% AUTHOR: 
%  Trevor McDougall, Claire Roberts-Thomson and Paul Barker. 
%                                                  [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
% McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%  Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
%  expression for the density of seawater in terms of Conservative 
%  Temperature, and related properties of seawater.  To be submitted 
%  to Ocean Science Discussions. 
%
% IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%  seawater - 2010: Calculation and use of thermodynamic properties.
%  Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%  UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
% Saunders, P. M., 1981: Practical conversion of pressure to depth. 
%  Journal of Physical Oceanography, 11, 573-574.
%
% Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.
%
% This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
    error('gsw_p_from_z: Requires two inputs, latitude and pressure')
end %if

[mz,nz] = size(z);
[ml,nl] = size(lat);

if (ml == 1) & (nl == 1)              % lat is a scalar - fill to size of z
    lat = lat*ones(size(z));
elseif (nl == nz) & (ml == 1)         % lat is row vector,
    lat = lat(ones(1,mz), :);              % copy down each column.
elseif (mz == ml) & (nl == 1)         % lat is column vector,
    lat = lat(:,ones(1,nz));               % copy across each row.
elseif (mz == ml) & (nz == nl)
    % ok
else
    error('gsw_p_from_z: Inputs array dimensions arguments do not agree')
end %if

if mz == 1
    z = z';
    lat = lat';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the caclulation
%--------------------------------------------------------------------------

db2Pa = 1e4; 
gamma = 2.26e-7;
DEG2RAD = pi/180;
X     = sin(lat*DEG2RAD);
sin2  = X.*X;
gs    = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2)).*sin2);

% get the first estimate of p from Saunders (1981)
c1 =  5.25e-3*sin2 + 5.92e-3;
p  = -2.*z./((1-c1) + sqrt((1-c1).*(1-c1) + 8.84e-6.*z));
% end of the first estimate from Saunders (1981)

df_dp = db2Pa * gsw_specvol_SSO_0_CT25(p); % initial value of the derivative of f

f     = gsw_enthalpy_SSO_0_CT25(p) + gs.*(z - 0.5*gamma*(z.*z));
p_old = p;
p     = p_old - f./df_dp;
pm    = 0.5*(p + p_old);
df_dp = db2Pa * gsw_specvol_SSO_0_CT25(pm);
p     = p_old - f./df_dp;

% After this one iteration through this modified Newton-Raphson iterative
% procedure, the remaining error in p is at computer machine precision,
% being no more than 1.6e-10 dbar. 

if transposed
    p = p';
end

end
