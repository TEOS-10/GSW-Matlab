function enthalpy_ice = gsw_enthalpy_ice(t,p)

% gsw_enthalpy_ice                                 specific enthalpy of ice
%==========================================================================
%
% USAGE:
%  enthalpy_ice = gsw_enthalpy_ice(t,p)
%
% DESCRIPTION:
%  Calculates the specific enthalpy of ice (h_Ih). 
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where t is MxN.
%
% OUTPUT:
%  enthalpy_ice  =  specific enthalpy of ice                       [ J/kg ]
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.06.12 (25th May, 2020)
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

if ~(nargin==2)
   error('gsw_enthalpy_ice:  Requires two inputs')
end %if

[mt,nt] = size(t);
[mp,np] = size(p);

if (mp == 1) & (np == 1)              % p scalar - fill to size of t
    p = p*ones(size(t));
elseif (nt == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,mt), :);              % copy down each column.
elseif (mt == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,nt));               % copy across each row.
elseif (nt == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,mt), :);                % copy down each column.
elseif (mt == mp) & (nt == np)
    % ok
else
    error('gsw_enthalpy_ice: Inputs array dimensions arguments do not agree')
end 

if mt == 1
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

rec_Pt = 1.634903221903779e-3;   % 1./Pt, where Pt = 611.657;  Experimental 
                                 % triple-point pressure in Pa.

Tt = 273.16;  % Triple-point temperature, kelvin (K).  
rec_Tt = 3.660858105139845e-3;   % 1/Tt = 3.660858105139845e-3; 

T = t + gsw_T0; % The input temperature t is in-situ temperature in
                % units of degrees Celcius.  T is the in-situ Absolute 
                % Temperature of the ice in degrees kelvin (K).  
tau = T.*rec_Tt;

db2Pa = 1e4;
dzi = db2Pa.*p.*rec_Pt;

g00 = -6.32020233335886e5;
g01 =  6.55022213658955e-1;
g02 = -1.89369929326131e-8;
g03 =  3.3974612327105304e-15;
g04 = -5.564648690589909e-22;

t1 = (3.68017112855051e-2 + 5.10878114959572e-2i);
t2 = (3.37315741065416e-1 + 3.35449415919309e-1i);

r1 = (4.47050716285388e1 + 6.56876847463481e1i);
r20	= (-7.25974574329220e1 - 7.81008427112870e1i);
r21	= (-5.57107698030123e-5 + 4.64578634580806e-5i);
r22	= (	2.34801409215913e-11 - 2.85651142904972e-11i);

g0 = g00 + dzi.*(g01 + dzi.*(g02 + dzi.*(g03 + g04.*dzi)));

r2 = r20 + dzi.*(r21 + r22.*dzi);

sqtau_t1 = (tau./t1).^2;
sqtau_t2 = (tau./t2).^2;

g = r1.*t1.*(log(1 - sqtau_t1) + sqtau_t1) ...
    + r2.*t2.*(log(1 - sqtau_t2) + sqtau_t2);

enthalpy_ice = g0 + Tt.*real(g); 
           
if transposed
    enthalpy_ice = enthalpy_ice.';
end

end