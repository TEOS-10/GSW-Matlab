function specvol_ice = gsw_specvol_from_pot_enthalpy_ice_poly(pot_enthalpy_ice,p)

% gsw_specvol_from_pot_enthalpy_ice_poly               specific volume from
%                                          potential enthalpy of ice (poly)                              
%==========================================================================
%
% USAGE:
%  specvol_ice = gsw_specvol_from_pot_enthalpy_ice_poly(pot_enthalpy_ice,p)
%
% DESCRIPTION:
%  Calculates the specific volume of ice from the potential enthalpy 
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar. 
%
%  Note that this is a comptationally efficient polynomial fit of the
%  programme gsw_specvol_from_pot_enthalpy_ice(pot_enthalpy_ice,p).
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OUTPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_specvol_from_pot_enthalpy_ice_poly:  Requires two inputs')
end

[mh,nh] = size(pot_enthalpy_ice);
[mp,np] = size(p);

if (mp == 1) & (np == 1)           
    p = p*ones(mh,nh);
elseif (nh == np) & (mp == 1)         
    p = p(ones(1,mh), :);             
elseif (mh == mp) & (np == 1)        
    p = p(:,ones(1,nh));               
elseif (nh == mp) & (np == 1)        
    p = p.';                     
    p = p(ones(1,mh), :);          
elseif (mh == mp) & (nh == np)
    % ok
else
    error('gsw_specvol_from_pot_enthalpy_ice_poly: Inputs array dimensions arguments do not agree')
end

if mh == 1
    pot_enthalpy_ice = pot_enthalpy_ice.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

v00 =  1.090843882189585e-3;
v01 = -2.993503247981132e-10;
v02 = -2.187697844100046e-15;
v03 = -6.597339965467078e-21;
v04 = -1.092035946035247e-26;
v05 = -9.164174131511919e-33;
v06 = -2.874852704079762e-39;
v07 = -1.558856898604863e-09;
v08 = -1.323880019968672e-15;
v09 = -2.685780994854041e-21;
v10 = -6.212785565052150e-27;
v11 = -4.797551830947803e-33;
v12 =  2.953039321948178e-15;
v13 = -2.403547886994993e-21;
v14 = -1.461069600352674e-20;

specvol_ice = v00 + pot_enthalpy_ice.*(v01 + pot_enthalpy_ice.*(v02 + pot_enthalpy_ice.*(v03 ...
    + pot_enthalpy_ice.*(v04 + pot_enthalpy_ice.*(v05 + v06.*pot_enthalpy_ice))))) + p.*(v07 ...
    + pot_enthalpy_ice.*(v08 + pot_enthalpy_ice.*(v09 + pot_enthalpy_ice.*(v10 + v11.*pot_enthalpy_ice))) ...
    + p.*(v12 + v13.*pot_enthalpy_ice + v14.*p));

if transposed
   specvol_ice = specvol_ice.';
end

end
