function pot_enthalpy_ice = gsw_pot_enthalpy_from_specvol_ice_poly(specvol_ice,p)

% gsw_pot_enthalpy_from_specvol_ice_poly            potential enthalpy from
%                                             specific volume of ice (poly)                             
%==========================================================================
%
% USAGE:
%  pot_enthalpy_ice = gsw_pot_enthalpy_from_specvol_ice_poly(specvol_ice,p)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice from the specific volume
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar. 
%
%  Note that this is a comptationally efficient polynomial version of 
%  specfic volume from potential enthalpy ice.
%
% INPUT:
%  specvol_ice =  specific volume                                [ m^3/kg ]
%  p           =  sea pressure                                     [ dbar ]
%              ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OUTPUT:
%  pot_enthalpy_ice =  potential enthalpy of ice                   [ J/kg ]
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
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_pot_enthalpy_from_specvol_ice_poly:  Requires two inputs')
end

[mh,nh] = size(specvol_ice);
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
    error('gsw_pot_enthalpy_from_specvol_ice_poly: Inputs array dimensions arguments do not agree')
end

if mh == 1
    specvol_ice = specvol_ice.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

h00 = -3.333548730778702e5;
h01 = -1.158533659785955e14;
h02 =  3.088110638980851e17;
h03 = -1.574635299523583e20;
h04 = -2.315779173662935e23;
h05 =  2.790387951606607e26;
h06 = -8.296646858006085e28;
h07 =  2.315408333461980e5;
h08 = -4.021450585146860e8;
h09 = -5.705978545141118e10;
h10 =  4.048693790458912e14;
h11 = -1.769028853661338e17;
h12 = -2.307016488735996e-2;
h13 =  2.135637957175861e1;
h14 =  8.593828063620491e-9;

pot_enthalpy_ice = h00 + specvol_ice.*(h01 + specvol_ice.*(h02 + specvol_ice.*(h03 ...
    + specvol_ice.*(h04 + specvol_ice.*(h05 + h06.*specvol_ice))))) + p.*(h07 ...
    + specvol_ice.*(h08 + specvol_ice.*(h09 + specvol_ice.*(h10 + h11.*specvol_ice))) ...
    + p.*(h12 + h13.*specvol_ice + h14.*p));

v_h_Ih = gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(pot_enthalpy_ice,p);
% initial estimate of v_h_Ih, the pot_enthalpy_ice derivative of specvol_ice

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure 
%--------------------------------------------------------------------------
 
pot_enthalpy_ice_old = pot_enthalpy_ice;
delta_v_ice = gsw_specvol_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old,p) - specvol_ice;
pot_enthalpy_ice = pot_enthalpy_ice_old - delta_v_ice./v_h_Ih ; % this is half way through the modified N-R method (McDougall and Wotherspoon, 2012)
pot_enthalpy_ice_mean = 0.5*(pot_enthalpy_ice + pot_enthalpy_ice_old);
v_h_Ih = gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(pot_enthalpy_ice_mean,p);
pot_enthalpy_ice = pot_enthalpy_ice_old - delta_v_ice./v_h_Ih; % this is the end of one loop.

pot_enthalpy_ice_old = pot_enthalpy_ice;
delta_v_ice = gsw_specvol_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old,p) - specvol_ice;
pot_enthalpy_ice = pot_enthalpy_ice_old - delta_v_ice./v_h_Ih ; % this is half way through the modified N-R method (McDougall and Wotherspoon, 2012)

% After 1.5 loops trhough the modified Newton-Raphson procedure the error
% in pot_enthalpy_ice is 8 x 10^-9 J/kg which is machine precission for 
% this calculation.

if transposed
   pot_enthalpy_ice = pot_enthalpy_ice.';
end

end

function v_h_Ih = gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(pot_enthalpy_ice,p)

%v00 =  1.090843882189585e-3;
v01 = -2.993503247981132e-10;
v02 = -2.187697844100046e-15;
v03 = -6.597339965467078e-21;
v04 = -1.092035946035247e-26;
v05 = -9.164174131511919e-33;
v06 = -2.874852704079762e-39;
%v07 = -1.558856898604863e-09;
v08 = -1.323880019968672e-15;
v09 = -2.685780994854041e-21;
v10 = -6.212785565052150e-27;
v11 = -4.797551830947803e-33;
%v12 =  2.953039321948178e-15;
v13 = -2.403547886994993e-21;
%v14 = -1.461069600352674e-20;

v_h_Ih = v01 + pot_enthalpy_ice.*(2*v02 + pot_enthalpy_ice.*(3*v03 ...
    + pot_enthalpy_ice.*(4*v04 + pot_enthalpy_ice.*(5*v05 ...
    + 6*v06.*pot_enthalpy_ice)))) + p.*(v08 + pot_enthalpy_ice.*(2*v09 ...
    + pot_enthalpy_ice.*(3*v10 + 4*v11.*pot_enthalpy_ice)) + v13.*p);

end
