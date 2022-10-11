function pot_enthalpy_ice = gsw_pot_enthalpy_from_pt_ice_poly(pt0_ice)

% gsw_pot_enthalpy_from_pt_ice_poly          potential enthalpy of ice from 
%                              potential temperature refered to the surface
%==========================================================================
%
% USAGE:
%  pot_enthalpy_ice = gsw_pot_enthalpy_from_pt_ice_poly(pt0_ice)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice from potential temperature of
%  ice (whose reference sea pressure is zero dbar).  This is a
%  compuationally efficient polynomial fit to the potential enthalpy of
%  ice.
%   
% INPUT:
%  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
%
% OUTPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
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
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_pot_enthalpy_from_pt_ice_poly:  Requires only one input')
end %if

[mt,nt] = size(pt0_ice);

if mt == 1
    pt0_ice = pt0_ice.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% initial estimate of the potential enthalpy.
p0 = -3.333601570157700e5;
p1 =  2.096693916810367e3;
p2 =  3.687110754043292;
p3 =  4.559401565980682e-4;
p4 = -2.516011957758120e-6;
p5 = -1.040364574632784e-8;
p6 = -1.701786588412454e-10;
p7 = -7.667191301635057e-13;
    
pot_enthalpy_ice = p0 + pt0_ice.*(p1 + pt0_ice.*(p2 + pt0_ice.*(p3 ...
                   + pt0_ice.*(p4 + pt0_ice.*(p5 + pt0_ice.*(p6 ...
                   + pt0_ice.*p7))))));
               
df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice);

for Iteration = 1:5
    pot_enthalpy_ice_old = pot_enthalpy_ice;
    f = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old) - pt0_ice;
    pot_enthalpy_ice = pot_enthalpy_ice_old - f./df_dt;
    pot_enthalpy_ice_mid = 0.5*(pot_enthalpy_ice + pot_enthalpy_ice_old);
    df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice_mid);
    pot_enthalpy_ice = pot_enthalpy_ice_old - f./df_dt;
end
               
% The error of this fit ranges between -6e-3 and 6e-3 J/kg over the potential 
% temperature range of -100 to 2 deg C, or the potential enthalpy range of 
% -5.7 x 10^5 to -3.3 x 10^5 J/kg. 
       
if transposed
   pot_enthalpy_ice = pot_enthalpy_ice.';
end

end

%##########################################################################

function dpt0_ice_dh = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice)

% gsw_pt_from_pot_enthalpy_ice_poly_dh                    derivative of pt0
%                                    w.r.t potential enthalpy of ice (poly)                              
%==========================================================================
%
% USAGE:
%  dpt0_ice_dh = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice)
%
% DESCRIPTION:
%  Calculates the derivative of potential temperature of ice with respect 
%  to potential enthalpy.  This is based on the compuationally-efficient 
%  polynomial fit to the potential enthalpy of ice. 
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%
% OUTPUT:
%  dpt0_ice_dh  =  derivative of potential temperature of ice 
%                  with respect to potential enthalpy             [ deg C ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
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

q1 = 2.594351081876611e-3;
r2 = 3.530155620427630e-8;
r3 = 2.330421169287162e-13;
r4 = 8.139369017110120e-19;
r5 = 1.610007265856420e-24;
r6 = 1.707103685781641e-30;
r7 = 7.658041152250651e-37;

dpt0_ice_dh = q1 + pot_enthalpy_ice.*(r2 + pot_enthalpy_ice.*(r3 ...
    + pot_enthalpy_ice.*(r4 + pot_enthalpy_ice.*(r5 + pot_enthalpy_ice.*(r6 ...
    + pot_enthalpy_ice.*r7)))));


end

