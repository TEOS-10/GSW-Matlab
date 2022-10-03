function pt0_ice = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice)

% gsw_pt_from_pot_enthalpy_ice_poly        potential temperature refered to
%                         the surface from potential enthalpy of ice (poly)                             
%==========================================================================
%
% USAGE:
%  pt0_ice = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice)
%
% DESCRIPTION:
%  Calculates the potential temperature of ice (whose reference sea 
%  pressure is zero dbar) from the potential enthalpy of ice.  This is a
%  compuationally efficient polynomial fit to the potential enthalpy of
%  ice.
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%
% OUTPUT:
%  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
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
   error('gsw_pt_from_pot_enthalpy_ice_poly:  Requires only one input')
end %if

[mh,nh] = size(pot_enthalpy_ice);

if mh == 1
    pot_enthalpy_ice = pot_enthalpy_ice.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

q0 = 2.533588268773218e2;
q1 = 2.594351081876611e-3;
q2 = 1.765077810213815e-8;
q3 = 7.768070564290540e-14;
q4 = 2.034842254277530e-19;
q5 = 3.220014531712841e-25;
q6 = 2.845172809636068e-31;
q7 = 1.094005878892950e-37;
    
pt0_ice = q0 + pot_enthalpy_ice.*(q1 + pot_enthalpy_ice.*(q2 + pot_enthalpy_ice.*(q3 ...
          + pot_enthalpy_ice.*(q4 + pot_enthalpy_ice.*(q5 + pot_enthalpy_ice.*(q6 ...
          + pot_enthalpy_ice.*q7))))));

% The error of this fit ranges between -5e-5 and 2e-4 deg C over the potential 
% temperature range of -100 to 2 deg C, or the potential enthalpy range of 
% -5.7 x 10^5 to -3.3 x 10^5 J/kg. 
      
if transposed
   pt0_ice = pt0_ice.';
end

end