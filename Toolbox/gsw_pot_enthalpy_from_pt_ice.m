function pot_enthalpy_ice = gsw_pot_enthalpy_from_pt_ice(pt0_ice)

% gsw_pot_enthalpy_from_pt_ice               potential enthalpy of ice from 
%                              potential temperature refered to the surface
%==========================================================================
%
% USAGE:
%  pot_enthalpy_ice = gsw_pot_enthalpy_from_pt_ice(pt0_ice)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice from potential temperature of
%  ice (whose reference sea pressure is zero dbar).  
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
   error('gsw_pot_enthalpy_from_pt_ice:  Requires only one input')
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

Tt = 273.16;                        % Triple-point temperature, kelvin (K).  
rec_Tt = 3.660858105139845e-3;               % 1/Tt = 3.660858105139845e-3; 

T = pt0_ice + gsw_T0;    % The input temperature t is potential temperature
                      % refered to 0 dbar (the surface) in units of degrees 
                   % Celcius.  T is the in-situ Absolute Temperature of the 
                                               % ice in degrees kelvin (K).  
tau = T.*rec_Tt;

g00 = -6.32020233335886e5;

t1 = (3.68017112855051e-2 + 5.10878114959572e-2i);
t2 = (3.37315741065416e-1 + 3.35449415919309e-1i);

r1 = (4.47050716285388e1 + 6.56876847463481e1i);
r20	= (-7.25974574329220e1 - 7.81008427112870e1i);

sqtau_t1 = (tau./t1).^2;
sqtau_t2 = (tau./t2).^2;

h0_part = r1.*t1.*(log(1 - sqtau_t1) + sqtau_t1) ...
          + r20.*t2.*(log(1 - sqtau_t2) + sqtau_t2);

pot_enthalpy_ice = g00 + Tt.*real(h0_part); 
              
if transposed
   pot_enthalpy_ice = pot_enthalpy_ice.';
end

end