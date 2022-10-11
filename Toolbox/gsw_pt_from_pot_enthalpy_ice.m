function pt0_ice = gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_ice)

% gsw_pt_from_pot_enthalpy_ice      potential temperature of ice refered to
%                                the surface from potential enthalpy of ice                              
%==========================================================================
%
% USAGE:
%  pt0_ice = gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_ice)
%
% DESCRIPTION:
%  Calculates the potential temperature of ice from the potential enthalpy 
%  of ice.  The reference sea pressure of both the potential temperature 
%  and the potential enthalpy is zero dbar. 
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_pt_from_pot_enthalpy_ice:  Requires only one input')
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

h00 = -6.320202333358860e5;           % this is gsw_enthalpy_ice(-gsw_T0,0)

pot_enthalpy_ice(pot_enthalpy_ice < h00) = h00;
   
p0 = zeros(size(pot_enthalpy_ice));

% The below polynomial is the starting polynomial for pt0 of ice.  
pt0_ice = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice);

% The variable "df_dt" below is the derivative of the above polynomial with
% respect to pot_enthalpy_ice.  This is the initial value of the derivative
% of the function f.
recip_df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice);

pt0_ice_old = pt0_ice;
f = gsw_pot_enthalpy_from_pt_ice(pt0_ice_old) - pot_enthalpy_ice;
pt0_ice = pt0_ice_old - f.*recip_df_dt ; % this is half way through the modified Newton's method (McDougall and Wotherspoon, 2012)
ptm_ice = 0.5.*(pt0_ice + pt0_ice_old);
recip_df_dt = 1./gsw_cp_ice(ptm_ice,p0);
pt0_ice = pt0_ice_old - f.*recip_df_dt;
%  This is the end of one full iteration of the modified Newton's Method. 
%  For input potential enthalpies greater than -5.1e-5, the above part of
%  the code gives the output potential temperature of ice accurate to 1e-13
%  degrees C.  The code below is only entered if pot_enthalpy_ice < -5.1e5,
%  that is, if the potential temperature of ice is less than about 
%  -100 deg C.  

if any(pot_enthalpy_ice < -5.1e5)  % equivalently, for pt0_ice less than about -100 deg C
    % These temperatures are less than those found in nature on planet Earth. 
    [I] = find(pot_enthalpy_ice < -5.1e5);
    
    pt0_cold_ice = gsw_pt0_cold_ice_poly(pot_enthalpy_ice(I)); 
    % the starting polynomial for pt0 of ice that has a potential enthalpy < -5.1e5
    
    df_dt = gsw_cp_ice(pt0_cold_ice+0.02,p0(I));% The heat capacity, cp, is
    % evaluated at 0.02 C greater than usual in order to avoid stability
    % issues and to ensure convergence near zero Absolute Temperature. 
    
    for Iteration = 1:6  
        pt0_cold_ice_old = pt0_cold_ice;
        f = gsw_pot_enthalpy_from_pt_ice(pt0_cold_ice_old) - pot_enthalpy_ice(I);
        pt0_cold_ice = pt0_cold_ice_old - f./df_dt ; % this is half way through the modified Newton's method
        ptm_cold_ice = 0.5.*(pt0_cold_ice + pt0_cold_ice_old);        
        df_dt = gsw_cp_ice(ptm_cold_ice+0.02,p0(I));% note the extra 0.02 C here as well    
        pt0_cold_ice = pt0_cold_ice_old - f./df_dt;
    end
    pt0_ice(I) = pt0_cold_ice;
end

%The potential temerature has a maximum error as listed in the table below.
%
%  potential temerature error (deg C)  |  @ potential temerature (deg C)
%--------------------------------------|--------------------------------
%                0.012                 |     -273.15 to -273.12
%              4 x 10^-4               |     -232.12 to -273.0
%             2.5 x 10^-6              |          -273
%              7 x 10^-9               |          -272
%            3.7 x 10^-10              |          -270
%              6 x 10^-11              |          -268
%             2.5 x 10^11              |          -266
%             3 x 10^-12               |          -260
%             7 x 10^-13               |          -250
%            2.2 x 10^-13              |          -220
%            1.7 x 10^-13              |         >= -160
%  
% Note.  The above errors in each temperature range are machine precissions 
% for this calculation.

if transposed
   pt0_ice = pt0_ice.';
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

q1 = 2.594351081876611e-3;
r2 = 3.530155620427630e-8;
r3 = 2.330421169287162e-13;
r4 = 8.139369017110120e-19;
r5 = 1.610007265856420e-24;
r6 = 1.707103685781641e-30;
r7 = 7.658041152250651e-37;

dpt0_ice_dh = q1 + pot_enthalpy_ice.*(r2 + pot_enthalpy_ice.*(r3 ...
              + pot_enthalpy_ice.*(r4 + pot_enthalpy_ice.*(r5  ...
              + pot_enthalpy_ice.*(r6 + pot_enthalpy_ice.*r7)))));


end

%##########################################################################

function pt0_cold_ice_poly = gsw_pt0_cold_ice_poly(pot_enthalpy_ice)

% gsw_pt0_cold_ice_poly                  a polynomial for pt0 for potential 
%                               enthalpy < -5.1e5  (pt0_ice < about -100 C)                             
%==========================================================================
%
% USAGE:
%  pt0_cold_poly = gsw_pt0_cold_ice_poly(pot_enthalpy_ice)
%
% DESCRIPTION:
%  Calculates an initial estimate of pt0_ice when it is less than about
%  -100 deg C. 
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%
% OUTPUT:
%  pt0_cold_ice_poly  =  initial estimate of potential temperatur 
%                        of very cold ice in dgress C (not K)     [ deg C ]                                        
%
% AUTHOR: 
%  Paul Barker                                         [ help@teos-10.org ]
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

h00 = -6.320202333358860e5;           % this is gsw_enthalpy_ice(-gsw_T0,0)

s0 =  1.493103204647916;
s1 =  2.372788609320607e-1;
s2 = -2.014996002119374e-3;
s3 =  2.640600197732682e-6;
s4 =  3.134706016844293e-5;
s5 =  2.733592344937913e-6;
s6 =  4.726828010223258e-8;
s7 = -2.735193883189589e-9;
s8 = -8.547714991377670e-11;

log_h_diff = log(pot_enthalpy_ice - h00);

log_Abs_theta0 = s0 + log_h_diff.*(s1 + log_h_diff.*(s2 + log_h_diff.*(s3 ...
                + log_h_diff.*(s4 + log_h_diff.*(s5 + log_h_diff.*(s6 ...
                + log_h_diff.*(s7 + log_h_diff.*s8)))))));

pt0_cold_ice_poly = exp(log_Abs_theta0) - 273.15;
    
end