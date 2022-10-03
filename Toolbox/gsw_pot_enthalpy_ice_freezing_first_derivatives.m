function [pot_enthalpy_ice_freezing_SA, pot_enthalpy_ice_freezing_P] = gsw_pot_enthalpy_ice_freezing_first_derivatives(SA,p)

% gsw_pot_enthalpy_ice_freezing_first_derivatives         first derivatives
%                        of the potential enthalpy of ice at the conditions  
%                   where ice and seawater are in thermodynamic equilibrium 
%==========================================================================
%
% USAGE:
%  [pot_enthalpy_ice_freezing_SA, pot_enthalpy_ice_freezing_P] = ...
%                     gsw_pot_enthalpy_ice_freezing_first_derivatives(SA,p)
%
% DESCRIPTION:
%  Calculates the first derivatives of the potential enthalpy of ice at
%  which seawater freezes, with respect to Absolute Salinity SA and
%  pressure P (in Pa).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  pot_enthalpy_ice_freezing_SA = the derivative of the potential enthalpy
%                  of ice at freezing (ITS-90) with respect to Absolute
%                  salinity at fixed pressure  [ K/(g/kg) ] i.e. [ K kg/g ]
%
%  pot_enthalpy_ice_freezing_P  = the derivative of the potential enthalpy
%                  of ice at freezing (ITS-90) with respect to pressure 
%                  (in Pa) at fixed Absolute Salinity              [ K/Pa ]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.  
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
   error('gsw_pot_enthalpy_ice_freezing_first_derivatives:  Requires two inputs')
end

[ms,ns] = size(SA);
[mp,np] = size(p);
saturation_fraction = zeros(ms,ns);

if (mp == 1) & (np == 1)
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)
    p = p(ones(1,ms), :);
elseif (ms == mp) & (np == 1)
    p = p(:,ones(1,ns));
elseif (ns == mp) & (np == 1)
    p = p.';
    p = p(ones(1,ms), :);
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_pot_enthalpy_ice_freezing_first_derivatives: Inputs array dimensions arguments do not agree')
end

if ms == 1
    SA = SA.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;
        
tf = gsw_t_freezing(SA,p,saturation_fraction);
pt_icef = gsw_pt0_from_t_ice(tf,p); 
ratio_temp = (gsw_T0 + pt_icef)./(gsw_T0 + tf);

cp_Ihf = gsw_cp_ice(tf,p);

[tf_SA, tf_P] = gsw_t_freezing_first_derivatives(SA,p);

pot_enthalpy_ice_freezing_SA = ratio_temp.*cp_Ihf.*tf_SA;

pot_enthalpy_ice_freezing_P = ratio_temp.*cp_Ihf.*tf_P ...
                             - (gsw_T0 + pt_icef).*gsw_gibbs_ice(1,1,tf,p);

% set any values that are out of range to be NaN. 
pot_enthalpy_ice_freezing_SA(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

pot_enthalpy_ice_freezing_P(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if transposed
    pot_enthalpy_ice_freezing_SA = pot_enthalpy_ice_freezing_SA.';
    pot_enthalpy_ice_freezing_P = pot_enthalpy_ice_freezing_P.';
end

end