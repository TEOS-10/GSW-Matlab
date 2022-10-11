function [pot_enthalpy_ice_freezing_SA, pot_enthalpy_ice_freezing_P] = gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(SA,p)

% gsw_pot_enthalpy_ice_freezing_first_derivatives_poly    first derivatives
%                  of the potential enthalpy of ice at the conditions where 
%                  ice and seawater are in thermodynamic equilibrium (poly)
%==========================================================================
%
% USAGE:
%  [pot_enthalpy_ice_freezing_SA, pot_enthalpy_ice_freezing_P] = ...
%                gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(SA,p)
%
% DESCRIPTION:
%  Calculates the first derivatives of the potential enthalpy of ice Ih at
%  which ice melts into seawater with Absolute Salinity SA and at pressure 
%  p.  This code uses the computationally efficient polynomial fit of the
%  freezing potential enthalpy of ice Ih (McDougall et al., 2015).
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
%                of ice at freezing (ITS-90) with respect to Absolute
%                salinity at fixed pressure  [ (J/kg)/(g/kg) ] i.e. [ J/g ]
%                                            
%  pot_enthalpy_ice_freezing_P  = the derivative of the potential enthalpy
%                of ice at freezing (ITS-90) with respect to pressure 
%                (in Pa) at fixed Absolute Salinity           [ (J/kg)/Pa ]
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
%  McDougall et al. 2015: A reference for this polynomial.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2) 
   error('gsw_pot_enthalpy_ice_freezing_first_derivatives_poly:  Requires two inputs')
end

[ms,ns] = size(SA);
[mp,np] = size(p);

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_pot_enthalpy_ice_freezing_first_derivatives_poly: Inputs array dimensions arguments do not agree')
end 

if ms == 1
    SA = SA.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

d1 = -1.249490228128056e4;
d2 =  1.336783910789822e4;
d3 = -4.811989517774642e4;
d4 =  8.044864276240987e4;
d5 = -7.124452125071862e4;
d6 =  2.280706828014839e4;
d7 =  0.315423710959628e3;
d8 = -3.592775732074710e2;
d9 =  1.644828513129230e3;
d10 = -4.809640968940840e3;
d11 =  2.901071777977272e3;
d12 = -9.218459682855746e2;
d13 =  0.379377450285737e3;
d14 = -2.672164989849465e3;
d15 =  5.044317489422632e3;
d16 = -2.631711865886377e3;
d17 = -0.160245473297112e3;
d18 =  4.029061696035465e2;
d19 = -3.682950019675760e2;

f1 = -2.034535061416256e4;
f2 =  0.315423710959628e3;
f3 = -0.239518382138314e3;
f4 =  0.822414256564615e3;
f5 = -1.923856387576336e3;
f6 =  0.967023925992424e3;
f7 = -0.263384562367307e3;
f8 = -5.051613740291480e3;
f9 =  7.587549005714740e2;
f10 = -3.562886653132620e3;
f11 =  5.044317489422632e3;
f12 = -2.105369492709102e3;
f13 =  6.387082316647800e2;
f14 = -4.807364198913360e2;
f15 =  8.058123392070929e2;
f16 = -5.524425029513641e2;

SA_r = SA.*1e-2;
x = sqrt(SA_r);
p_r = p.*1e-4;

pot_enthalpy_ice_freezing_SA = (d1 + x.*(d2  + x.*(d3  + x.*(d4  + x.*(d5  + d6*x)))) ...
       + p_r.*(d7 + x.*(d8 + x.*(d9 + x.*(d10 + x.*(d11 + d12*x)))) ... 
       + p_r.*(d13 + x.*(d14 + x.*(d15 + d16*x)) ...
       + p_r.*(d17 + x.*(d18 + d19*x))))).*1e-2;

pot_enthalpy_ice_freezing_P = (f1 + SA_r.*(f2 + x.*(f3 + x.*(f4 + x.*(f5 + x.*(f6 + f7*x))))) ...
     + p_r.*(f8 + SA_r.*(f9 + x.*(f10 + x.*(f11 + f12*x))) ...
     + p_r.*(f13 + SA_r.*(f14 + x.*(f15 + f16*x))))).*1e-8;

% set any values that are out of range to be NaN. 
pot_enthalpy_ice_freezing_SA(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

pot_enthalpy_ice_freezing_P(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

%  The maximum error of the derivatives pot_enthalpy_ice_freezing_SA over
%  the full range of input values is -0.4 to 1 and for
%  pot_enthalpy_ice_freezing_P it is -2.5 to 2 x 10^-7.

if transposed
    pot_enthalpy_ice_freezing_SA = pot_enthalpy_ice_freezing_SA.';
    pot_enthalpy_ice_freezing_P = pot_enthalpy_ice_freezing_P.';
end

end