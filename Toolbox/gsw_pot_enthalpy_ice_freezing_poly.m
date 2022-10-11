function pot_enthalpy_ice_freezing = gsw_pot_enthalpy_ice_freezing_poly(SA,p)

% gsw_pot_enthalpy_ice_freezing_poly              potential enthalpy of ice
%                                          at which seawater freezes (poly)
%==========================================================================
%
% USAGE:
%  pot_enthalpy_ice_freezing = gsw_pot_enthalpy_ice_freezing_poly(SA,p)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice at which seawater freezes.
%  The error of this fit ranges between -2.5 and 1 J/kg with an rms of 
%  1.07, between SA of 0 and 120 g/kg and p between 0 and 10,000 dbar (the
%  error in the fit is between -0.7 and 0.7 with an rms of
%  0.3, between SA of 0 and 120 g/kg and p between 0 and 5,000 dbar) when
%  compared with the potential enthalpy calculated from the exact in-situ 
%  freezing temperature which is found by a Newton-Raphson iteration of the 
%  equality of the chemical potentials of water in seawater and in ice.  
%  Note that the potential enthalpy at freezing can be found
%  by this exact method using the function gsw_pot_enthalpy_ice_freezing.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  p may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  pot_enthalpy_ice_freezing = potential enthalpy of ice at freezing 
%                              of seawater                         [ J/kg ]
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

if ~(nargin == 2 ) 
   error('gsw_pot_enthalpy_ice_freezing_poly:  Requires two inputs')
end

[ms,ns] = size(SA);
[mp,np] = size(p);

if (mp == 1) & (np == 1) 
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)
    p = p(ones(1,ms), :);
elseif (ms == mp) & (np == 1)
    p = p(:,ones(1,ns));
elseif (ns == mp) & (np == 1)
    p = p.';
    p = p(ones(1,ms), :);
elseif (ms == np) & (ns == mp)
    p = p.';
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_pot_enthalpy_ice_freezing_poly: Inputs array dimensions arguments do not agree')
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
SA(SA<0) = 0; 

c0 = -3.333548730778702e5;
%
c1 =  -1.249490228128056e4;
c2 =   0.891189273859881e4;
c3 =  -2.405994758887321e4;
c4 =   3.217945710496395e4;
c5 =  -2.374817375023954e4;
c6 =   0.651630522289954e4;
%
c7 =  -2.034535061416256e4;
c8 =  -0.252580687014574e4;
c9 =   0.021290274388826e4;
%  
c10 =   0.315423710959628e3;
c11 =  -0.239518382138314e3;
c12 =   0.379377450285737e3;
c13 =   0.822414256564615e3;
c14 =  -1.781443326566310e3;
c15 =  -0.160245473297112e3;
c16 =  -1.923856387576336e3;
c17 =   2.522158744711316e3;
c18 =   0.268604113069031e3;
c19 =   0.967023925992424e3;
c20 =  -1.052684746354551e3;
c21 =  -0.184147500983788e3;
c22 =  -0.263384562367307e3;

SA_r = SA.*1e-2;
x = sqrt(SA_r);
p_r = p.*1e-4;

pot_enthalpy_ice_freezing = c0 + SA_r.*(c1 + x.*(c2 + x.*(c3 + x.*(c4 + x.*(c5 + c6.*x))))) ...
            + p_r.*(c7 + p_r.*(c8 + c9.*p_r)) + SA_r.*p_r.*(c10 + p_r.*(c12 ...
            + p_r.*(c15 + c21.*SA_r)) + SA_r.*(c13 + c17.*p_r + c19.*SA_r) ...
            + x.*(c11 + p_r.*(c14 + c18.*p_r) + SA_r.*(c16 + c20.*p_r + c22.*SA_r)));

% set any values that are out of range to be NaN. 
pot_enthalpy_ice_freezing(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if transposed
    pot_enthalpy_ice_freezing = pot_enthalpy_ice_freezing.';
end

end