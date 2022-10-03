function latentheat_evap = gsw_latentheat_evap_CT(SA,CT)

% gsw_latentheat_evap_CT                         latent heat of evaporation 
%==========================================================================
%
% USAGE: 
%  latentheat_evap = gsw_latentheat_evap_CT(SA,CT)
%
% DESCRIPTION:
%  Calculates latent heat, or enthalpy, of evaporation at p = 0 (the 
%  surface).  It is defined as a function of Absolute Salinity, SA, and
%  Conservative Temperature, CT, and is valid in the ranges 
%  0 < SA < 42 g/kg and 0 < CT < 40 deg C.  The errors range between 
%  -0.4 and 0.6 J/kg.
%  
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  latentheat_evap = latent heat of evaporation                    [ J/kg ]
%
% AUTHOR:  
%  Paul Barker, Trevor McDougall & Rainer Feistel      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See section 3.39 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_latentheat_evap_CT:  Requires two input arguments')
end %if
[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_latentheat_evap_CT: SA and CT must have same dimensions')
end

if ms == 1
    SA = SA.';
    CT = CT.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
 
c0 =   2.499065844825125e6;
c1 =  -1.544590633515099e-1;
c2 =  -9.096800915831875e4;
c3 =   1.665513670736000e2;
c4 =   4.589984751248335e1;
c5 =   1.894281502222415e1;
c6 =   1.192559661490269e3;
c7 =  -6.631757848479068e3;
c8 =  -1.104989199195898e2;
c9 =  -1.207006482532330e3;
c10 = -3.148710097513822e3;
c11 =  7.437431482069087e2;
c12 =  2.519335841663499e3;
c13 =  1.186568375570869e1;
c14 =  5.731307337366114e2;
c15 =  1.213387273240204e3;
c16 =  1.062383995581363e3;
c17 = -6.399956483223386e2;
c18 = -1.541083032068263e3;
c19 =  8.460780175632090e1;
c20 = -3.233571307223379e2;
c21 = -2.031538422351553e2;
c22 =  4.351585544019463e1;
c23 = -8.062279018001309e2;
c24 =  7.510134932437941e2;
c25 =  1.797443329095446e2;
c26 = -2.389853928747630e1;
c27 =  1.021046205356775e2;

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
x = sqrt(SA.*sfac);
y = CT.*0.025;

latentheat_evap = c0 + x.*(c1 + c4*y + x.*(c3 + y.*(c7 + c12*y) + x.*(c6 ...
    + y.*(c11 + y.*(c17 + c24*y)) + x.*(c10 + y.*(c16 + c23*y) + x.*(c15 ...
    + c22*y + c21*x))))) + y.*(c2 + y.*(c5 + c8*x + y.*(c9 + x.*(c13 ...
    + c18*x) + y.*(c14 + x.*(c19 + c25*x) + y.*(c20 + c26*x + c27*y)))));

if transposed
    latentheat_evap = latentheat_evap.';
end

end
