function latentheat_melting = gsw_latentheat_melting(SA,p)

% gsw_latentheat_melting                             latent heat of melting 
%==========================================================================
%
% USAGE: 
%  latentheat_melting = gsw_latentheat_melting(SA,p)
%
% DESCRIPTION:
%  Calculates latent heat, or enthalpy, of melting.  It is defined in terms
%  of Absolute Salinity, SA, and sea pressure, p, and is valid in the  
%  ranges 0 < SA < 42 g kg^-1 and 0 < p < 10,000 dbar.  The errors range 
%  between -0.4 and 0.3 J kg^-1.
%  
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & p need to have the same dimensions.
%
% OUTPUT:
%  latentheat_melting  =  latent heat of melting                   [ J/kg ]     
%
% AUTHOR:  
%  Paul Barker, Trevor McDougall & Rainer Feistel      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See section 3.34 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_latentheat_melting:  Requires two input arguments')
end %if

[ms,ns] = size(SA);
[mp,np] = size(p);

if (mp ~= ms | np ~= ns)
    error('gsw_latentheat_melting: SA and p must have same dimensions')
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

c0 =  3.334265169240710e5;
c1 = -2.789444646733159;
c2 = -1.822150156453350e4;
c3 = -4.984585692734338e3;
c4 = -7.371966528571920e1;
c5 = -7.605802553358546e3;
c6 =  1.195857305019339e3;
c7 =  1.233720336206392e3;
c8 =  2.294798676591890e2;
c9 =  9.655751370889338e2;
c10 = -5.792068522727968e2;
c11 = -1.649446955902331e3;
c12 = -1.029021448430547e3;
c13 = -3.171558017172501e2;
c14 = -1.751401389905041e2;
c15 =  6.836527214265952e2;
c16 =  1.078283734113611e3;
c17 =  5.613896351265648e2;
c18 =  6.968934948667265e2;
c19 =  1.793032021946783e2;
c20 =  8.692558481134256e1;
c21 = -2.371103254714944e2;
c22 = -5.775033277201674e2;
c23 = -3.019749254648732e2;
c24 = -6.420420579160927e2;
c25 = -2.657570848596042e2;
c26 = -1.646738151143109e1;
c27 =  4.618228988300871;

S_u = 40*(35.16504/35);
x = sqrt(SA./S_u);
y = p.*1e-4;

latentheat_melting = c0 + x.*(c1 + c4*y + x.*(c3   ...
    + y.*(c7 + c12*y) + x.*(c6 + y.*(c11 + y.*(c17 + c24*y)) ...
    + x.*(c10  + y.*(c16 + c23*y) + x.*(c15 + c22*y + c21*x)))))  ...
    + y.*(c2 + y.*(c5 + c8*x + y.*(c9 + x.*(c13 + c18*x) ...
    + y.*(c14 + x.*(c19 + c25*x) + y.*(c20 + c26*x + c27*y)))));

% Note that the computed latent heat of melting from this function has 
% errors which range between -0.4 and 0.3 J kg^-1, when compared with the 
% latent heats of melting derived from the Gibbs functions of ice and of 
% seawater (using the SIA code of TEOS-10), however, the underlying data to
% the Gibbs function contains uncertainities of 200 J kg^-1 (IOC et al., 2010).  

if transposed
    latentheat_melting = latentheat_melting.';
end

end
