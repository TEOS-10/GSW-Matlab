function spiciness2 = gsw_spiciness2(SA,CT)

% gsw_spiciness2                                 spiciness at p = 2000 dbar
%                                                        (75-term equation)
%==========================================================================
% 
% USAGE:  
%  spiciness2 = gsw_spiciness2(SA,CT,p)
%
% DESCRIPTION:
%  Calculates spiciness from Absolute Salinity and Conservative 
%  Temperature at a pressure of 2000 dbar, as described by McDougall and 
%  Krzysik (2015).  This routine is based on the computationally-efficient 
%  expression for specific volume in terms of SA, CT and p (Roquet et al., 
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  spiciness2  =  spiciness referenced to a pressure of 2000 dbar 
%                                                                [ kg/m^3 ]
%
% AUTHOR: 
%  Oliver Krzysik and Trevor McDougall                 [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  McDougall, T.J., and O.A. Krzysik, 2015: Spiciness. Journal of Marine 
%   Research, 73, 141-152.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%   
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_spiciness2:  Requires two inputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_spiciness2: SA and CT must have same dimensions')
end

if ms == 1
    SA = SA';
    CT = CT';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

%deltaS = 24;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
offset = 5.971840214030754e-1;                      % offset = deltaS*sfac.

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;

s01 = -9.17327320732265e1;
s02 = -1.31200235147912e1;
s03 =  2.49574345782503e1;
s04 = -2.41678075247398e1;
s05 =  3.61654631402053e1;
s06 = -3.22582164667710e1;
s07 =  1.45092623982509e1;
s08 =  2.87776645983195e2;
s09 =  3.13902307672447e1;
s10 =  1.69777467534459;
s11 = -5.69630115740438;
s12 = -7.97586359017987e1;
s13 =  1.07507460387751e2;
s14 = -5.58234404964787e1;
s15 = -6.41708068766557e2;
s16 = -2.53494801286161e1;
s17 = -9.86755437385364e1;
s18 =  1.52406930795842e2;
s19 =  4.23888258264105e1;
s20 = -1.60118811141438e2;
s21 =  9.67497898053989e1;
s22 =  8.27674355478637e2;
s23 =  5.27561234412133e-1;
s24 =  1.87440206992396e2;
s25 = -2.83295392345171e2;
s26 =  5.14485994597635e1;
s27 =  1.29975755062696e2;
s28 = -9.36526588377456e1;
s29 = -5.74911728972948e2;
s30 =  1.91175851862772e1;
s31 = -1.59347231968841e2;
s32 =  2.33884725744938e2;
s33 = -7.87744010546157e1;
s34 = -6.04757235443685e1;
s35 =  5.27869695599657e1;
s36 =  2.12517758478878e2;
s37 = -1.24351794740528e1;
s38 =  6.53904308937490e1;
s39 = -9.44804080763788e1;
s40 =  3.93874257887364e1;
s41 =  1.49425448888996e1;
s42 = -1.62350721656367e1;
s43 = -3.25936844276669e1;
s44 =  2.44035700301595;
s45 = -1.05079633683795e1;
s46 =  1.51515796259082e1;
s47 = -7.06609886460683;
s48 = -1.48043337052968;
s49 =  2.10066653978515;
 
spiciness2 = s01 + ys.*(s02 + ys.*(s03 + ys.*(s04 + ys.*(s05 + ys.*(s06 + s07*ys))))) ...
    + xs.*(s08 + ys.*(s09 + ys.*(s10 + ys.*(s11 + ys.*(s12 + ys.*(s13 + s14*ys)))))...
    + xs.*(s15 + ys.*(s16 + ys.*(s17 + ys.*(s18 + ys.*(s19 + ys.*(s20 + s21*ys))))) ...
    + xs.*(s22 + ys.*(s23 + ys.*(s24 + ys.*(s25 + ys.*(s26 + ys.*(s27 + s28*ys))))) ...
    + xs.*(s29 + ys.*(s30 + ys.*(s31 + ys.*(s32 + ys.*(s33 + ys.*(s34 + s35*ys))))) ...
    + xs.*(s36 + ys.*(s37 + ys.*(s38 + ys.*(s39 + ys.*(s40 + ys.*(s41 + s42*ys))))) ...
    + xs.*(s43 + ys.*(s44 + ys.*(s45 + ys.*(s46 + ys.*(s47 + ys.*(s48 + s49*ys)))))))))));

if transposed
    spiciness2 = spiciness2.';
end

end