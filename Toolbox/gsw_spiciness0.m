function spiciness0 = gsw_spiciness0(SA,CT)

% gsw_spiciness0                                    spiciness at p = 0 dbar
%                                                        (75-term equation)
%==========================================================================
% 
% USAGE:  
%  spiciness0 = gsw_spiciness0(SA,CT)
%
% DESCRIPTION:
%  Calculates spiciness from Absolute Salinity and Conservative 
%  Temperature at a pressure of 0 dbar, as described by McDougall and 
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
%  spiciness0  =  spiciness referenced to a pressure of 0 dbar, 
%                 i.e. the surface                               [ kg/m^3 ]
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
   error('gsw_spiciness0:  Requires two inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_spiciness0: SA and CT must have same dimensions')
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

s01 = -9.22982898371678e1;
s02 = -1.35727873628866e1;
s03 =  1.87353650994010e1;
s04 = -1.61360047373455e1;
s05 =  3.76112762286425e1;
s06 = -4.27086671461257e1;
s07 =  2.00820111041594e1;
s08 =  2.87969717584045e2;
s09 =  1.13747111959674e1;
s10 =  6.07377192990680e1;
s11 = -7.37514033570187e1;
s12 = -7.51171878953574e1;
s13 =  1.63310989721504e2;
s14 = -8.83222751638095e1;
s15 = -6.41725302237048e2;
s16 =  2.79732530789261e1;
s17 = -2.49466901993728e2;
s18 =  3.26691295035416e2;
s19 =  2.66389243708181e1;
s20 = -2.93170905757579e2;
s21 =  1.76053907144524e2;
s22 =  8.27634318120224e2;
s23 = -7.02156220126926e1;
s24 =  3.82973336590803e2;
s25 = -5.06206828083959e2;
s26 =  6.69626565169529e1;
s27 =  3.02851235050766e2;
s28 = -1.96345285604621e2;
s29 = -5.74040806713526e2;
s30 =  7.03285905478333e1;
s31 = -2.97870298879716e2;
s32 =  3.88340373735118e2;
s33 = -8.29188936089122e1;
s34 = -1.87602137195354e2;
s35 =  1.27096944425793e2;
s36 =  2.11671167892147e2;
s37 = -3.15140919876285e1;
s38 =  1.16458864953602e2;
s39 = -1.50029730802344e2;
s40 =  3.76293848660589e1;
s41 =  6.47247424373200e1;
s42 = -4.47159994408867e1;
s43 = -3.23533339449055e1;
s44 =  5.30648562097667;
s45 = -1.82051249177948e1;
s46 =  2.33184351090495e1;
s47 = -6.22909903460368;
s48 = -9.55975464301446;
s49 =  6.61877073960113;
 
spiciness0 = s01 + ys.*(s02 + ys.*(s03 + ys.*(s04 + ys.*(s05 + ys.*(s06 + s07*ys))))) ...
    + xs.*(s08 + ys.*(s09 + ys.*(s10 + ys.*(s11 + ys.*(s12 + ys.*(s13 + s14*ys)))))...
    + xs.*(s15 + ys.*(s16 + ys.*(s17 + ys.*(s18 + ys.*(s19 + ys.*(s20 + s21*ys))))) ...
    + xs.*(s22 + ys.*(s23 + ys.*(s24 + ys.*(s25 + ys.*(s26 + ys.*(s27 + s28*ys))))) ...
    + xs.*(s29 + ys.*(s30 + ys.*(s31 + ys.*(s32 + ys.*(s33 + ys.*(s34 + s35*ys))))) ...
    + xs.*(s36 + ys.*(s37 + ys.*(s38 + ys.*(s39 + ys.*(s40 + ys.*(s41 + s42*ys))))) ...
    + xs.*(s43 + ys.*(s44 + ys.*(s45 + ys.*(s46 + ys.*(s47 + ys.*(s48 + s49*ys)))))))))));

if transposed
    spiciness0 = spiciness0.';
end

end