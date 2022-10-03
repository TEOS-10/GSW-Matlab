function [v_SA, v_CT, v_P] = gsw_specvol_first_derivatives(SA,CT,p)

% gsw_specvol_first_derivatives                     first order derivatives
%                                     of specific volume (75-term equation)
% =========================================================================
%
% USAGE:
%  [v_SA, v_CT, v_P] = gsw_specvol_first_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three first-order derivatives of specific
%  volume (v),
%   (1) v_SA, first-order derivative with respect to Absolute Salinity 
%       at constant CT & p.
%   (2) v_CT, first-order derivative with respect to CT at 
%       constant SA & p. 
%   (3) v_P, first-order derivative with respect to P at constant SA 
%       and CT. 
%
%  Note that this function uses the using the computationally-efficient
%  75-term expression for specific volume (Roquet et al., 2015).  There is 
%  an alternative to calling this function, namely 
%  gsw_specvol_first_derivatives_CT_exact(SA,CT,p) which uses the full 
%  Gibbs function (IOC et al., 2010).   
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described 
%  in McDougall et al. (2010).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  v_SA  =  The first derivative of specific volume with respect to 
%           Absolute Salinity at constant CT & p.     [ (m^3/kg)(g/kg)^-1 ]
%  v_CT  =  The first derivative of specific volume with respect to 
%           CT at constant SA and p.                         [ m^3/(K kg) ]
%  v_P   =  The first derivative of specific volume with respect to 
%           P at constant SA and CT.                        [ m^3/(Pa kg) ]
%
% AUTHOR:   
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_specvol_first_derivatives:  Requires three inputs')
end %if

if ~(nargout == 3)
   error('gsw_specvol_first_derivatives:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_specvol_first_derivatives: SA and CT do not have the same dimensions')
end %if

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_specvol_first_derivatives: The dimensions of p do not agree')
end %if

if ms == 1
    SA = SA.';
    CT = CT.';
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

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
offset = 5.971840214030754e-1;                      % offset = deltaS*sfac.

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;
z = p.*1e-4;

a000 = -1.5649734675e-5; 
a001 =  1.8505765429e-5; 
a002 = -1.1736386731e-6; 
a003 = -3.6527006553e-7; 
a004 =  3.1454099902e-7; 
a010 =  5.5524212968e-5; 
a011 = -2.3433213706e-5; 
a012 =  4.2610057480e-6; 
a013 =  5.7391810318e-7; 
a020 = -4.9563477777e-5; 
a021 =  2.37838968519e-5; 
a022 = -1.38397620111e-6; 
a030 =  2.76445290808e-5; 
a031 = -1.36408749928e-5; 
a032 = -2.53411666056e-7; 
a040 = -4.0269807770e-6; 
a041 =  2.5368383407e-6; 
a050 =  1.23258565608e-6; 
a100 =  3.5009599764e-5; 
a101 = -9.5677088156e-6; 
a102 = -5.5699154557e-6; 
a103 = -2.7295696237e-7; 
a110 = -7.4871684688e-5; 
a111 = -4.7356616722e-7; 
a112 =  7.8274774160e-7; 
a120 =  7.2424438449e-5; 
a121 = -1.03676320965e-5; 
a122 =  2.32856664276e-8; 
a130 = -3.50383492616e-5; 
a131 =  5.1826871132e-6; 
a140 = -1.6526379450e-6; 
a200 = -4.3592678561e-5; 
a201 =  1.1100834765e-5; 
a202 =  5.4620748834e-6; 
a210 =  7.1815645520e-5; 
a211 =  5.8566692590e-6; 
a212 = -1.31462208134e-6; 
a220 = -4.3060899144e-5; 
a221 =  9.4965918234e-7; 
a230 =  1.74814722392e-5; 
a300 =  3.4532461828e-5; 
a301 = -9.8447117844e-6; 
a302 = -1.3544185627e-6; 
a310 = -3.7397168374e-5; 
a311 = -9.7652278400e-7; 
a320 =  6.8589973668e-6; 
a400 = -1.1959409788e-5; 
a401 =  2.5909225260e-6; 
a410 =  7.7190678488e-6; 
a500 =  1.3864594581e-6; 

b000 = -3.1038981976e-4; 
b003 =  3.6310188515e-7; 
b004 = -1.1147125423e-7; 
b010 =  3.5009599764e-5; 
b013 = -2.7295696237e-7; 
b020 = -3.7435842344e-5; 
b030 =  2.4141479483e-5; 
b040 = -8.7595873154e-6; 
b050 = -3.3052758900e-7; 
b100 =  1.33856134076e-3; 
b103 =  3.3492607560e-8; 
b110 = -8.7185357122e-5; 
b120 =  7.1815645520e-5; 
b130 = -2.8707266096e-5; 
b140 =  8.7407361196e-6; 
b200 = -2.55143801811e-3; 
b210 =  1.03597385484e-4; 
b220 = -5.6095752561e-5; 
b230 =  6.8589973668e-6; 
b300 =  2.32344279772e-3; 
b310 = -4.7837639152e-5; 
b320 =  1.54381356976e-5; 
b400 = -1.05461852535e-3; 
b410 =  6.9322972905e-6; 
b500 =  1.9159474383e-4; 
b001 =  2.4262468747e-5; 
b011 = -9.5677088156e-6; 
b021 = -2.3678308361e-7; 
b031 = -3.4558773655e-6; 
b041 =  1.2956717783e-6; 
b101 = -6.9584921948e-5; 
b111 =  2.2201669530e-5; 
b121 =  5.8566692590e-6; 
b131 =  6.3310612156e-7; 
b201 =  1.12412331915e-4; 
b211 = -2.95341353532e-5; 
b221 = -1.4647841760e-6; 
b301 = -6.9288874448e-5; 
b311 =  1.0363690104e-5; 
b401 =  1.54637136265e-5; 
b002 = -5.8484432984e-7; 
b012 = -5.5699154557e-6; 
b022 =  3.9137387080e-7; 
b032 =  7.7618888092e-9; 
b102 = -9.62445031940e-6; 
b112 =  1.09241497668e-5; 
b122 = -1.31462208134e-6; 
b202 =  1.47789320994e-5; 
b212 = -4.0632556881e-6; 
b302 = -7.1247898908e-6; 

c000 = -6.0799143809e-5; 
c001 =  1.99712338438e-5; 
c002 = -3.3928084311e-6; 
c003 =  4.2124612320e-7; 
c004 = -6.3236306430e-8; 
c005 =  1.1768102358e-8; 
c010 =  1.8505765429e-5; 
c011 = -2.3472773462e-6; 
c012 = -1.09581019659e-6; 
c013 =  1.25816399608e-6; 
c020 = -1.1716606853e-5; 
c021 =  4.2610057480e-6; 
c022 =  8.6087715477e-7; 
c030 =  7.9279656173e-6; 
c031 = -9.2265080074e-7; 
c040 = -3.4102187482e-6; 
c041 = -1.26705833028e-7; 
c050 =  5.0736766814e-7; 
c100 =  2.4262468747e-5; 
c101 = -1.16968865968e-6; 
c102 =  1.08930565545e-6; 
c103 = -4.4588501692e-7; 
c110 = -9.5677088156e-6; 
c111 = -1.11398309114e-5; 
c112 = -8.1887088711e-7; 
c120 = -2.3678308361e-7; 
c121 =  7.8274774160e-7; 
c130 = -3.4558773655e-6; 
c131 =  1.55237776184e-8; 
c140 =  1.2956717783e-6; 
c200 = -3.4792460974e-5; 
c201 = -9.6244503194e-6; 
c202 =  5.0238911340e-8; 
c210 =  1.1100834765e-5; 
c211 =  1.09241497668e-5; 
c220 =  2.9283346295e-6; 
c221 = -1.31462208134e-6; 
c230 =  3.1655306078e-7; 
c300 =  3.7470777305e-5; 
c301 =  9.8526213996e-6; 
c310 = -9.8447117844e-6; 
c311 = -2.7088371254e-6; 
c320 = -4.8826139200e-7; 
c400 = -1.7322218612e-5; 
c401 = -3.5623949454e-6; 
c410 =  2.5909225260e-6; 
c500 =  3.0927427253e-6; 

v_CT_part = a000 + xs.*(a100 + xs.*(a200 + xs.*(a300 + xs.*(a400 + a500.*xs)))) ...
    + ys.*(a010 + xs.*(a110 + xs.*(a210 + xs.*(a310 + a410.*xs))) ...
    + ys.*(a020 + xs.*(a120 + xs.*(a220 + a320.*xs)) + ys.*(a030 ...
    + xs.*(a130 + a230.*xs) + ys.*(a040 + a140.*xs + a050.*ys )))) ...
    + z.*(a001 + xs.*(a101 + xs.*(a201 + xs.*(a301 + a401.*xs))) ...
    + ys.*(a011 + xs.*(a111 + xs.*(a211 + a311.*xs)) + ys.*(a021 ...
    + xs.*(a121 + a221.*xs) + ys.*(a031 + a131.*xs + a041.*ys))) ...
    + z.*(a002 + xs.*(a102 + xs.*(a202 + a302.*xs)) + ys.*(a012 ...
    + xs.*(a112 + a212.*xs) + ys.*(a022 + a122.*xs + a032.*ys)) ...
    + z.*(a003 + a103.*xs + a013.*ys + a004.*z))) ;

v_CT = 0.025.*v_CT_part;

v_SA_part = b000 + xs.*(b100 + xs.*(b200 + xs.*(b300 + xs.*(b400 + b500.*xs)))) ...
    + ys.*(b010 + xs.*(b110 + xs.*(b210 + xs.*(b310 + b410.*xs))) ...
    + ys.*(b020 + xs.*(b120 + xs.*(b220 + b320.*xs)) + ys.*(b030 ...
    + xs.*(b130 + b230.*xs) + ys.*(b040 + b140.*xs + b050.*ys)))) ...
    + z.*(b001 + xs.*(b101 + xs.*(b201 + xs.*(b301 + b401.*xs))) ...
    + ys.*(b011 + xs.*(b111 + xs.*(b211 + b311.*xs)) + ys.*(b021 ...
    + xs.*(b121 + b221.*xs) + ys.*(b031 + b131.*xs + b041.*ys))) ...
    + z.*(b002 + xs.*(b102 + xs.*(b202 + b302.*xs))+ ys.*(b012 ...
    + xs.*(b112 + b212.*xs) + ys.*(b022 + b122.*xs + b032.*ys)) ...
    + z.*(b003 +  b103.*xs + b013.*ys + b004.*z)));

v_SA = 0.5.*sfac.*v_SA_part./xs;

v_p_part = c000 + xs.*(c100 + xs.*(c200 + xs.*(c300 + xs.*(c400 + c500.*xs)))) ...
    + ys.*(c010 + xs.*(c110 + xs.*(c210 + xs.*(c310 + c410.*xs))) + ys.*(c020 ...
    + xs.*(c120 + xs.*(c220 + c320.*xs)) + ys.*(c030 + xs.*(c130 + c230*xs) ...
    + ys.*(c040 + c140.*xs + c050.*ys)))) + z.*(c001 + xs.*(c101 + xs.*(c201 ...
    + xs.*(c301 + c401.*xs))) + ys.*(c011 + xs.*(c111 + xs.*(c211 + c311.*xs)) ...
    + ys.*(c021 + xs.*(c121 + c221.*xs) + ys.*(c031 + c131.*xs + c041.*ys))) ...
    + z.*( c002 + xs.*(c102 + c202.*xs) + ys.*(c012 + c112.*xs + c022.*ys) ...
    + z.*(c003 + c103.*xs + c013.*ys + z.*(c004 + c005.*z))));

v_P = 1e-8.*v_p_part;

if transposed
    v_SA = v_SA.';
    v_CT = v_CT.';
    v_P = v_P.';
end

end
