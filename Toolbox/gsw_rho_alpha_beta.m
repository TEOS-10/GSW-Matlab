function [rho, alpha, beta] = gsw_rho_alpha_beta(SA,CT,p)

% gsw_rho_alpha_beta            in-situ density, thermal expansion & saline 
%                                contraction coefficient (75-term equation)
%==========================================================================
% 
% USAGE:  
% [rho, alpha, beta] = gsw_rho_alpha_beta(SA,CT,p)
%
% DESCRIPTION:
%  Calculates in-situ density, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from 
%  Absolute Salinity and Conservative Temperature.  This function uses the
%  computationally-efficient expression for specific volume in terms of 
%  SA, CT and p (Roquet et al., 2015).
%
%  Note that potential density (pot_rho) with respect to reference pressure
%  p_ref is obtained by calling this function with the pressure argument 
%  being p_ref as in [pot_rho, ~, ~] = gsw_rho_alpha_beta(SA,CT,p_ref).
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
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho    =  in-situ density                                       [ kg/m ]
%  alpha  =  thermal expansion coefficient                          [ 1/K ]
%            with respect to Conservative Temperature
%  beta   =  saline (i.e. haline) contraction                      [ kg/g ]
%            coefficient at constant Conservative Temperature
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual. 
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_rho_alpha_beta:  Requires three inputs')
end 

if ~(nargout == 3)
   error('gsw_rho_alpha_beta:  Requires three outputs')
end 

p = gsw_resize(p,size(SA));

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_rho_alpha_beta: SA and CT must have same dimensions')
end
% 
% if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
%     p = p*ones(size(SA));
% elseif (ns == np) & (mp == 1)         % p is row vector,
%     p = p(ones(1,ms), :);              % copy down each column.
% elseif (ms == mp) & (np == 1)         % p is column vector,
%     p = p(:,ones(1,ns));               % copy across each row.
% elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
%     p = p.';                              % transposed then
%     p = p(ones(1,ms), :);                % copy down each column.
% elseif (ms == mp) & (ns == np)
%     % ok
% else
%     error('gsw_rho_alpha_beta: Inputs array dimensions arguments do not agree')
% end %if

% if ms == 1
%     SA = SA.';
%     CT = CT.';
%     p = p.';
%     transposed = 1;
% else
%     transposed = 0;
% end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

%deltaS = 24;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
offset = 5.971840214030754e-1;                      % offset = deltaS*sfac.

x2 = sfac.*SA;
clear SA
xs = sqrt(x2 + offset);
ys = CT.*0.025;
clear CT
z = p.*1e-4;
clear p

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
b001 =  2.4262468747e-5; 
b002 = -5.8484432984e-7; 
b003 =  3.6310188515e-7; 
b004 = -1.1147125423e-7; 
b010 =  3.5009599764e-5; 
b011 = -9.5677088156e-6; 
b012 = -5.5699154557e-6; 
b013 = -2.7295696237e-7; 
b020 = -3.7435842344e-5; 
b021 = -2.3678308361e-7; 
b022 =  3.9137387080e-7; 
b030 =  2.4141479483e-5; 
b031 = -3.4558773655e-6; 
b032 =  7.7618888092e-9; 
b040 = -8.7595873154e-6; 
b041 =  1.2956717783e-6; 
b050 = -3.3052758900e-7; 
b100 =  1.33856134076e-3; 
b101 = -6.9584921948e-5; 
b102 = -9.62445031940e-6; 
b103 =  3.3492607560e-8; 
b110 = -8.7185357122e-5; 
b111 =  2.2201669530e-5; 
b112 =  1.09241497668e-5; 
b120 =  7.1815645520e-5; 
b121 =  5.8566692590e-6; 
b122 = -1.31462208134e-6; 
b130 = -2.8707266096e-5; 
b131 =  6.3310612156e-7; 
b140 =  8.7407361196e-6; 
b200 = -2.55143801811e-3; 
b201 =  1.12412331915e-4; 
b202 =  1.47789320994e-5; 
b210 =  1.03597385484e-4; 
b211 = -2.95341353532e-5; 
b212 = -4.0632556881e-6; 
b220 = -5.6095752561e-5; 
b221 = -1.4647841760e-6; 
b230 =  6.8589973668e-6; 
b300 =  2.32344279772e-3; 
b301 = -6.9288874448e-5; 
b302 = -7.1247898908e-6; 
b310 = -4.7837639152e-5; 
b311 =  1.0363690104e-5; 
b320 =  1.54381356976e-5; 
b400 = -1.05461852535e-3; 
b401 =  1.54637136265e-5; 
b410 =  6.9322972905e-6; 
b500 =  1.9159474383e-4; 

v000 =  1.0769995862e-3; 
v001 = -6.0799143809e-5; 
v002 =  9.9856169219e-6; 
v003 = -1.1309361437e-6; 
v004 =  1.0531153080e-7; 
v005 = -1.2647261286e-8; 
v006 =  1.9613503930e-9; 
v010 = -1.5649734675e-5; 
v011 =  1.8505765429e-5; 
v012 = -1.1736386731e-6; 
v013 = -3.6527006553e-7; 
v014 =  3.1454099902e-7; 
v020 =  2.7762106484e-5; 
v021 = -1.1716606853e-5; 
v022 =  2.1305028740e-6; 
v023 =  2.8695905159e-7; 
v030 = -1.6521159259e-5; 
v031 =  7.9279656173e-6; 
v032 = -4.6132540037e-7; 
v040 =  6.9111322702e-6; 
v041 = -3.4102187482e-6; 
v042 = -6.3352916514e-8; 
v050 = -8.0539615540e-7; 
v051 =  5.0736766814e-7; 
v060 =  2.0543094268e-7;
v100 = -3.1038981976e-4; 
v101 =  2.4262468747e-5; 
v102 = -5.8484432984e-7; 
v103 =  3.6310188515e-7; 
v104 = -1.1147125423e-7;
v110 =  3.5009599764e-5; 
v111 = -9.5677088156e-6; 
v112 = -5.5699154557e-6; 
v113 = -2.7295696237e-7; 
v120 = -3.7435842344e-5; 
v121 = -2.3678308361e-7; 
v122 =  3.9137387080e-7; 
v130 =  2.4141479483e-5; 
v131 = -3.4558773655e-6; 
v132 =  7.7618888092e-9; 
v140 = -8.7595873154e-6; 
v141 =  1.2956717783e-6; 
v150 = -3.3052758900e-7; 
v200 =  6.6928067038e-4; 
v201 = -3.4792460974e-5; 
v202 = -4.8122251597e-6; 
v203 =  1.6746303780e-8; 
v210 = -4.3592678561e-5; 
v211 =  1.1100834765e-5; 
v212 =  5.4620748834e-6; 
v220 =  3.5907822760e-5; 
v221 =  2.9283346295e-6; 
v222 = -6.5731104067e-7; 
v230 = -1.4353633048e-5; 
v231 =  3.1655306078e-7; 
v240 =  4.3703680598e-6; 
v300 = -8.5047933937e-4; 
v301 =  3.7470777305e-5; 
v302 =  4.9263106998e-6; 
v310 =  3.4532461828e-5; 
v311 = -9.8447117844e-6; 
v312 = -1.3544185627e-6; 
v320 = -1.8698584187e-5; 
v321 = -4.8826139200e-7; 
v330 =  2.2863324556e-6;
v400 =  5.8086069943e-4; 
v401 = -1.7322218612e-5; 
v402 = -1.7811974727e-6; 
v410 = -1.1959409788e-5; 
v411 =  2.5909225260e-6; 
v420 =  3.8595339244e-6; 
v500 = -2.1092370507e-4; 
v501 =  3.0927427253e-6; 
v510 =  1.3864594581e-6;
v600 =  3.1932457305e-5; 
      
v = v000 + xs.*(v100 + xs.*(v200 + xs.*(v300 + xs.*(v400 + xs.*(v500 ...
    + v600.*xs))))) + ys.*(v010 + xs.*(v110 + xs.*(v210 + xs.*(v310 + xs.*(v410 ...
    + v510.*xs)))) + ys.*(v020 + xs.*(v120 + xs.*(v220 + xs.*(v320 + v420.*xs))) ...
    + ys.*(v030 + xs.*(v130 + xs.*(v230 + v330.*xs)) + ys.*(v040 + xs.*(v140 ...
    + v240*xs) + ys.*(v050 + v150.*xs + v060.*ys))))) + z.*(v001 + xs.*(v101 ...
    + xs.*(v201 + xs.*(v301 + xs.*(v401 + v501.*xs)))) + ys.*(v011 + xs.*(v111 ...
    + xs.*(v211 + xs.*(v311 + v411.*xs))) + ys.*(v021 + xs.*(v121 + xs.*(v221 ...
    + v321.*xs)) + ys.*(v031 + xs.*(v131 + v231.*xs) + ys.*(v041 + v141.*xs ...
    + v051.*ys)))) + z.*(v002 + xs.*(v102 + xs.*(v202 + xs.*(v302 + v402.*xs))) ...
    + ys.*(v012 + xs.*(v112 + xs.*(v212 + v312.*xs)) + ys.*(v022 + xs.*(v122 ...
    + v222.*xs) + ys.*(v032 + v132.*xs + v042.*ys))) + z.*(v003 + xs.*(v103 ...
    + v203.*xs) + ys.*(v013 + v113.*xs + v023.*ys) + z.*(v004 + v104.*xs + v014.*ys ...
    + z.*(v005 + v006.*z)))));

rho = 1./v;

v_CT = a000 + xs.*(a100 + xs.*(a200 + xs.*(a300 + xs.*(a400 + a500.*xs)))) ...
    + ys.*(a010 + xs.*(a110 + xs.*(a210 + xs.*(a310 + a410.*xs))) ...
    + ys.*(a020 + xs.*(a120 + xs.*(a220 + a320.*xs)) + ys.*(a030 ...
    + xs.*(a130 + a230.*xs) + ys.*(a040 + a140.*xs + a050.*ys )))) ...
    + z.*(a001 + xs.*(a101 + xs.*(a201 + xs.*(a301 + a401.*xs))) ...
    + ys.*(a011 + xs.*(a111 + xs.*(a211 + a311.*xs)) + ys.*(a021 ...
    + xs.*(a121 + a221.*xs) + ys.*(a031 + a131.*xs + a041.*ys))) ...
    + z.*(a002 + xs.*(a102 + xs.*(a202 + a302.*xs)) + ys.*(a012 ...
    + xs.*(a112 + a212.*xs) + ys.*(a022 + a122.*xs + a032.*ys)) ...
    + z.*(a003 + a103.*xs + a013.*ys + a004.*z))) ;

alpha = 0.025.*v_CT./v;
clear v_CT

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
 
clear ys z

beta = -v_SA_part.*0.5.*sfac./(v.*xs);

%--------------------------------------------------------------------------
% This function calculates rho, alpha and beta using the computationally-
% efficient expression for specific volume in terms of SA, CT and p.  If 
% one wanted to compute rho, alpha and beta with the full TEOS-10 Gibbs 
% function expression for specific volume, the following lines of code will
% do this.
%
%    t = gsw_t_from_CT(SA,CT,p);
%    rho = gsw_rho_t_exact(SA,t,p);
%    alpha = gsw_alpha_wrt_CT_t_exact(SA,t,p);
%    beta = gsw_beta_const_CT_t_exact(SA,t,p);
%
%    or call the following, it is identical to the lines above.
%
%   [rho, alpha, beta] = gsw_rho_alpha_beta_CT_exact(SA,CT,p)
%
%--------------This is the end of the alternative code---------------------

% if transposed
%     rho = rho.';
%     alpha = alpha.';
%     beta = beta.';
% end

end
