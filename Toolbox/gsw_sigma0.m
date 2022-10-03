function sigma0 = gsw_sigma0(SA,CT)

% gsw_sigma0                       potential density anomaly with reference
%                                 sea pressure of 0 dbar (75-term equation)
%==========================================================================
% 
% USAGE:  
%  sigma0 = gsw_sigma0(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 0 dbar,
%  this being this particular potential density minus 1000 kg/m^3.  This
%  function has inputs of Absolute Salinity and Conservative Temperature.
%  This function uses the computationally-efficient expression for 
%  specific volume in terms of SA, CT and p (Roquet et al., 2015).
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
%  sigma0  =  potential density anomaly with                     [ kg/m^3 ]
%             respect to a reference pressure of 0 dbar,   
%             that is, this potential density - 1000 kg/m^3.
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
%    See Eqn. (A.30.1) of this TEOS-10 Manual. 
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
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_sigma0:  Requires two inputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_sigma0: SA and CT must have same dimensions')
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

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

%deltaS = 24;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
offset = 5.971840214030754e-1;                      % offset = deltaS*sfac.

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;

v000 =  1.0769995862e-3; 
% v001 = -6.0799143809e-5; 
% v002 =  9.9856169219e-6; 
% v003 = -1.1309361437e-6; 
% v004 =  1.0531153080e-7; 
% v005 = -1.2647261286e-8; 
% v006 =  1.9613503930e-9; 
v010 = -1.5649734675e-5; 
% v011 =  1.8505765429e-5; 
% v012 = -1.1736386731e-6; 
% v013 = -3.6527006553e-7; 
% v014 =  3.1454099902e-7; 
v020 =  2.7762106484e-5; 
% v021 = -1.1716606853e-5; 
% v022 =  2.1305028740e-6; 
% v023 =  2.8695905159e-7; 
v030 = -1.6521159259e-5; 
% v031 =  7.9279656173e-6; 
% v032 = -4.6132540037e-7; 
v040 =  6.9111322702e-6; 
% v041 = -3.4102187482e-6; 
% v042 = -6.3352916514e-8; 
v050 = -8.0539615540e-7; 
% v051 =  5.0736766814e-7; 
v060 =  2.0543094268e-7;
v100 = -3.1038981976e-4; 
% v101 =  2.4262468747e-5; 
% v102 = -5.8484432984e-7; 
% v103 =  3.6310188515e-7; 
% v104 = -1.1147125423e-7;
v110 =  3.5009599764e-5; 
% v111 = -9.5677088156e-6; 
% v112 = -5.5699154557e-6; 
% v113 = -2.7295696237e-7; 
v120 = -3.7435842344e-5; 
% v121 = -2.3678308361e-7; 
% v122 =  3.9137387080e-7; 
v130 =  2.4141479483e-5; 
% v131 = -3.4558773655e-6; 
% v132 =  7.7618888092e-9; 
v140 = -8.7595873154e-6; 
% v141 =  1.2956717783e-6; 
v150 = -3.3052758900e-7; 
v200 =  6.6928067038e-4; 
% v201 = -3.4792460974e-5; 
% v202 = -4.8122251597e-6; 
% v203 =  1.6746303780e-8; 
v210 = -4.3592678561e-5; 
% v211 =  1.1100834765e-5; 
% v212 =  5.4620748834e-6; 
v220 =  3.5907822760e-5; 
% v221 =  2.9283346295e-6; 
% v222 = -6.5731104067e-7; 
v230 = -1.4353633048e-5; 
% v231 =  3.1655306078e-7; 
v240 =  4.3703680598e-6; 
v300 = -8.5047933937e-4; 
% v301 =  3.7470777305e-5; 
% v302 =  4.9263106998e-6; 
v310 =  3.4532461828e-5; 
% v311 = -9.8447117844e-6; 
% v312 = -1.3544185627e-6; 
v320 = -1.8698584187e-5; 
% v321 = -4.8826139200e-7; 
v330 =  2.2863324556e-6;
v400 =  5.8086069943e-4; 
% v401 = -1.7322218612e-5; 
% v402 = -1.7811974727e-6; 
v410 = -1.1959409788e-5; 
% v411 =  2.5909225260e-6; 
v420 =  3.8595339244e-6; 
v500 = -2.1092370507e-4; 
% v501 =  3.0927427253e-6; 
v510 =  1.3864594581e-6;
v600 =  3.1932457305e-5; 


vp0 = v000 + xs.*(v100 + xs.*(v200 + xs.*(v300 + xs.*(v400 + xs.*(v500 ...
    + v600.*xs))))) + ys.*(v010 + xs.*(v110 + xs.*(v210 + xs.*(v310 + xs.*(v410 ...
    + v510.*xs)))) + ys.*(v020 + xs.*(v120 + xs.*(v220 + xs.*(v320 + v420.*xs))) ...
    + ys.*(v030 + xs.*(v130 + xs.*(v230 + v330.*xs)) + ys.*(v040 + xs.*(v140 ...
    + v240*xs) + ys.*(v050 + v150.*xs + v060.*ys)))));

sigma0 = 1./vp0 - 1000;

%--------------------------------------------------------------------------
% This function calculates sigma0 using the computationally-efficient 
% expression for specific volume in terms of SA, CT and p.  If one wanted 
% to compute sigma0 with the full TEOS-10 Gibbs function expression for 
% specific volume, the following lines of code will enable this.
%
%  pr0 = zeros(size(SA));
%  sigma0 = gsw_rho_CT_exact(SA,CT,pr0) - 1000;
%
%---------------This is the end of the alternative code--------------------

if transposed
    sigma0 = sigma0.';
end

% The output, being potential density anomaly, has units of kg/m^3 and is 
% potential density with 1000 kg/m^3 subtracted from it. 

end
