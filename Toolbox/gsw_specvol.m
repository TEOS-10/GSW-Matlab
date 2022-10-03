function specvol = gsw_specvol(SA,CT,p)

% gsw_specvol                            specific volume (75-term equation)
%==========================================================================
% 
% USAGE:  
%  specvol = gsw_specvol(SA,CT,p)
% 
% DESCRIPTION:
%  Calculates specific volume from Absolute Salinity, Conservative 
%  Temperature and pressure, using the computationally-efficient 75-term 
%  polynomial expression for specific volume (Roquet et al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if 
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
%  specvol  =  specific volume                                   [ m^3/kg ]
%
% AUTHOR: 
%  Fabien Roquet, Gurvan Madec, Trevor McDougall & Paul Barker
%                                                      [ help@teos-10.org ]
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
%   temperature and density of seawater.  J. Atmos. Ocean. Tech., 20, 
%   730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
% 
% ==========================================================================

% --------------------------------------------------------------------------
% Check variables and resize if necessary
% --------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_specvol:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_specvol: SA and CT do not have the same dimensions')
end %if

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p.*ones(size(SA));
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
    error('gsw_specvol: The dimensions of p do not agree')
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

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35))
                                             % deltaSA = 24 g/kg
offset = 5.971840214030754e-1;               % offset = deltaSA*sfac

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;
z = p.*1e-4;

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

specvol = v000 + xs.*(v100 + xs.*(v200 + xs.*(v300 + xs.*(v400 + xs.*(v500 + xs.*v600))))) ...
   + ys.*(v010 + xs.*(v110 + xs.*(v210 + xs.*(v310 + xs.*(v410 + xs.*v510)))) ...
   + ys.*(v020 + xs.*(v120 + xs.*(v220 + xs.*(v320 + xs.*v420))) ...
   + ys.*(v030 + xs.*(v130 + xs.*(v230 + xs.*v330)) ...
   + ys.*(v040 + xs.*(v140 + xs.*v240) ...
   + ys.*(v050 + xs.*v150 ...
   + ys.*v060))))) ...
+ z.*(    v001 + xs.*(v101 + xs.*(v201 + xs.*(v301 + xs.*(v401 + xs.*v501)))) ...
   + ys.*(v011 + xs.*(v111 + xs.*(v211 + xs.*(v311 + xs.*v411))) ...
   + ys.*(v021 + xs.*(v121 + xs.*(v221 + xs.*v321)) ....
   + ys.*(v031 + xs.*(v131 + xs.*v231) ...
   + ys.*(v041 + xs.*v141 ...
   + ys.*v051)))) ...
+ z.*(    v002 + xs.*(v102 + xs.*(v202 + xs.*(v302 + xs.*v402))) ...
   + ys.*(v012 + xs.*(v112 + xs.*(v212 + xs.*v312)) ...
   + ys.*(v022 + xs.*(v122 + xs.*v222) ...
   + ys.*(v032 + xs.*v132 ...
   + ys.*v042))) ...
+ z.*(    v003 + xs.*(v103 + xs.*v203) ...
   + ys.*(v013 + xs.*v113 ...
   + ys.*v023) ...
+ z.*(    v004 + xs.*v104 ...
   + ys.*v014 ...  
+ z.*(    v005 ...
+ z.*v006)))));

%--------------------------------------------------------------------------
% This function calculates specific volume using the computationally
% efficient 75-term expression for specific volume in terms of SA, CT and
% p. If one wanted to compute specvol from SA, CT, and p with the full 
% TEOS-10 Gibbs function, the following lines of code will enable this.
%
%    t = gsw_t_from_CT(SA,CT,p);
%    specvol = gsw_specvol_t_exact(SA,t,p);
%
%    or call the following, it is identical to the lines above.
%
%    specvol = gsw_specvol_CT_exact(SA,CT,p)
%
%-----------------This is the end of the alternative code------------------

if transposed
     specvol = specvol.';
end

end
