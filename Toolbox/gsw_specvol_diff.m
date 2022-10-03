function specvol_diff = gsw_specvol_diff(SA,CT,p_shallow,p_deep)

% gsw_specvol_diff         specific volume difference between two pressures
%                                                        (75-term equation)
%==========================================================================
% 
% USAGE:  
%  specvol_diff = gsw_specvol_diff(SA,CT,p_shallow,p_deep)
%
% DESCRIPTION:
%  Calculates the difference of the specific volume of seawater between 
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  This 
%  function uses the computationally-efficient expression for specific 
%  volume in terms of SA, CT and p (Roquet et al., 2015).  The output
%  (specvol_diff) is the specific volume evaluated at (SA,CT,p_deep)
%  minus the specific volume at (SA,CT,p_shallow). 
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA     =  Absolute Salinity                                     [ g/kg ]
%  CT     =  Conservative Temperature (ITS-90)                    [ deg C ]
%  p_shallow  =  upper sea pressure                                [ dbar ]
%                ( i.e. shallower absolute pressure - 10.1325 dbar ) 
%  p_deep     =  lower sea pressure                                [ dbar ]
%                ( i.e. deeper absolute pressure - 10.1325 dbar )
%
%  SA & CT  need to have the same dimensions.
%  p_shallow and p_deep may have dimensions Mx1 or 1xN or MxN, 
%  where SA and CT are MxN.
%
% OUTPUT:
%  specvol_diff  =  difference of specific volume                [ m^3/kg ]
%                       (deep minus shallow)                      
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
%    See Eqn. (3.7.3) of this TEOS-10 Manual. 
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
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
    error('gsw_specvol_diff: requires four inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mpu,npu] = size(p_shallow);
[mpl,npl] = size(p_deep);

if (ms~=mt) | (ns~=nt)
    error('gsw_specvol_diff: SA & CT need to have the same dimensions')
end

if (mpu == 1) & (npu == 1)                     
    p_shallow = p_shallow*ones(size(SA));
elseif (ns == npu) & (mpu == 1)              
    p_shallow = p_shallow(ones(1,ms), :);         
elseif (ms == mpu) & (npu == 1)              
    p_shallow = p_shallow(:,ones(1,ns));        
elseif (ns == mpu) & (npu == 1)          
    p_shallow = p_shallow.';                           
    p_shallow = p_shallow(ones(1,ms), :);               
elseif (ms == mpu) & (ns == npu)
    % ok
end

if (mpl == 1) & (npl == 1)                      
    p_deep = p_deep*ones(size(SA));
elseif (ns == npl) & (mpl == 1)                    
    p_deep = p_deep(ones(1,ms), :);               
elseif (ms == mpl) & (npl == 1)                  
    p_deep = p_deep(:,ones(1,ns));                 
elseif (ns == mpl) & (npl == 1)         
    p_deep = p_deep.';                              
    p_deep = p_deep(ones(1,ms), :);            
elseif (ms == mpl) & (ns == npl)
    % ok
else
    error('gsw_specvol_diff: Inputs array dimensions arguments do not agree')
end

if ms == 1
    SA = SA.';
    CT = CT.';
    p_shallow = p_shallow.';
    p_deep = p_deep.';
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
z_shallow = p_shallow.*1e-4;
z_deep = p_deep.*1e-4;

%v000 =  1.0769995862e-3; 
v001 = -6.0799143809e-5; 
v002 =  9.9856169219e-6; 
v003 = -1.1309361437e-6; 
v004 =  1.0531153080e-7; 
v005 = -1.2647261286e-8; 
v006 =  1.9613503930e-9; 
%v100 = -3.1038981976e-4; 
v101 =  2.4262468747e-5; 
v102 = -5.8484432984e-7; 
v103 =  3.6310188515e-7; 
v104 = -1.1147125423e-7; 
%v200 =  6.6928067038e-4; 
v201 = -3.4792460974e-5; 
v202 = -4.8122251597e-6; 
v203 =  1.6746303780e-8; 
%v300 = -8.5047933937e-4; 
v301 =  3.7470777305e-5; 
v302 =  4.9263106998e-6; 
%v400 =  5.8086069943e-4; 
v401 = -1.7322218612e-5; 
v402 = -1.7811974727e-6; 
%v500 = -2.1092370507e-4; 
v501 =  3.0927427253e-6; 
%v600 =  3.1932457305e-5; 
%v010 = -1.5649734675e-5; 
v011 =  1.8505765429e-5; 
v012 = -1.1736386731e-6; 
v013 = -3.6527006553e-7; 
v014 =  3.1454099902e-7; 
%v110 =  3.5009599764e-5; 
v111 = -9.5677088156e-6; 
v112 = -5.5699154557e-6; 
v113 = -2.7295696237e-7; 
%v210 = -4.3592678561e-5; 
v211 =  1.1100834765e-5; 
v212 =  5.4620748834e-6; 
%v310 =  3.4532461828e-5; 
v311 = -9.8447117844e-6; 
v312 = -1.3544185627e-6; 
%v410 = -1.1959409788e-5; 
v411 =  2.5909225260e-6; 
%v510 =  1.3864594581e-6;
%v020 =  2.7762106484e-5; 
v021 = -1.1716606853e-5; 
v022 =  2.1305028740e-6; 
v023 =  2.8695905159e-7; 
%v120 = -3.7435842344e-5; 
v121 = -2.3678308361e-7; 
v122 =  3.9137387080e-7; 
%v220 =  3.5907822760e-5; 
v221 =  2.9283346295e-6; 
v222 = -6.5731104067e-7; 
%v320 = -1.8698584187e-5; 
v321 = -4.8826139200e-7; 
%v420 =  3.8595339244e-6; 
%v030 = -1.6521159259e-5; 
v031 =  7.9279656173e-6; 
v032 = -4.6132540037e-7; 
%v130 =  2.4141479483e-5; 
v131 = -3.4558773655e-6; 
v132 =  7.7618888092e-9; 
%v230 = -1.4353633048e-5; 
v231 =  3.1655306078e-7; 
%v330 =  2.2863324556e-6;
%v040 =  6.9111322702e-6; 
v041 = -3.4102187482e-6; 
v042 = -6.3352916514e-8; 
%v140 = -8.7595873154e-6; 
v141 =  1.2956717783e-6; 
%v240 =  4.3703680598e-6; 
%v050 = -8.0539615540e-7; 
v051 =  5.0736766814e-7; 
%v150 = -3.3052758900e-7; 
%v060 =  2.0543094268e-7;

xy_part_1 = v001 + xs.*(v101 + xs.*(v201 + xs.*(v301 + xs.*(v401 + v501.*xs)))) ...
  + ys.*(v011 + xs.*(v111 + xs.*(v211 + xs.*(v311 + v411.*xs))) + ys.*(v021 ...
  + xs.*(v121 + xs.*(v221 + v321.*xs)) + ys.*(v031 + xs.*(v131 + v231.*xs) ...
  + ys.*(v041 + v141.*xs + v051.*ys))));

xy_part_2 = v002 + xs.*(v102 + xs.*(v202 + xs.*(v302 + v402.*xs))) ...
    + ys.*(v012 + xs.*(v112 + xs.*(v212 + v312.*xs)) + ys.*(v022 + xs.*(v122 ...
    + v222.*xs) + ys.*(v032 + v132.*xs + v042.*ys)));

xy_part_3 = v003 + xs.*(v103 + v203.*xs) + ys.*(v013 + v113.*xs + v023.*ys);
    
xy_part_4 = v004 + v104.*xs + v014.*ys;

dz = z_deep - z_shallow;

z_part_2 = z_deep + z_shallow;

z_deep2 = z_deep.*z_deep;
z_shallow2 = z_shallow.*z_shallow;
z_deep_times_z_shallow = z_deep.*z_shallow;
z_deep2_plus_z_shallow2 = z_deep2 + z_shallow2;

z_part_3 = z_deep2_plus_z_shallow2 + z_deep_times_z_shallow;

z_part_4 = z_part_2.*z_deep2_plus_z_shallow2;

z_shallow3 = z_shallow.*z_shallow2;

z_part_5 = z_deep2.*z_part_3 + z_shallow3.*z_part_2;

z_part_6 = z_part_2.*z_part_3.*(z_deep2_plus_z_shallow2 - z_deep_times_z_shallow);

specvol_diff = dz.*(xy_part_1 + z_part_2.*xy_part_2 + z_part_3.*xy_part_3 ...
    + z_part_4.*xy_part_4 + z_part_5.*v005 + z_part_6.*v006);

%--------------------------------------------------------------------------
% This function calculates specvol_diff using the computationally
% efficient expression for specific volume in terms of SA, CT and p. 
% If one wanted to compute the specific volume difference using the full  
% TEOS-10 Gibbs function, the following lines of code will enable this.
%
%    pt = gsw_pt_from_CT(SA,CT);
%    pr0 = zeros(size(SA)); 
%    t_shallow = gsw_pt_from_t(SA,pt,pr0,p_shallow);
%    t_deep = gsw_pt_from_t(SA,pt,pr0,p_deep);
%    specvol_diff = gsw_specvol_t_exact(SA,t_deep,p_deep) - ...
%                    gsw_specvol_t_exact(SA,t_shallow,p_shallow);
%
%    or call the following, it is identical to the lines above.
%
%    specvol_diff = gsw_specvol_diff_CT_exact(SA,CT,p_shallow,p_deep)
%
%-----------------This is the end of the alternative code------------------

if transposed
    specvol_diff = specvol_diff.';
end

end
