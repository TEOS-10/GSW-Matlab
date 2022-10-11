function enthalpy_diff = gsw_enthalpy_diff(SA,CT,p_shallow,p_deep)

% gsw_enthalpy_diff            difference of enthalpy between two pressures
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  enthalpy_diff = gsw_enthalpy_diff(SA,CT,p_shallow,p_deep)
% 
% DESCRIPTION:
%  Calculates the difference of the specific enthalpy of seawater between 
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  This 
%  function uses the computationally-efficient expression for specific 
%  volume in terms of SA, CT and p (Roquet et al., 2015).  The output
%  (enthalpy_diff) is the specific enthalpy evaluated at (SA,CT,p_deep)
%  minus the specific enthalpy at (SA,CT,p_shallow). 
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA         =  Absolute Salinity                                 [ g/kg ]
%  CT         =  Conservative Temperature (ITS-90)                [ deg C ]
%  p_shallow  =  upper sea pressure                                [ dbar ]
%                ( i.e. shallower absolute pressure - 10.1325 dbar ) 
%  p_deep     =  lower sea pressure                                [ dbar ]
%                ( i.e. deeper absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p_shallow and p_deep may have dimensions Mx1 or 1xN or MxN, 
%  where SA and CT are MxN.
%
% OUTPUT:
%  enthalpy_diff  =  difference of specific enthalpy               [ J/kg ]
%                       (deep minus shallow)
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker.                     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (30th October, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) and (A.30.6) of this TEOS-10 Manual. 
%
%  McDougall, T.J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%    See Eqns. (18) and (22)
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

if ~(nargin == 4)
    error('gsw_enthalpy_diff: requires four inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mpu,npu] = size(p_shallow);
[mpl,npl] = size(p_deep);

if (ms~=mt) | (ns~=nt)
    error('gsw_enthalpy_diff: SA & CT need to have the same dimensions')
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
    error('gsw_enthalpy_diff: Inputs array dimensions arguments do not agree')
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

% % Set lower temperature limit for water that is approx 1C colder than the
% % freezing temperature.
% CT_frozen = -1.79e-2 - 6.7157e-2.*SA - 9.2708e-4.*p_deep;
% 
% if any(CT(:) < (CT_frozen(:) - 1))
%     [Icold] = find(CT < CT_frozen - 1);
%     CT(Icold) = NaN;
% end

%db2Pa = 1e4;                           % factor to convert from dbar to Pa

sfac = 0.0248826675584615;                    % sfac = 1/(40*(35.16504/35))
                                                        % deltaSA = 24 g/kg
offset = 5.971840214030754e-1;                       % offset = deltaS*sfac

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;
z_shallow = p_shallow.*1e-4;
z_deep = p_deep.*1e-4;

h001 =  1.0769995862e-3; 
h002 = -3.0399571905e-5; 
h003 =  3.3285389740e-6; 
h004 = -2.8273403593e-7; 
h005 =  2.1062306160e-8; 
h006 = -2.1078768810e-9; 
h007 =  2.8019291329e-10; 
h011 = -1.5649734675e-5; 
h012 =  9.2528827145e-6; 
h013 = -3.9121289103e-7; 
h014 = -9.1317516383e-8; 
h015 =  6.2908199804e-8; 
h021 =  2.7762106484e-5; 
h022 = -5.8583034265e-6; 
h023 =  7.1016762467e-7; 
h024 =  7.1739762898e-8; 
h031 = -1.6521159259e-5; 
h032 =  3.9639828087e-6; 
h033 = -1.5377513346e-7; 
h042 = -1.7051093741e-6; 
h043 = -2.1117638838e-8; 
h041 =  6.9111322702e-6; 
h051 = -8.0539615540e-7; 
h052 =  2.5368383407e-7; 
h061 =  2.0543094268e-7; 
h101 = -3.1038981976e-4; 
h102 =  1.21312343735e-5; 
h103 = -1.9494810995e-7; 
h104 =  9.0775471288e-8; 
h105 = -2.2294250846e-8; 
h111 =  3.5009599764e-5; 
h112 = -4.7838544078e-6; 
h113 = -1.8566384852e-6; 
h114 = -6.8239240593e-8; 
h121 = -3.7435842344e-5; 
h122 = -1.18391541805e-7; 
h123 =  1.3045795693e-7; 
h131 =  2.4141479483e-5; 
h132 = -1.72793868275e-6; 
h133 =  2.5872962697e-9; 
h141 = -8.7595873154e-6; 
h142 =  6.4783588915e-7; 
h151 = -3.3052758900e-7;
h201 =  6.6928067038e-4; 
h202 = -1.7396230487e-5; 
h203 = -1.6040750532e-6; 
h204 =  4.1865759450e-9; 
h211 = -4.3592678561e-5; 
h212 =  5.5504173825e-6; 
h213 =  1.8206916278e-6; 
h221 =  3.5907822760e-5; 
h222 =  1.46416731475e-6; 
h223 = -2.1910368022e-7; 
h231 = -1.4353633048e-5; 
h232 =  1.5827653039e-7; 
h241 =  4.3703680598e-6;
h301 = -8.5047933937e-4; 
h302 =  1.87353886525e-5; 
h303 =  1.6421035666e-6; 
h311 =  3.4532461828e-5; 
h312 = -4.9223558922e-6; 
h313 = -4.5147285423e-7; 
h321 = -1.8698584187e-5; 
h322 = -2.4413069600e-7; 
h331 =  2.2863324556e-6;
h401 =  5.8086069943e-4; 
h402 = -8.6611093060e-6; 
h403 = -5.9373249090e-7; 
h411 = -1.1959409788e-5; 
h421 =  3.8595339244e-6; 
h412 =  1.2954612630e-6;
h501 = -2.1092370507e-4; 
h502 =  1.54637136265e-6; 
h511 =  1.3864594581e-6; 
h601 =  3.1932457305e-5; 

xy_part_1 = h001 + xs.*(h101 + xs.*(h201 + xs.*(h301 + xs.*(h401 ...
    + xs.*(h501 + h601.*xs))))) + ys.*(h011 + xs.*(h111 + xs.*(h211 + xs.*(h311 ...
    + xs.*(h411 + h511.*xs)))) + ys.*(h021 + xs.*(h121 + xs.*(h221 + xs.*(h321 ...
    + h421.*xs))) + ys.*(h031 + xs.*(h131 + xs.*(h231 + h331.*xs)) + ys.*(h041 ...
    + xs.*(h141 + h241.*xs) + ys.*(h051 + h151.*xs + h061.*ys)))));

xy_part_2 = h002 + xs.*(h102 + xs.*(h202 + xs.*(h302 + xs.*(h402 + h502.*xs)))) ...
    + ys.*(h012 + xs.*(h112 + xs.*(h212 + xs.*(h312 + h412.*xs))) + ys.*(h022 ...
    + xs.*(h122 + xs.*(h222 + h322.*xs)) + ys.*(h032 + xs.*(h132 + h232.*xs) ...
    + ys.*(h042 + h142.*xs + h052.*ys))));

xy_part_3 = h003 + xs.*(h103 + xs.*(h203 + xs.*(h303 + h403.*xs))) + ys.*(h013 ...
    + xs.*(h113 + xs.*(h213 + h313.*xs)) + ys.*(h023 + xs.*(h123 + h223.*xs) ...
    + ys.*(h033 + h133.*xs + h043.*ys)));

xy_part_4 = h004 + xs.*(h104 + h204.*xs) + ys.*(h014 + h114.*xs + h024.*ys);

xy_part_5 = h005 + h105.*xs + h015.*ys;

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

z_part_7 = z_part_5.*z_deep2 + z_shallow2.*z_shallow3.*z_part_2;

enthalpy_diff = dz.*(xy_part_1 + z_part_2.*xy_part_2 + z_part_3.*xy_part_3 + z_part_4.*xy_part_4 ...
    + z_part_5.*xy_part_5 + z_part_6.*h006 + z_part_7.*h007).*1e8;

%--------------------------------------------------------------------------
% This function calculates enthalpy_diff using the computationally
% efficient expression for specific volume in terms of SA, CT and p.  If 
% one wanted to compute the enthalpy difference using the full TEOS-10 
% Gibbs function, the following lines of code will enable this.
%
%    t_shallow = gsw_t_from_CT(SA,CT,p_shallow);
%    t_deep = gsw_t_from_CT(SA,CT,p_deep);
%    enthalpy_diff = gsw_enthalpy_t_exact(SA,t_deep,p_deep) - ...
%                    gsw_enthalpy_t_exact(SA,t_shallow,p_shallow);
%
%    or call the following, it is identical to the lines above.
%
%    enthalpy_diff = gsw_enthalpy_diff_CT_exact(SA,CT,p_shallow,p_deep)
%
%-----------------This is the end of the alternative code------------------

if transposed
    enthalpy_diff = enthalpy_diff.';
end

end
