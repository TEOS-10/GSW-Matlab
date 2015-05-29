function enthalpy_diff = gsw_enthalpy_diff(SA,CT,p_shallow,p_deep)

% gsw_enthalpy_diff                 difference of enthalpy at two pressures
%                                                        (48-term equation)
%==========================================================================
%
% USAGE:
%  enthalpy_diff = gsw_enthalpy_diff(SA,CT,p_shallow,p_deep)
% 
% DESCRIPTION:
%  Calculates the difference of the specific enthalpy of seawater between 
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  This 
%  function uses the computationally-efficient 48-term expression for
%  density in terms of SA, CT and p (IOC et al., 2010).  The output
%  (enthalpy_diff_CT) is the specific enthalpy evaluated at (SA,CT,p_deep)
%  minus the specific enthalpy at (SA,CT,p_shallow). 
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2010).  The GSW library function 
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
%  p_shallow and p_deep may have dimensions Mx1 or 1xN or MxN, 
%  where SA and CT are MxN.
%
% OUTPUT:
%  enthalpy_diff_CT  =  difference of specific enthalpy            [ J/kg ]
%                       (deep minus shallow)
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker.                     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) and (A.30.6) of this TEOS-10 Manual. 
%
%  McDougall, T. J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%    See Eqns. (18) and (22)
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

if (mpu == 1) & (npu == 1)                           % p_shallow is a scalar 
    p_shallow = p_shallow*ones(size(SA));
elseif (ns == npu) & (mpu == 1)                  % p_shallow is row vector,
    p_shallow = p_shallow(ones(1,ms), :);          % copy down each column.
elseif (ms == mpu) & (npu == 1)               % p_shallow is column vector,
    p_shallow = p_shallow(:,ones(1,ns));            % copy across each row.
elseif (ns == mpu) & (npu == 1)          % p_shallow is a transposed row vector,
    p_shallow = p_shallow.';                              % transposed then
    p_shallow = p_shallow(ones(1,ms), :);                % copy down each column.
elseif (ms == mpu) & (ns == npu)
    % ok
end

if (mpl == 1) & (npl == 1)                             % p_deep is a scalar  
    p_deep = p_deep*ones(size(SA));
elseif (ns == npl) & (mpl == 1)                     % p_deep is row vector,
    p_deep = p_deep(ones(1,ms), :);                % copy down each column.
elseif (ms == mpl) & (npl == 1)                  % p_deep is column vector,
    p_deep = p_deep(:,ones(1,ns));                  % copy across each row.
elseif (ns == mpl) & (npl == 1)          % p_deep is a transposed row vector,
    p_deep = p_deep.';                              % transposed then
    p_deep = p_deep(ones(1,ms), :);                % copy down each column.
elseif (ms == mpl) & (ns == npl)
    % ok
else
    error('gsw_enthalpy_diff: Inputs array dimensions arguments do not agree')
end %if

if ~isempty(find(p_shallow - p_deep > 0));
 error('gsw_enthalpy_diff: p_deep needs to be greater than or equal to p_shallow')
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

db2Pa = 1e4;                      % factor to convert from dbar to Pa
%cp0 = 3991.86795711963;          % from Eqn. (3.3.3) of IOC et al. (2010).

v01 =  9.998420897506056e+2;
v02 =  2.839940833161907;
v03 = -3.147759265588511e-2;
v04 =  1.181805545074306e-3;
v05 = -6.698001071123802;
v06 = -2.986498947203215e-2;
v07 =  2.327859407479162e-4;
v08 = -3.988822378968490e-2;
v09 =  5.095422573880500e-4;
v10 = -1.426984671633621e-5;
v11 =  1.645039373682922e-7;
v12 = -2.233269627352527e-2;
v13 = -3.436090079851880e-4;
v14 =  3.726050720345733e-6;
v15 = -1.806789763745328e-4;
v16 =  6.876837219536232e-7;
v17 = -3.087032500374211e-7;
v18 = -1.988366587925593e-8;
v19 = -1.061519070296458e-11;
v20 =  1.550932729220080e-10;
v21 =  1.0;
v22 =  2.775927747785646e-3;
v23 = -2.349607444135925e-5;
v24 =  1.119513357486743e-6;
v25 =  6.743689325042773e-10;
v26 = -7.521448093615448e-3;
v27 = -2.764306979894411e-5;
v28 =  1.262937315098546e-7;
v29 =  9.527875081696435e-10;
v30 = -1.811147201949891e-11;
v31 = -3.303308871386421e-5;
v32 =  3.801564588876298e-7;
v33 = -7.672876869259043e-9;
v34 = -4.634182341116144e-11;
v35 =  2.681097235569143e-12;
v36 =  5.419326551148740e-6;
v37 = -2.742185394906099e-5;
v38 = -3.212746477974189e-7;
v39 =  3.191413910561627e-9;
v40 = -1.931012931541776e-12;
v41 = -1.105097577149576e-7;
v42 =  6.211426728363857e-10;
v43 = -1.119011592875110e-10;
v44 = -1.941660213148725e-11;
v45 = -1.864826425365600e-14;
v46 =  1.119522344879478e-14;
v47 = -1.200507748551599e-15;
v48 =  6.057902487546866e-17;

sqrtSA = sqrt(SA);

 a0 =  v21 + CT.*(v22 + CT.*(v23 + CT.*(v24 + v25*CT))) ...
     + SA.*(v26 + CT.*(v27 + CT.*(v28 + CT.*(v29 + v30*CT))) + v36*SA ...
 + sqrtSA.*(v31 + CT.*(v32 + CT.*(v33 + CT.*(v34 + v35*CT)))));

a1 = v37 + CT.*(v38 + CT.*(v39 + v40*CT)) + SA.*(v41 + v42*CT);

a2 = v43 + CT.*(v44 + v45*CT + v46*SA);

a3 = v47 + v48*CT;

b0 = v01 + CT.*(v02 + CT.*(v03 + v04*CT))  ...
         + SA.*(v05 + CT.*(v06 + v07*CT) ...
     + sqrtSA.*(v08 + CT.*(v09 + CT.*(v10 + v11*CT))));
 
b1 = 0.5*(v12 + CT.*(v13 + v14*CT) + SA.*(v15 + v16*CT));

b2 = v17 + CT.*(v18 + v19*CT) + v20*SA;

b1sq = b1.*b1; 
sqrt_disc = sqrt(b1sq - b0.*b2);

N = a0 + (2*a3.*b0.*b1./b2 - a2.*b0)./b2;

M = a1 + (4*a3.*b1sq./b2 - a3.*b0 - 2*a2.*b1)./b2;

A = b1 - sqrt_disc;
B = b1 + sqrt_disc;
delta_p = p_deep - p_shallow;
p_sum = p_deep + p_shallow;
part1 = b0 + p_shallow.*(2*b1 + b2.*p_shallow);

part2 = (B + b2.*p_deep).*(A + b2.*p_shallow);

part3 = (N.*b2 - M.*b1)./(b2.*(B - A));

enthalpy_diff = db2Pa.*(delta_p.*(a2 - 2*a3.*b1./b2 + 0.5*a3.*p_sum)./b2 + ...
                     (M./(2*b2)).*log(1 + delta_p.*(2*b1 + b2.*p_sum)./part1) + ... 
                     part3.*log(1 + delta_p.*b2.*(B - A)./part2));

%--------------------------------------------------------------------------
% This function calculates enthalpy_diff using the computationally
% efficient 48-term expression for density in terms of SA, CT and p.  If 
% one wanted to compute the enthalpy difference using the full TEOS-10 
% Gibbs function, the following lines of code will enable this.
%
%    pt = gsw_pt_from_CT(SA,CT);
%    pr0 = zeros(size(SA)); 
%    t_shallow = gsw_pt_from_t(SA,pt,pr0,p_shallow);
%    t_deep = gsw_pt_from_t(SA,pt,pr0,p_deep);
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
