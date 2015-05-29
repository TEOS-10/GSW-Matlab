function enthalpy = gsw_enthalpy(SA,CT,p)

% gsw_enthalpy                                specific enthalpy of seawater
%                                                        (48-term equation)
%==========================================================================
%
% USAGE:
%  enthalpy = gsw_enthalpy(SA,CT,p)
%
% DESCRIPTION:
%  Calculates specific enthalpy of seawater using the computationally-
%  efficient 48-term expression for density in terms of SA, CT and p
%  (IOC et al., 2010)
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2010).  The GSW library function 
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
%  enthalpy  =  specific enthalpy                                  [ J/kg ]
%
% AUTHOR: 
%  Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.    
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.6) of this TEOS-10 Manual. 
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

if ~(nargin == 3)
    error('gsw_enthalpy: requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_enthalpy: SA and CT must have same dimensions')
end

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
    error('gsw_enthalpy_: Inputs array dimensions arguments do not agree')
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

db2Pa = 1e4;                      % factor to convert from dbar to Pa
cp0 = 3991.86795711963;           % from Eqn. (3.3.3) of IOC et al. (2010).

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

a0 = v21 + CT.*(v22 + CT.*(v23 + CT.*(v24 + v25*CT))) ...
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

part = (N.*b2 - M.*b1)./(b2.*(B - A));

enthalpy = cp0*CT + ...
           db2Pa.*(p.*(a2 - 2*a3.*b1./b2 + 0.5*a3.*p)./b2 + ...
           (M./(2*b2)).*log(1 + p.*(2*b1 + b2.*p)./b0) + ...
           part.*log(1 + (b2.*p.*(B - A))./(A.*(B + b2.*p))));

%--------------------------------------------------------------------------
% This function calculates enthalpy using the computationally-efficient
% 48-term expression for density in terms of SA, CT and p.  If one wanted to
% compute enthalpy from SA, CT, and p with the full TEOS-10 Gibbs function,  
% the following lines of code will enable this.
%
%    pt = gsw_pt_from_CT(SA,CT);
%    pr0 = zeros(size(SA)); 
%    t = gsw_pt_from_t(SA,pt,pr0,p);
%    enthalpy = gsw_enthalpy_t_exact(SA,t,p);
%
%    or call the following, it is identical to the lines above.
%
%    enthalpy = gsw_enthalpy_CT_exact(SA,CT,p)
%
%-----------------This is the end of the alternative code------------------

if transposed
    enthalpy = enthalpy.';
end

end
