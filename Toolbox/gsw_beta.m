function beta = gsw_beta(SA,CT,p)

% gsw_beta                       saline contraction coefficient at constant
%                               Conservative Temperature (48-term equation)
%==========================================================================
%
% USAGE:  
%  beta = gsw_beta(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the saline (i.e. haline) contraction coefficient of seawater  
%  at constant Conservative Temperature using the computationally-efficient
%  48-term expression for density in terms of SA, CT and p 
%  (McDougall et al., 2010).
%   
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2010).  The GSW library function 
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
%  beta  =  saline contraction coefficient                         [ kg/g ]
%           at constant Conservative Temperature
%    
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2012)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.19.3) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_beta:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_beta: SA and CT must have same dimensions')
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
    error('gsw_beta: Inputs array dimensions arguments do not agree')
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

b01 = -6.698001071123802;
b02 = -2.986498947203215e-2;
b03 =  2.327859407479162e-4;
b04 = -5.983233568452735e-2;
b05 =  7.643133860820750e-4;
b06 = -2.140477007450431e-5;
b07 =  2.467559060524383e-7;
b08 = -1.806789763745328e-4;
b09 =  6.876837219536232e-7;
b10 =  1.550932729220080e-10;
b11 = -7.521448093615448e-3;
b12 = -2.764306979894411e-5;
b13 =  1.262937315098546e-7;
b14 =  9.527875081696435e-10;
b15 = -1.811147201949891e-11;
b16 = -4.954963307079632e-5;
b17 =  5.702346883314446e-7;
b18 = -1.150931530388857e-8;
b19 = -6.951273511674217e-11;
b20 =  4.021645853353715e-12;
b21 =  1.083865310229748e-5;
b22 = -1.105097577149576e-7;
b23 =  6.211426728363857e-10;
b24 =  1.119522344879478e-14;

sqrtSA = sqrt(SA);

v_hat_denominator = v01 + CT.*(v02 + CT.*(v03 + v04*CT))  ...
             + SA.*(v05 + CT.*(v06 + v07*CT) ...
         + sqrtSA.*(v08 + CT.*(v09 + CT.*(v10 + v11*CT)))) ...
              + p.*(v12 + CT.*(v13 + v14*CT) + SA.*(v15 + v16*CT) ...
              + p.*(v17 + CT.*(v18 + v19*CT) + v20*SA));
          
v_hat_numerator = v21 + CT.*(v22 + CT.*(v23 + CT.*(v24 + v25*CT))) ...
           + SA.*(v26 + CT.*(v27 + CT.*(v28 + CT.*(v29 + v30*CT))) + v36*SA ...
       + sqrtSA.*(v31 + CT.*(v32 + CT.*(v33 + CT.*(v34 + v35*CT)))))  ...
            + p.*(v37 + CT.*(v38 + CT.*(v39 + v40*CT))  ...
           + SA.*(v41 + v42*CT) ...
            + p.*(v43 + CT.*(v44 + v45*CT + v46*SA) ...
            + p.*(v47 + v48*CT)));
       
dvhatden_dSA = b01 + CT.*(b02 + b03*CT) ...
    + sqrtSA.*(b04 + CT.*(b05 + CT.*(b06 + b07*CT))) ...
         + p.*(b08 + b09*CT + b10*p) ;

dvhatnum_dSA = b11 + CT.*(b12 + CT.*(b13 + CT.*(b14 + b15*CT))) ...
    + sqrtSA.*(b16 + CT.*(b17 + CT.*(b18 + CT.*(b19 + b20*CT)))) + b21*SA ...
         + p.*(b22 + CT.*(b23 + b24*p));

beta = (v_hat_numerator.*dvhatden_dSA - v_hat_denominator.*dvhatnum_dSA)./ ...
           (v_hat_numerator.*v_hat_denominator);

if transposed
    beta = beta.';
end

end
