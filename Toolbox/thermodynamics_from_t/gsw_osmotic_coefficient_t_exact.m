function osmotic_coefficient_t_exact = gsw_osmotic_coefficient_t_exact(SA,t,p)

% gsw_osmotic_coefficient_t_exact                       osmotic coefficient
%==========================================================================
%
% USAGE:
%  osmotic_coefficient_t_exact = gsw_osmotic_coefficient_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the osmotic coefficient of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  osmotic_coefficient_t_exact  =  osmotic coefficient of seawater    
%                                                              [ unitless ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_osmotic_coefficient_t_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_osmotic_coefficient_t_exact: SA and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % p is a scalar - fill to size of SA
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
    error('gsw_osmotic_coefficient_t_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    t = t.';
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

%R = 8.314472;                                         % molar gas constant


%M_S = 0.0314038218;  % mole-weighted average atomic weight of the elements 
                      % of Reference-Composition sea salt, in units of 
                      % kg mol^-1.  Strictly speaking, the formula below 
                      % applies only to seawater of Reference Composition. 
                      % If molality and the osmotic coefficienit is 
                      % required to an accuracy of better than 0.1% we 
                      % suggest you contact the authors for further 
                      % guidance.


k = 3.777007343340624e-3;                                      % k = M_S/R
                                         
part = k.*(1000 - SA)./(gsw_T0 + t); 

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).

x2 = sfac.*SA;
x = sqrt(x2);
y = t.*0.025;
z = p.*1e-4; %Note that the input pressure (p) is sea pressure in units of dbar.

oc01 =  7.231916621570606e1;
oc02 =  1.059039593127674e1;
oc03 = -3.025914794694813e1;
oc04 =  5.040733670521486e1;
oc05 = -4.074543321119333e1;
oc06 =  1.864215613820487e1;
oc07 = -3.022566485046178;
oc08 = -6.138647522851840;
oc09 =  1.353207379758663e1;
oc10 = -7.316560781114737;
oc11 =  1.829232499785750;
oc12 = -5.358042980767074e-1;
oc13 = -1.705887283375562;
oc14 = -1.246962174707332e-1;
oc15 =  1.228376913546017;
oc16 =  1.089364009088042e-2;
oc17 = -4.264828939262248e-1;
oc18 =  6.213127679460041e-2;
oc19 =  2.481543497315280;
oc20 = -1.363368964861909;
oc21 = -5.640491627443773e-1;
oc22=   1.344724779893754;
oc23 = -2.180866793244492;
oc24 =  4.765753255963401;
oc25 = -5.726993916772165;
oc26 =  2.918303792060746;
oc27 = -6.506082399183509e-1;
oc28 = -1.015695507663942e-1;
oc29 =  1.035024326471108;
oc30 = -6.742173543702397e-1;
oc31 =  8.465642650849419e-1;
oc32 = -7.508472135244717e-1;
oc33 = -3.668086444057845e-1;
oc34 =  3.189939162107803e-1;
oc35 = -4.245629194309487e-2;
 
tl = oc01 + oc02*y ...
    + x.*(oc03 + x.*(oc04 + x.*(oc05 + x.*(oc06 + oc07*x))) ...
    + y.*(oc08 + x.*(oc09 + x.*(oc10 + oc11*x))...
    + y.*(oc12 + oc13*x + y.*(oc14 + oc15*x + y.*(oc16 + x.*(oc17 + oc18*y))))) ...
    + z.*(oc19 + x.*(oc20 + oc21*y + oc22*x) + y.*(oc23 + y.*(oc24 + y.*(oc25 + oc26*y))) ...
    + z.*(oc27 + oc28*x + y.*(oc29 + oc30*y) ...
    + z.*(oc31 + oc32*x + y.*(oc33 + oc34*y) + oc35*z))));

osmotic_coefficient_t_exact = tl.*part;

if transposed
    osmotic_coefficient_t_exact = osmotic_coefficient_t_exact.';
end

end
