function t_freezing = gsw_t_freezing(SA,p,saturation_fraction)

% gsw_t_freezing              in-situ temperature at which seawater freezes
%==========================================================================
%
% USAGE:
%  t_freezing = gsw_t_freezing(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the in-situ temperature at which seawater freezes. The 
%  in-situ temperature freezing point is calculated from the exact 
%  in-situ freezing temperature which is found by a modified Newton-Raphson
%  iteration (McDougall and Wotherspoon, 2013) of the equality of the 
%  chemical potentials of water in seawater and in ice.
%
%  An alternative GSW function, gsw_t_freezing_poly, it is based on a 
%  computationally-efficient polynomial, and is accurate to within -5e-4 K 
%  and 6e-4 K, when compared with this function.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in 
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default 
%    is 0, air free) 
%
%  p & saturation_fraction (if provided) may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
%               (ITS-90)                
%
% AUTHOR: 
%  Trevor McDougall, Paul Barker and Rainer Feistal    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.  
%
%  McDougall T.J., and S.J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3) 
   error('gsw_t_freezing: Requires either two or three inputs')
end %if

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end
    
if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_t_freezing: saturation fraction MUST be between zero and one.')
end
    
[ms,ns] = size(SA);
[mp,np] = size(p);
[msf,nsf] = size(saturation_fraction);

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_t_freezing: Inputs array dimensions arguments do not agree')
end %if

if (msf == 1) & (nsf == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(SA));         % fill to size of SA
elseif (ns == nsf) & (msf == 1)                        % saturation_fraction is row vector,
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (nsf == 1)                     % saturation_fraction is column vector,
    saturation_fraction = saturation_fraction(:,ones(1,ns));        % copy across each row.
elseif (ns == msf) & (nsf == 1)           % saturation_fraction is a transposed row vector,
    saturation_fraction = saturation_fraction.';                           % transposed then
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (ns == nsf)
    % ok
else
    error('gsw_t_freezing: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

%   The following code gives a rather accurate polynomial-based expression 
%   for the freezing temperature, adjsted for the saturation fraction.  
%   This is the value that is used as the seed for the modified Newton's
%   method. The error of the following polynomial ranges between -8e-4 K 
%   and 3e-4 K when compared with the outout of this function.  

c0 = 0.002519;

c1 = -5.946302841607319;
c2 =  4.136051661346983;
c3 = -1.115150523403847e1;
c4 =  1.476878746184548e1;
c5 = -1.088873263630961e1;
c6 =  2.961018839640730;
    
c7 = -7.433320943962606;
c8 = -1.561578562479883;
c9 =  4.073774363480365e-2;

c10 =  1.158414435887717e-2;
c11 = -4.122639292422863e-1;
c12 = -1.123186915628260e-1;
c13 =  5.715012685553502e-1;
c14 =  2.021682115652684e-1;
c15 =  4.140574258089767e-2;
c16 = -6.034228641903586e-1;
c17 = -1.205825928146808e-2;
c18 = -2.812172968619369e-1;
c19 =  1.877244474023750e-2;
c20 = -1.204395563789007e-1;
c21 =  2.349147739749606e-1;
c22 =  2.748444541144219e-3;

SA_r = SA.*1e-2;
x = sqrt(SA_r);
p_r = p.*1e-4;

%  The initial value of t_freezing_exact (for air-free seawater)
tf = c0 ...
     + SA_r.*(c1 + x.*(c2 + x.*(c3 + x.*(c4 + x.*(c5 + c6.*x))))) ...
     + p_r.*(c7 + p_r.*(c8 + c9.*p_r)) ...
     + SA_r.*p_r.*(c10 + p_r.*(c12 + p_r.*(c15 + c21.*SA_r)) + SA_r.*(c13 + c17.*p_r + c19.*SA_r) ...
     + x.*(c11 + p_r.*(c14 + c18.*p_r)  + SA_r.*(c16 + c20.*p_r + c22.*SA_r)));

df_dt = 1000.*gsw_t_deriv_chem_potential_water_t_exact(SA,tf,p) - gsw_gibbs_ice(1,0,tf,p);
%  df_dt here is the initial value of the derivative of the function f whose
%  zero (f = 0) we are finding (see Eqn. (3.33.2) of IOC et al (2010)).  

tf_old = tf;
f = 1000.*gsw_chem_potential_water_t_exact(SA,tf_old,p) - gsw_gibbs_ice(0,0,tf_old,p);
tf = tf_old - f./df_dt ; % this is half way through the modified method (McDougall and Wotherspoon, 2013)
tfm = 0.5.*(tf + tf_old);
df_dt = 1000.*gsw_t_deriv_chem_potential_water_t_exact(SA,tfm,p) - gsw_gibbs_ice(1,0,tfm,p);
tf = tf_old - f./df_dt; % this is the end of one iteration of the modified Newton method

tf_old = tf;
f = 1000.*gsw_chem_potential_water_t_exact(SA,tf_old,p) - gsw_gibbs_ice(0,0,tf_old,p);
tf = tf_old - f./df_dt ; % this is half way through the modified method (McDougall and Wotherspoon, 2013)

% Adjust for the effects of dissolved air
t_freezing = tf - saturation_fraction.*(1e-3).*(2.4 - SA./70.33008); 

%  If the data is outside the range of applicability of TEOS-10, set the
%  output to NaN.  
t_freezing(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if transposed
    t_freezing = t_freezing.';
end

% The maximum error is 2x10^-13 degrees C for this code which has one and a
% half iterations of the modified Newton's method.  This is machine precision
% for this calculation.  This is the maximum error over the whole 
% (SA,p) domain with SA varying between 0 and 42 g/kg and p varying 
% independently between 0 and 10,000 dbar. 

end