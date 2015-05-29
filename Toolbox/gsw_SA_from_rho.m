function SA = gsw_SA_from_rho(rho,CT,p)

% gsw_SA_from_rho                            Absolute Salinity from density
% =========================================================================
%
% USAGE:
%  SA = gsw_SA_from_rho(rho,CT,p)
% 
% DESCRIPTION:
%  Calculates the Absolute Salinity of a seawater sample, for given values
%  of its density, Conservative Temperature and sea pressure (in dbar). 
%  This function uses the computationally-efficient 48-term expression for 
%  density in terms of SA, CT and p (IOC et al., 2013).
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2013).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  rho =  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it. 
%     That is, it is 'density', not 'density anomaly'.
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & CT are MxN.
%
% OUTPUT:
%  SA  =  Absolute Salinity.                                       [ g/kg ]
%   Note. This is expressed on the Reference-Composition Salinity
%     Scale of Millero et al. (2008). 
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.5 of this TEOS-10 Manual. 
%
%  Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: 
%   The composition of Standard Seawater and the definition of the 
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
   error('gsw_SA_from_rho:  Requires three inputs')
end %if

[md,nd] = size(rho);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= md | nt ~= nd)
    error('gsw_SA_from_rho: rho and CT must have same dimensions')
end

if (mp == 1) & (np == 1)               % p scalar - fill to size of rho
    p = p*ones(size(rho));
elseif (nd == np) & (mp == 1)          % p is row vector,
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (np == 1)          % p is column vector,
    p = p(:,ones(1,nd));               % copy across each row.
elseif (nd == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                            % transposed then
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
    error('gsw_SA_from_rho: Inputs array dimensions arguments do not agree')
end %if

if md == 1
    rho = rho.';
    CT = CT.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

v_lab = ones(size(rho))./rho;
v_0 = gsw_specvol(zeros(size(rho)),CT,p);
v_50 = gsw_specvol(50*ones(size(rho)),CT,p);
 
SA = 50*(v_lab - v_0)./(v_50 - v_0);            % initial estimate of SA.
SA(SA < 0 | SA > 50) = NaN;

v_SA = (v_50 - v_0)./50; %initial estimate of v_SA, the SA derivative of v

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure 
%--------------------------------------------------------------------------

for Number_of_iterations = 1:2 
    SA_old = SA;
    delta_v = gsw_specvol(SA_old,CT,p) - v_lab;
    SA = SA_old - delta_v./v_SA ; % this is half way through the modified N-R method (McDougall and Wotherspoon, 2012)
    SA_mean = 0.5*(SA + SA_old);
    [rho,alpha,beta] = gsw_rho_alpha_beta(SA_mean,CT,p);
    v_SA = - beta./rho; 
    SA = SA_old - delta_v./v_SA;
    SA(SA < 0 | SA > 50) = NaN; 
end

% After two iterations of this modified Newton-Raphson iteration,
% the error in SA is no larger than 8x10^-13 g/kg, which 
% is machine precision for this calculation. 
 
if transposed
    SA = SA.';
end

end
