function SA = gsw_SA_from_rho_CT_exact(rho,CT,p)

% gsw_SA_from_rho_CT_exact                   Absolute Salinity from density 
% =========================================================================
%
% USAGE:
%  SA = gsw_SA_from_rho_CT_exact(rho,CT,p)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of a seawater sample, for given values
%  of its density, Conservative Temperature and sea pressure (in dbar). 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_SA_from_rho(rho,CT,p), which uses the computationally 
%  efficient 75-term expression for density in terms of SA, CT and p 
%  (Roquet et al., 2015).   
%
% INPUT:
%  rho =  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it. 
%     That is, it is 'density', NOT 'density anomaly'.
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
%  Trevor McDougall & Paul Barker                     [ help_gsw@csiro.au ]
%      
% VERSION NUMBER: 3.06.12 (25th May, 2020)
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

if ~(nargin==3)
   error('gsw_SA_from_rho_CT_exact:  Requires three inputs')
end %if

[md,nd] = size(rho);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= md | nt ~= nd)
    error('gsw_SA_from_rho_CT_exact: rho and CT must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of rho
    p = p*ones(size(rho));
elseif (nd == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,nd));               % copy across each row.
elseif (nd == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,md), :);                          % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
    error('gsw_SA_from_rho_CT_exact: Inputs array dimensions arguments do not agree')
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
v_0 = gsw_specvol_CT_exact(zeros(size(rho)),CT,p);
v_120 = gsw_specvol_CT_exact(120*ones(size(rho)),CT,p);
 
SA = 120*(v_lab - v_0)./(v_120 - v_0);            % initial estimate of SA.
SA(SA < 0 | SA > 120) = NaN;

v_SA = (v_120 - v_0)./120; %initial estimate of v_SA, the SA derivative of v

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure 
%--------------------------------------------------------------------------

for Number_of_iterations = 1:2 
    SA_old = SA;
    delta_v = gsw_specvol_CT_exact(SA_old,CT,p) - v_lab;
    SA = SA_old - delta_v./v_SA ; % this is half way through the modified N-R method (McDougall and Wotherspoon, 2013)
    SA_mean = 0.5*(SA + SA_old);
    [v_SA, dummy, dummy] = gsw_specvol_first_derivatives_CT_exact(SA_mean,CT,p);
    SA = SA_old - delta_v./v_SA;
    SA(SA < 0 | SA > 120) = NaN; 
end

% After two iterations of this modified Newton-Raphson iteration,
% the error in SA is no larger than 8x10^-13 g kg^-1, which 
% is machine precision for this calculation. 
 
if transposed
    SA = SA.';
end

end
