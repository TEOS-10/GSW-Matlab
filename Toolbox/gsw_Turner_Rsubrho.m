function [Tu, Rsubrho, p_mid] = gsw_Turner_Rsubrho(SA,CT,p)

% gsw_Turner_Rsubrho              Turner angle & Rsubrho (75-term equation)
%==========================================================================
% 
% USAGE:  
%  [Tu, Rsubrho, p_mid] = gsw_Turner_Rsubrho(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the Turner angle and the Rsubrho as a function of pressure 
%  down a vertical water column.  These quantities express the relative 
%  contributions of the vertical gradients of Conservative Temperature 
%  and Absolute Salinity to the vertical stability (the square of the 
%  Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at 
%  the mid pressure between the individual data points in the vertical.  
%  This function uses computationally-efficient 75-term expression for 
%  specific volume in terms of SA, CT and p (Roquet et al., 2015).  Note
%  that in the double-diffusive literature, papers concerned with the 
%  "diffusive" form of double-diffusive convection often define the 
%  stability ratio as the reciprocal of what is defined here as the 
%  stability ratio. 
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:  
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions, 
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  Tu       =  Turner angle, on the same (M-1)xN grid as p_mid.
%              Turner angle has units of:           [ degrees of rotation ]
%  Rsubrho  =  Stability Ratio, on the same (M-1)xN grid as p_mid.
%              Rsubrho is dimensionless.                       [ unitless ]
%  p_mid    =  mid pressure between the indivual points of the p grid. 
%              That is, p_mid is on a (M-1)xN grid in the vertical.  
%              p_mid has units of:                                 [ dbar ]
%
% AUTHOR:  
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (30th July, 2021)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.15.1) and (3.16.1) of this TEOS-10 Manual. 
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
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_Turner_Rsubrho:  Requires three inputs')
end 

if ~(nargout == 3 )
   error('gsw_Turner_Rsubrho:  Requires three outputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns )
    error('gsw_Turner_Rsubrho: SA and CT must have same dimensions')
end

if (ms*ns == 1)
    error('There must be at least 2 values')
end

if (mp == 1) & (np == 1)              % p is a scalar - must be two bottles
    error('There must be at least 2 pressure values')
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (np == ms) & (np == 1)          % p is a transposed row vector,
    p = p.';                                         % transposed then
    p = p(ones(1,ms), :);                    % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_Turner_Rsubrho: Inputs array dimensions arguments do not agree')
end 

[mp,np] = size(p);

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

Ishallow = 1:(mp-1);
Ideep = 2:mp;
p_mid = 0.5*(p(Ishallow,:) + p(Ideep,:));
SA_mid = 0.5*(SA(Ishallow,:) + SA(Ideep,:));
CT_mid = 0.5*(CT(Ishallow,:) + CT(Ideep,:));

dSA = SA(Ishallow,:) - SA(Ideep,:);
dCT = CT(Ishallow,:) - CT(Ideep,:);

[dummy, alpha, beta] = gsw_specvol_alpha_beta(SA_mid,CT_mid,p_mid);

%--------------------------------------------------------------------------
% This function evaluates Tu and Rsubrho using the computationally-efficient
% 75-term expression for specific volume in terms of SA, CT and p. If one 
% wanted to compute Tu and Rsubrho using the full TEOS-10 Gibbs function 
% expression, the following lines of code would do that.  
%
%    t_mid = gsw_t_from_CT(SA_mid,CT_mid,p_mid);
%    beta = gsw_beta_const_CT_t_exact(SA_mid,t_mid,p_mid);
%    alpha = gsw_alpha_wrt_CT_t_exact(SA_mid,t_mid,p_mid);
%
% --------------This is the end of the alternative code--------------------

Tu = atan2((alpha.*dCT + beta.*dSA),(alpha.*dCT - beta.*dSA));
Tu = Tu.*(180/pi);

Rsubrho = nan(size(dSA));
Rsubrho(dSA ~= 0) = (alpha(dSA ~= 0).*dCT(dSA ~= 0))./(beta(dSA ~= 0).*dSA(dSA ~= 0));

end
