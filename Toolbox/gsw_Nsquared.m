function [N2, p_mid] = gsw_Nsquared(SA,CT,p,lat)

% gsw_Nsquared             buoyancy (Brunt-Vaisala) frequency squared (N^2)
%                                                        (75-term equation)
%==========================================================================
% 
% USAGE:  
%  [N2, p_mid] = gsw_Nsquared(SA,CT,p,{lat})
%
% DESCRIPTION:
%  Calculates the buoyancy frequency squared (N^2)(i.e. the Brunt-Vaisala 
%  frequency squared) at the mid pressure from the equation,
%
%           2      2     beta.dSA - alpha.dCT
%         N   =  g  . -------------------------
%                         specvol_local.dP
%
%  The pressure increment, dP, in the above formula is in Pa, so that it is
%  10^4 times the pressure increment dp in dbar. 
%
%  Note. This routine uses specvol from "gsw_specvol", which is based on 
%  the computationally efficient expression for specific volume in terms of 
%  SA, CT and p (Roquet et al., 2015).
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
% OPTIONAL:
%  lat  =  latitude in decimal degrees north                [ -90 ... +90 ]
%  Note. If lat is not supplied, a default gravitational acceleration
%    of 9.7963 m/s^2 (Griffies, 2004) will be applied.
%
%  SA & CT need to have the same dimensions. 
%  p & lat (if provided) may have dimensions 1x1 or Mx1 or 1xN or MxN, 
%  where SA & CT are MxN.
%
% OUTPUT:
%  N2     =  Brunt-Vaisala Frequency squared  (M-1xN)         [ rad^2/s^2 ]
%  p_mid  =  Mid pressure between p grid      (M-1xN)              [ dbar ]
%
%  The units of N2 are radians^2 s^-2 however in may textbooks this is
%  abreviated to s^-2 as radians does not have a unit.  To convert the
%  frequency to hertz, cycles sec^-1, divide the frequency by pi, 
%  i.e. N/(2 x pi).
%
% AUTHOR:  
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (4th August, 2021)
%
% REFERENCES:
%  Griffies, S.M., 2004: Fundamentals of Ocean Climate Models. Princeton, 
%   NJ: Princeton University Press, 518 pp + xxxiv.
%   
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.10 and Eqn. (3.10.2) of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3 | nargin == 4)
   error('gsw_Nsquared:  Requires three or four inputs')
end 

if ~(nargout == 2)
   error('gsw_Nsquared:  Requires two outputs')
end 

p = gsw_resize(p,size(SA));

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_Nsquared: SA and CT must have same dimensions')
end

if (ms*ns == 1)
    error('gsw_Nsquared: There must be at least 2 bottles')
end

if (mp == 1) & (np == 1)              % p is a scalar - must be two bottles
    error('gsw_Nsquared:  There must be at least 2 bottles')
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
    error('gsw_Nsquared: Inputs array dimensions arguments do not agree')
end 

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

[mp,np] = size(p);

if exist('lat','var')
%    lat = gsw_resize(lat,size(SA));
    grav = gsw_grav(lat,p);
else
    grav = 9.7963*ones(size(p));             % (Griffies, 2004)
end 

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

db2Pa = 1e4;
Ishallow = 1:(mp-1);
Ideep = 2:mp;

grav_local = 0.5*(grav(Ishallow,:,:) + grav(Ideep,:,:));

dSA = (SA(Ideep,:,:) - SA(Ishallow,:,:));
SA_mid = 0.5*(SA(Ishallow,:,:) + SA(Ideep,:,:));
dCT = (CT(Ideep,:,:) - CT(Ishallow,:,:));
CT_mid = 0.5*(CT(Ishallow,:,:) + CT(Ideep,:,:));
dp = (p(Ideep,:,:) - p(Ishallow,:,:));
p_mid = 0.5*(p(Ishallow,:,:) + p(Ideep,:,:));

[specvol_mid, alpha_mid, beta_mid] = gsw_specvol_alpha_beta(SA_mid,CT_mid,p_mid);

%--------------------------------------------------------------------------
% This function calculates rho, alpha & beta using the computationally
% efficient expression for specific volume in terms of SA, CT and p.  If 
% one wanted to use the full TEOS-10 Gibbs function expression for specific
% volume, the following lines of code will enable this.
%
%    specvol_mid = gsw_specvol_CT_exact(SA_mid,CT_mid,p_mid);
%    alpha_mid = gsw_alpha_CT_exact(SA_mid,CT_mid,p_mid);
%    beta_mid = gsw_beta_CT_exact(SA_mid,CT_mid,p_mid);
%
%--This is the end of the alternative code to evaluate rho, alpha & beta---

N2 = ((grav_local.*grav_local)./(specvol_mid.*db2Pa.*dp)).*(beta_mid.*dSA - alpha_mid.*dCT);

if transposed
    N2 = N2.';
    p_mid = p_mid.';
end

end