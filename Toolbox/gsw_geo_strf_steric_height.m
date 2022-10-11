function steric_height = gsw_geo_strf_steric_height(SA,CT,p,p_ref)

% gsw_geo_strf_steric_height                          steric height anomaly
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:  
%  steric_height = gsw_geo_strf_steric_height(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates steric height anomaly as the pressure integral of specific 
%  volume anomaly from the pressure p of the “bottle” to the reference 
%  pressure p_ref, divided by the constant value of the gravitational 
%  acceleration, 9.7963 m s^-2.  That is, this function returns the dynamic
%  height anomaly divided by 9.7963 m s^-2; this being  the gravitational 
%  acceleration averaged over the surface of the global ocean (see page 46 
%  of Griffies, 2004).  Hence, steric_height is the steric height anomaly 
%  with respect to a given reference pressure p_ref.  
%
%  Dynamic height anomaly is the geostrophic streamfunction for the 
%  difference between the horizontal velocity at the pressure concerned, p, 
%  and the horizontal velocity at p_ref.  Dynamic height anomaly is the 
%  exact geostrophic streamfunction in isobaric surfaces even though the 
%  gravitational acceleration varies with latitude and pressure.  Steric 
%  height anomaly, being simply proportional to dynamic height anomaly, is 
%  also an exact geostrophic streamfunction in an isobaric surface (up to 
%  the constant of proportionality, 9.7963 m s^-2). 
%
%  Note however that steric_height is not exactly the height (in metres)
%  of an isobaric surface above a geopotential surface.  It is tempting to 
%  divide dynamic height anomaly by the local value of the gravitational
%  acceleration, but doing so robs the resulting quantity of either being
%     (i)  an exact geostrophic streamfunction, or 
%     (ii) exactly the height of an isobaric surface above a geopotential 
%          surface.
%  By using a constant value of the gravitational acceleration, we have
%  retained the first of these two properties.  So it should be noted that 
%  becasue of the variation of the gravitational acceleration with
%  latitude, steric_height does not exactly represent the height of an
%  isobaric surface above a geopotential surface under the assumption of
%  geostropy.  
%
%  The reference values used for the specific volume anomaly are 
%  SSO = 35.16504 g/kg and CT = 0 deg C.  This function calculates 
%  specific volume anomaly using the computationally efficient 75-term 
%  expression for specific volume of Roquet et al. (2015). Note that the 
%  75-term equation has been fitted in a restricted range of parameter 
%  space, and is most accurate inside the "oceanographic funnel" described
%  in McDougall et al. (2003).  For dynamical oceanography we may take the
%  75-term rational function expression for density as essentially
%  reflecting the full accuracy of TEOS-10.  The GSW internal library 
%  function "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to
%  test if some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref =  reference pressure                                     [ dbar ]
%           ( i.e. reference absolute pressure - 10.1325 dbar )
%
%  SA & CT  need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%  p_ref needs to be a single value, it can have dimensions 1x1 or Mx1 or  
%  1xN or MxN.
%
% OUTPUT:
%  steric_height  =  dynamic height anomaly divided by 9.7963 m s^-2  [ m ]
%   Note. If p_ref exceeds the pressure of the deepest “bottle” on a 
%     vertical profile, the steric height anomaly for each “bottle” 
%     on the whole vertical profile is returned as NaN.
%
% AUTHOR:  
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Griffies, S. M., 2004: Fundamentals of Ocean Climate Models. Princeton, 
%   NJ: Princeton University Press, 518 pp + xxxiv.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) and section 3.27 of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall and P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, 29-43.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_geo_strf_steric_height:  Requires four inputs')
end %if

unique_p_ref = unique(p_ref);
if ~isscalar(unique_p_ref)
    error('gsw_geo_strf_steric_height: The reference pressure p_ref must be unique')
end
clear p_ref
p_ref = unique_p_ref;

if p_ref < 0
    error('gsw_geo_strf_steric_height: The reference pressure p_ref must be positive')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_steric_height: SA & CT need to have the same dimensions')
elseif (ms*ns == 1)
    error('gsw_geo_strf_steric_height: There must be at least 2 values')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    error('gsw_geo_strf_steric_height: need more than one pressure')
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
    error('gsw_geo_strf_steric_height: Inputs array dimensions arguments do not agree')
end %if

[Inan] = find(isnan(SA + CT + p));
SA(Inan) = NaN;
CT(Inan) = NaN;
p(Inan) = NaN;

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

if max(p(:)) < p_ref
    error('gsw_geo_strf_steric_height: The reference pressure p_ref is deeper than all bottles')
end

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

dynamic_height_anomaly = gsw_geo_strf_dyn_height(SA,CT,p,p_ref);
const_grav = 9.7963;             % (Griffies, 2004);
steric_height = dynamic_height_anomaly./const_grav;  

if transposed
   steric_height = steric_height.';
end 

end
