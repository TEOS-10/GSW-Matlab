function [N2, N2_p, N2_specvol, N2_alpha, N2_beta, dSA, dCT, dp] = gsw_Nsquared_min(SA,CT,p,lat)

% gsw_Nsquared_min                         minimum buoyancy (Brunt-Vaisala)
%                               frequency squared (N^2)  (75-term equation)
%==========================================================================
% 
% USAGE:  
%  [N2, N2_p, N2_specvol, N2_alpha, N2_beta, dSA, dCT, dp] = ...
%                                           gsw_Nsquared_min(SA,CT,p,{lat})
%
% DESCRIPTION:
%  Calculates the minimum buoyancy frequency squared (N^2)(i.e. the 
%  Brunt-Vaisala frequency squared) between two bottles from the equation,
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
%  N2         =  minimum Brunt-Vaisala Frequency squared      [ rad^2/s^2 ]
%  N2_p       =  pressure of minimum N2                            [ dbar ]
%  N2_specvol =  specific volume at the minimum N2                [ kg/m3 ]
%  N2_alpha   =  thermal expansion coefficient with respect         [ 1/K ]
%                to Conservative Temperature at the minimum N2
%  N2_beta    =  saline contraction coefficient at constant        [ kg/g ]
%                Conservative Temperature at the minimum N2
%  dSA        =  difference in salinity between bottles            [ g/kg ]
%  dCT        =  difference in Conservative Temperature between   [ deg C ]
%                bottles
%  dp         =  difference in pressure between bottles            [ dbar ]
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
%  Griffies, S. M., 2004: Fundamentals of Ocean Climate Models. Princeton, 
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
   error('gsw_Nsquared_min:  Requires three or four inputs')
end 
if ~(nargout >= 2)
   error('gsw_Nsquared_min:  Requires at least two outputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_Nsquared_min: SA and CT must have same dimensions')
end

if (ms*ns == 1)
    error('gsw_Nsquared_min: There must be at least 2 bottles')
end

if (mp == 1) & (np == 1)
    error('gsw_Nsquared_min:  There must be at least 2 bottles')
elseif (ns == np) & (mp == 1)
    p = p(ones(1,ms), :);
elseif (ms == mp) & (np == 1)
    p = p(:,ones(1,ns));
elseif (ns == mp) & (np == 1)
    p = p.'; 
    p = p(ones(1,ms), :);
elseif (ms == np) & (mp == 1)
     p = p.';  
     p = p(:,ones(1,ns));
elseif (ms == np) & (ns == mp)
     p = p.';   
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_Nsquared_min: Inputs array dimensions arguments do not agree')
end 

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

[mp,number_profiles] = size(p); 

if exist('lat','var')
    if transposed
        lat = lat.';
    end
    [mL,nL] = size(lat);
    if (mL == 1) & (nL == 1)
        lat = lat*ones(mp,number_profiles);
    elseif (number_profiles == nL) & (mL == 1)
        lat = lat(ones(1,mp), :);
    elseif (mp == mL) & (nL == 1)
        lat = lat(:,ones(1,number_profiles));
    elseif (number_profiles == mL) & (nL == 1)
        lat = lat.';
        lat = lat(ones(1,mp), :);
    elseif (mp == nL) & (mL == 1)
        lat = lat.';
        lat = lat(:,ones(1,nL));
    elseif (mp == nL) & (number_profiles == mL)
        lat = lat.';
    elseif (mp == mL) & (number_profiles == nL)
        % ok
    else
        error('gsw_Nsquared_min: Inputs array dimensions arguments do not agree')
    end
    grav = gsw_grav(lat,p);
else
    grav = 9.7963*ones(mp,number_profiles);             % (Griffies, 2004)
end 

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

db2Pa = 1e4;
Ishallow = 1:(mp-1);
Ideep = 2:mp;

dp = p(Ideep,:) - p(Ishallow,:);

dSA = SA(Ideep,:) - SA(Ishallow,:);
dCT = CT(Ideep,:) - CT(Ishallow,:);

[specvol_bottle, alpha_bottle, beta_bottle] = gsw_specvol_alpha_beta(SA,CT,p);

N2_shallow = (grav(Ishallow,:).*grav(Ishallow,:)./(specvol_bottle(Ishallow,:).*db2Pa.*dp)).*(beta_bottle(Ishallow,:).*dSA - alpha_bottle(Ishallow,:).*dCT);
N2_deep = (grav(Ideep,:).*grav(Ideep,:)./(specvol_bottle(Ideep,:).*db2Pa.*dp)).*(beta_bottle(Ideep,:).*dSA - alpha_bottle(Ideep,:).*dCT);

N2 = nan(mp-1,number_profiles);
N2_p = nan(mp-1,number_profiles);
N2_specvol = nan(mp-1,number_profiles);
N2_alpha = nan(mp-1,number_profiles);
N2_beta = nan(mp-1,number_profiles);

for Iprofile = 1:number_profiles
    dummy_N2 = [N2_shallow(:,Iprofile),N2_deep(:,Iprofile)];
    dummy_p = [p(Ishallow,Iprofile),p(Ideep,Iprofile)];
    dummy_specvol = [specvol_bottle(Ishallow,Iprofile),specvol_bottle(Ideep,Iprofile)];
    dummy_alpha = [alpha_bottle(Ishallow,Iprofile),alpha_bottle(Ideep,Iprofile)];
    dummy_beta = [beta_bottle(Ishallow,Iprofile),beta_bottle(Ideep,Iprofile)];

    [N2(:,Iprofile),IN2] = min(dummy_N2,[],2);
    
    for Ibottle = 1:mp-1
        N2_p(Ibottle,Iprofile) = dummy_p(Ibottle,IN2(Ibottle));
        N2_specvol(Ibottle,Iprofile) = dummy_specvol(Ibottle,IN2(Ibottle));
        N2_alpha(Ibottle,Iprofile) = dummy_alpha(Ibottle,IN2(Ibottle));
        N2_beta(Ibottle,Iprofile) = dummy_beta(Ibottle,IN2(Ibottle));
    end
end

if transposed
    N2 = N2.';
    N2_p = N2_p.';
    N2_specvol = N2_specvol.';
    N2_alpha = N2_alpha.';
    N2_beta = N2_beta.';
    dSA = dSA.';
    dCT = dCT.';
    dp = dp.';
end

end