function [N2, p_mid] = gsw_Nsquared(SA,CT,p,lat)

% gsw_Nsquared             buoyancy (Brunt-Vaisala) frequency squared (N^2)
%                                                        (48-term equation)
%==========================================================================
% 
% USAGE:  
%  [N2, p_mid] = gsw_Nsquared(SA,CT,p,{lat})
%
% DESCRIPTION:
%  Calculates the buoyancy frequency squared (N^2)(i.e. the Brunt-Vaisala 
%  frequency squared) at the mid pressure from the equation,
%
%           2      2     d(rho_local)
%         N   =  g   x  --------------
%                           dP
%
%  Note. This routine uses rho from "gsw_rho", which is the computationally
%    efficient 48-term expression for density in terms of SA, CT and p.      
%  Note also that the pressure increment, dP, in the above formula is in 
%    Pa, so that it is 10^4 times the pressure increment dp in dbar. 
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2011).  The GSW library function 
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
%  N2     =  Brunt-Vaisala Frequency squared  (M-1xN)             [ 1/s^2 ]
%  p_mid  =  Mid pressure between p grid      (M-1xN)              [ dbar ]
%
% AUTHOR:  
%  Trevor McDougall and Paul Barker                   [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 3.0 (22nd March, 2011)
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
%  McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2011:  A 
%   computationally efficient 48-term expression for the density of 
%   seawater in terms of Conservative Temperature, and related properties
%   of seawater.  To be submitted to Ocean Science Discussions. 
%
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3 | nargin == 4)
   error('gsw_Nsquared:  Requires three or four inputs')
end %if
if ~(nargout == 2)
   error('gsw_Nsquared:  Requires two outputs')
end %if

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
end %if

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
    if transposed
        lat = lat.';
    end
    [mL,nL] = size(lat);
    [ms,ns] = size(SA);
    if (mL == 1) & (nL == 1)              % lat scalar - fill to size of SA
        lat = lat*ones(size(SA));
    elseif (ns == nL) & (mL == 1)         % lat is row vector,
        lat = lat(ones(1,ms), :);          % copy down each column.
    elseif (ms == mL) & (nL == 1)         % lat is column vector,
        lat = lat(:,ones(1,ns));           % copy across each row.
    elseif (ms == mL) & (ns == nL)
        % ok
    else
        error('gsw_Nsquared: Inputs array dimensions arguments do not agree')
    end %if
    grav = gsw_grav(lat,p);
else
    grav = 9.7963*ones(size(p));             % (Griffies, 2004)
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

db2Pa = 1e4;
Ishallow = 1:(mp-1);
Ideep = 2:mp;
p_mid = (p(Ishallow,:) + p(Ideep,:))/2;

d_rho_local_deep = gsw_rho(SA(Ideep,:),CT(Ideep,:),p_mid);
d_rho_local_shallow = gsw_rho(SA(Ishallow,:),CT(Ishallow,:),p_mid);
d_rho_local = d_rho_local_deep - d_rho_local_shallow;

%--------------------------------------------------------------------------
% This function calculates d_rho_local using the computationally-efficient 
% 48-term expression for density in terms of SA, CT and p. If one wanted to
% compute d_rho_local with the full TEOS-10 Gibbs function expression for 
% density, the following lines of code will enable this.
%
%    pt = gsw_pt_from_CT(SA,CT);
%    pr0 = zeros(size(SA)); 
%    t = gsw_pt_from_t(SA,pt,pr0,p);
%    d_rho_local = gsw_rho_t_exact(SA(Ideep,:),t(Ideep,:),p_mid) - ...
%                    gsw_rho_t_exact(SA(Ishallow,:),t(Ishallow,:),p_mid);
%
%----This is the end of the alternative code to evaluate d_rho_local-------

grav_local = (grav(Ishallow,:) + grav(Ideep,:))/2;
d_p = (p(Ideep,:) - p(Ishallow,:) );

% the pressure difference is converted from dbar to Pa 
N2 = (grav_local.*grav_local).*(d_rho_local)./(db2Pa.*d_p);

if transposed
    N2 = N2.';
    p_mid = p_mid.';
end

end