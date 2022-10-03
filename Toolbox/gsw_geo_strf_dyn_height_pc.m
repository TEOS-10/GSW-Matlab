function [geo_strf_dyn_height_pc, p_mid] = gsw_geo_strf_dyn_height_pc(SA,CT,delta_p)

% gsw_geo_strf_dyn_height_pc                     dynamic height anomaly for
%                            piecewise constant profiles (75-term equation)
%==========================================================================
%
% USAGE:  
%  [geo_strf_dyn_height_pc, p_mid] = gsw_geo_strf_dyn_height_pc(SA,CT,delta_p)
%
% DESCRIPTION:
%  Calculates dynamic height anomaly as the integral of specific volume 
%  anomaly from the the sea surface pressure (0 Pa) to the pressure p.
%  This function, gsw_geo_strf_dyn_height_pc, is to used when the 
%  Absolute Salinity and Conservative Temperature are piecewise constant in 
%  the vertical over sucessive pressure intervals of delta_p (such as in
%  a forward "z-coordinate" ocean model).  "geo_strf_dyn_height_pc" is
%  the dynamic height anomaly with respect to the sea surface.  That is, 
%  "geo_strf_dyn_height_pc" is the geostrophic streamfunction for the 
%  difference between the horizontal velocity at the pressure concerned, p,
%  and the horizontal velocity at the sea surface.  Dynamic height anomaly 
%  is the geostrophic streamfunction in an isobaric surface.  The reference
%  values used for the specific volume anomaly are SA = SSO = 35.16504 g/kg
%  and CT = 0 deg C.  The output values of geo_strf_dyn_height_pc are 
%  given at the mid-point pressures, p_mid, of each layer in which SA and 
%  CT are vertically piecewice constant (pc).  This function calculates 
%  enthalpy using the computationally-efficient 75-term expression for 
%  specific volume of Roquet et al., (2015). 
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  CT       =  Conservative Temperature (ITS-90)                  [ deg C ]
%  delta_p  =  difference in sea pressure between the deep and     [ dbar ]
%              shallow extents of each layer in which SA and CT
%              are vertically constant. delta_p must be positive.
%              
%  Note. sea pressure is absolute pressure minus 10.1325 dbar.
%
%  SA & CT need to have the same dimensions.
%  delta_p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  geo_strf_dyn_height_pc =  dynamic height anomaly             [ m^2/s^2 ]
%  p_mid                  =  mid-point pressure in each layer      [ dbar ]
%
% AUTHOR: 
%  Trevor McDougall and Claire Roberts-Thomson         [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) and (A.30.6) of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  McDougall, T.J. and A. Klocker, 2010: An approximate geostrophic
%   streamfunction for use in density surfaces. Ocean Modelling, 32,
%   105-117.  
%    See section 8 of this paper.  
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

if ~(nargin == 3)
   error('gsw_geo_strf_dyn_height_pc:  Requires three inputs')
end 

if ~(nargout == 2)
   error('gsw_geo_strf_dyn_height_pc:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mdp,ndp] = size(delta_p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_dyn_height_pc: SA & CT need to have the same dimensions')
elseif (ms*ns == 1)
    error('gsw_geo_strf_dyn_height_pc: There must be at least 2 values')
end

if (mdp == 1) & (ndp == 1)              % delta_p scalar - fill to size of SA
    error('gsw_geo_strf_dyn_height_pc: need more than one pressure')
elseif (ns == ndp) & (mdp == 1)         % delta_p is row vector,
    delta_p = delta_p(ones(1,ms), :);   % copy down each column.
elseif (ms == mdp) & (ndp == 1)          % delta_p is column vector,
    delta_p = delta_p(:,ones(1,ns));    % copy across each row.
elseif (ms == mdp) & (ns == ndp)
    % ok
else
    error('gsw_geo_strf_dyn_height_pc: Inputs array dimensions arguments do not agree')
end %if

transposed = 0;
if ms == 1 
   delta_p  =  delta_p(:);
   CT  =  CT(:);
   SA  =  SA(:);
   transposed = 1;
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

if any(delta_p < 0)
    error('gsw_geo_strf_dyn_height_pc: pressure must be monotonic')
end

p_deep = cumsum(delta_p); 
p_shallow = p_deep - delta_p;

delta_h = gsw_enthalpy_diff(SA,CT,p_shallow,p_deep);
dyn_height_deep = -cumsum(delta_h);
%            This is Phi minus Phi_0 of Eqn. (3.32.2) of IOC et al. (2010).

p_mid = 0.5*(p_shallow  + p_deep);
delta_h_half = gsw_enthalpy_diff(SA,CT,p_mid,p_deep);

geo_strf_dyn_height_pc = gsw_enthalpy_SSO_0(p_mid) + ...
                           dyn_height_deep + delta_h_half;

%--------------------------------------------------------------------------
% This function calculates dynamic height anomaly piecewise constant using 
% the computationally-efficient 75-term expression for specific volume in
% terms of SA, CT and p. If one wanted to compute dynamic height anomaly 
% with the full TEOS-10 Gibbs function expression for specific volume, the 
% following lines of code will enable this.
%
%   delta_h = gsw_enthalpy_diff_CT_exact(SA,CT,p_shallow,p_deep);
%   dyn_height_deep = -cumsum(delta_h);
%   p_mid = 0.5*(p_shallow  + p_deep);
%   delta_h_half = gsw_enthalpy_diff_CT_exact(SA,CT,p_mid,p_deep);
%   SSO = gsw_SSO*ones(size(SA));
%   CT_0 = zeros(size(SSO));
%   geo_strf_dyn_height_pc = gsw_enthalpy_CT_exact(SSO,CT_0,p_mid) + ...
%                               dyn_height_deep + delta_h_half;
%
%---------------This is the end of the alternative code--------------------

if transposed
    geo_strf_dyn_height_pc = geo_strf_dyn_height_pc.';
    p_mid = p_mid.';
end

end
