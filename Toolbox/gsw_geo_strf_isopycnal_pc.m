function [geo_strf_isopycnal_pc, p_mid] = gsw_geo_strf_isopycnal_pc(SA,CT,delta_p,gamma_n,layer_index,A)

% gsw_geo_strf_isopycnal_pc                     McDougall-Klocker piecewise 
%                    constant geostrophic streamfunction (75-term equation)
%==========================================================================
%
% USAGE:
%  [geo_strf_isopycnal_pc, p_mid] = gsw_geo_strf_isopycnal_pc(SA,CT,delta_p,gamma_n,layer_index,A)
%
% DESCRIPTION:
%  Calculates the McDougall-Klocker geostrophic streamfunction (see Eqn.
%  (3.30.1) of IOC et al. (2010).  This function is to used when the 
%  Absolute Salinity and Conservative Temperature are piecewise constant in
%  the vertical over sucessive pressure intervals of delta_p (such as in a 
%  forward "z-coordinate" ocean model, and in isopycnal layered ocean 
%  models).  The McDougall-Klocker geostrophic streamfunction is designed 
%  to be used as the geostrophic streamfunction in an approximately neutral
%  surface (such as a Neutral Density surface, a potential density surface 
%  or an omega surface (Klocker et al. (2009)).  Reference values of 
%  Absolute Salinity, Conservative Temperature and pressure are found by 
%  interpolation of a one-dimensional look-up table, with the interpolating
%  variable being Neutral Density (gamma_n) or sigma_2.  This function 
%  calculates specific volume anomaly using the computationally efficient 
%  75-term expression for specific volume of Roquet et al. (2015).
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  CT       =  Conservative Temperature (ITS-90)                  [ deg C ]
%  delta_p  =  difference in sea pressure between the deep and shallow
%              extents of each layer in which SA and CT are vertically 
%              constant.  delta_p must be positive.                [ dbar ]
%  Note. Sea pressure is absolute pressure minus 10.1325 dbar.
%
%  gamma_n     = Neutral Density anomaly                         [ kg/m^3 ]
%                ( i.e. Neutral Density minus 1000 kg/m^3 )
%  layer_index = Index of the layers of the gamma_n surfaces
%  A           = if nothing is entered the programme defaults to "Neutral
%                Density" as the vertical interpolating variable. 
%              = 's2' or 'sigma2', for sigma_2 as the vertical 
%                interpolating variable. 
%
%  SA, CT & delta_p need to have the same dimensions.
%  gamma_n & layer_indx need to have the same dimensions, there should be
%  only one "gamma_n" or "sigma_2" value per level of interest.
%  A needs to be 1x1.
%
% OUTPUT:
%  geo_strf_isopycnal_pc =  McDougall & Klocker (2010)          [ m^2/s^2 ]
%                           geostrophic streamfunction          
%  p_mid                 =  mid-point pressure in each layer       [ dbar ]
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
%    See section 3.30 of this TEOS-10 Manual.
%
%  Jackett, D. R. and T. J. McDougall, 1997: A neutral density variable
%   for the world’s oceans. Journal of Physical Oceanography, 27, 237-263.
%
%  Klocker, A., T. J. McDougall and D. R. Jackett, 2009: A new method 
%   for forming approximately neutral surfaces.  Ocean Sci., 5, 155-172. 
%
%  McDougall, T. J. and A. Klocker, 2010: An approximate geostrophic
%   streamfunction for use in density surfaces.  Ocean Modelling, 32,
%   105-117.  
%    The McDougall-Klocker geostrophic streamfunction is defined in
%    Eqn. (62) of this paper.
%    See section 8 of this paper for a discussion of this piecewise-
%    constant version of the McDougall-Klocker geostrophic streamfunction. 
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

if ~(nargin == 5 | nargin == 6)
    error('gsw_geo_strf_isopycnal_pc:  Requires five or six inputs')
end %if

if ~(nargout == 2)
    error('gsw_geo_strf_isopycnal_pc:  Requires two outputs')
end %if

if ~exist('A','var')
    A = 'gn';
elseif ~ischar(A)
    A = 'gn';
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(delta_p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_isopycnal_pc: SA & CT need to have the same dimensions')
end

if (mp == 1) & (np == 1)              % delta_p scalar - fill to size of SA
    error('gsw_geo_strf_isopycnal_pc: needs more than one pressure');
elseif (ns == np) & (mp == 1)         % delta_p is row vector,
    delta_p = delta_p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % delta_p is column vector,
    delta_p = delta_p(:,ones(1,ns));               % copy across each row.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_geo_strf_isopycnal_pc: Inputs array dimensions arguments do not agree')
end %if

[mgn,ngn] = size(gamma_n);
[mli,nli] = size(layer_index);

if (mgn~=mli) | (ngn~=nli)
    error('gsw_geo_strf_isopycnal_pc: Inputs array (layers) dimensions arguments do not agree')
end

if ms == 1  % row vector
    delta_p = delta_p(:);
    CT = CT(:);
    SA = SA(:);
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

db2Pa = 1e4;

SA_iref_cast = NaN(size(gamma_n));
CT_iref_cast = SA_iref_cast;
p_iref_cast = SA_iref_cast;

[Inn] = find(~isnan(gamma_n));
[SA_iref_cast(Inn),CT_iref_cast(Inn),p_iref_cast(Inn)] = gsw_interp_ref_cast(gamma_n(Inn),A);

[dyn_height_pc, p_mid_fine] = gsw_geo_strf_dyn_height_pc(SA,CT,delta_p);

p_mid_fine(p_mid_fine == 0) = NaN;
p_mid = p_mid_fine(layer_index,:);

SA_iref_cast_nd = SA_iref_cast.* ones(size(p_mid));
CT_iref_cast_nd = CT_iref_cast.* ones(size(p_mid));

part1 = 0.5*db2Pa*(p_mid_fine(layer_index,:) - p_iref_cast).*(gsw_specvol(SA(layer_index,:),CT(layer_index,:),p_mid_fine(layer_index,:)) - ...
                                   gsw_specvol(SA_iref_cast_nd,CT_iref_cast_nd,p_mid_fine(layer_index,:)));

part2 = -0.225e-15*db2Pa*db2Pa*(CT(layer_index,:) - CT_iref_cast).*...
    (p_mid_fine(layer_index,:) - p_iref_cast).*(p_mid_fine(layer_index,:) - p_iref_cast);

part3 = dyn_height_pc(layer_index,:) - gsw_enthalpy_SSO_0(p_mid_fine(layer_index,:)) + ...
        gsw_enthalpy(SA_iref_cast_nd,CT_iref_cast_nd,p_mid_fine(layer_index,:)) - gsw_cp0*CT_iref_cast;
    
%--------------------------------------------------------------------------
% This function calculates the McDougall-Klocker piecewise constant 
% streamfunction using the computationally efficient 75-term expression for
% specific volume in terms of SA, CT and p. If one wanted to compute this
% with the full TEOS-10 Gibbs function expression for specific volume, the 
% following lines of code will enable this. Note that dynamic height will
% also need to be evaluated using the full Gibbs function.
%
%    part1 = 0.5*db2Pa*(p_mid_fine(layer_index,:) - p_iref_cast).*(gsw_specvol_CT_exact(SA(layer_index,:),CT(layer_index,:),p_mid_fine(layer_index,:)) - ...
%                                    gsw_specvol_CT_exact(SA_iref_cast,CT_iref_cast,p_mid_fine(layer_index,:))); 
%    part2 = -0.225e-15*db2Pa*db2Pa*(CT(layer_index,:) - CT_iref_cast).*...
%     (p_mid_fine(layer_index,:) - p_iref_cast).*(p_mid_fine(layer_index,:) - p_iref_cast);
%    SSO = gsw_SSO*ones(size(SA));
%    CT_0 = zeros(size(CT));
%    part3 = dyn_height_pc(layer_index,:) - gsw_enthalpy_CT_exact(SSO,CT_0,p_mid_fine(layer_index,:)) + ...
%         gsw_enthalpy_CT_exact(SA_iref_cast_nd,CT_iref_cast_nd,p_mid_fine(layer_index,:)) - gsw_cp0*CT_iref_cast;
%
%---------------This is the end of the alternative code--------------------
 
geo_strf_isopycnal_pc = part1 + part2 + part3;

if transposed
    geo_strf_isopycnal_pc = geo_strf_isopycnal_pc.';
    p_mid = p_mid.';
end %if

end
