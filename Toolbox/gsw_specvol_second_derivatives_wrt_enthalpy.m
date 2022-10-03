function [v_SA_SA_wrt_h, v_SA_h, v_h_h] = gsw_specvol_second_derivatives_wrt_enthalpy(SA,CT,p)

% gsw_specvol_second_derivatives_wrt_enthalpy      second order derivatives
%                               of volume specific with respect to enthalpy
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  [v_SA_SA_wrt_h, v_SA_h, v_h_h] = ...
%                   gsw_specvol_second_derivatives_wrt_enthalpy(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three first-order derivatives of specific
%  volume (v) with respect to enthalpy,
%   (1) v_SA_SA_wrt_h, second-order derivative with respect to Absolute Salinity 
%       at constant h & p.
%   (2) v_SA_h, second-order derivative with respect to SA & h at 
%       constant p. 
%   (3) v_h_h, second-order derivative with respect to h at 
%       constant SA & p. 
%
%  Note that this function uses the using the computationally-efficient
%  75 term expression for specific volume (Roquet et al., 2015).  There is
%  an alternative to calling this function, namely 
%  gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p) which uses 
%  the full Gibbs function (IOC et al., 2010).   
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described 
%  in McDougall et al. (2010).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  v_SA_SA_wrt_h = The second-order derivative of specific volume with
%                  respect to Absolute Salinity at constant h & p.
%                                                    [ (m^3/kg)(g/kg)^-2  ]
%  v_SA_h  = The second-order derivative of specific volume with respect to 
%            SA and h at constant p.        [ (m^3/kg)(g/kg)^-1 (J/kg)^-1 ]
%  v_h_h   = The second-order derivative with respect to h at 
%            constant SA & p.                         [ (m^3/kg)(J/kg)^-2 ]
%
% AUTHOR:   
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.15 (26th May, 2022)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
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
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_specvol_second_derivatives_wrt_enthalpy:  Requires three inputs')
end

if ~(nargout == 3)
   error('gsw_specvol_second_derivatives_wrt_enthalpy:  Requires three outputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_specvol_second_derivatives_wrt_enthalpy: SA and CT do not have the same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
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
    error('gsw_specvol_second_derivatives_wrt_enthalpy: The dimensions of p do not agree')
end

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

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

[dummy, v_CT, dummy] = gsw_specvol_first_derivatives(SA,CT,p);
[h_SA, h_CT] = gsw_enthalpy_first_derivatives(SA,CT,p);
[vCT_SA_SA, vCT_SA_CT, vCT_CT_CT, dummy, dummy] = gsw_specvol_second_derivatives(SA,CT,p);
[h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives(SA,CT,p);

rec_h_CT = 1./h_CT;
rec_h_CT2 = rec_h_CT.^2;

v_h_h = (vCT_CT_CT.*h_CT - h_CT_CT.*v_CT).*(rec_h_CT2.*rec_h_CT);

v_SA_h = (vCT_SA_CT.*h_CT - v_CT.*h_SA_CT).*rec_h_CT2 - h_SA.*v_h_h;
          
v_SA_SA_wrt_h = vCT_SA_SA - (h_CT.*(vCT_SA_CT.*h_SA + v_CT.*h_SA_SA) ...
              - v_CT.*h_SA.*h_SA_CT).*rec_h_CT2 - h_SA.*v_SA_h;

if transposed
    v_SA_SA_wrt_h = v_SA_SA_wrt_h.';
    v_SA_h = v_SA_h.';
    v_h_h = v_h_h.';
end

end
