function [v_SA, v_h] = gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(SA,CT,p)

% gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact     first derivatives
%                               of specific volume with respect to enthalpy
% =========================================================================
%
% USAGE:
%  [v_SA, v_h] = ...
%              gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following two first-order derivatives of specific
%  volume (v),
%   (1) v_SA, first-order derivative with respect to Absolute Salinity 
%       at constant h & p.
%   (2) v_h, first-order derivative with respect to h at constant SA & p. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_specvol_first_derivatives_wrt_enthalpy(SA,CT,p) which uses the 
%  computationally efficient 75 term expression for specific volume in 
%  terms of SA, CT and p (Roquet et al., 2015).   
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
%  v_SA =  The first derivative of specific volume with respect to 
%          Absolute Salinity at constant h & p.         [ J/(kg (g/kg)^2) ]
%  v_h  =  The first derivative of specific volume with respect to 
%          h at constant SA & p.                         [ J/(kg K(g/kg)) ]
%
% AUTHOR:   
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
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
   error('gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact:  Requires three inputs')
end %if

if ~(nargout == 2)
   error('gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact: SA and CT do not have the same dimensions')
end %if

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
    error('gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact: The dimensions of p do not agree')
end %if

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

[vCT_SA, vCT_CT, dummy] = gsw_specvol_first_derivatives_CT_exact(SA,CT,p);
[h_SA, h_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT,p);

rec_h_CT = 1./h_CT;

v_SA = vCT_SA - (vCT_CT.*h_SA).*rec_h_CT;

v_h = vCT_CT.*rec_h_CT;

if transposed
    v_SA = v_SA.';
    v_h = v_h.';
end

end
