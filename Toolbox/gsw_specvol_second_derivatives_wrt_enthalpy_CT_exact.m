function [v_SA_SA, v_SA_h, v_h_h] = gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p)

% gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact   second derivatives
%                               of specific volume with respect to enthalpy
% =========================================================================
%
% USAGE:
%  [v_SA_SA, v_SA_h, v_h_h] = ...
%             gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of specific
%  volume (v) with respect to enthalpy,
%   (1) v_SA_SA, second-order derivative with respect to Absolute Salinity 
%       at constant h & p.
%   (2) v_SA_h, second-order derivative with respect to SA & h at 
%       constant p. 
%   (3) v_h_h, second-order derivative with respect to h at 
%       constant SA & p. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_specvol_first_derivatives_wrt_enthalpy(SA,CT,p) which uses the 
%  computationally efficient expression for specific volume in terms of SA,
%  CT and p (Roquet et al., 2015).   
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
%  v_SA_SA = The second-order derivative of specific volume with respect to 
%            Absolute Salinity at constant h & p.       [ J/(kg (g/kg)^2) ]
%  v_SA_h  = The second-order derivative of specific volume with respect to 
%            SA and h at constant p.                     [ J/(kg K(g/kg)) ]
%  v_h_h   = The second-order derivative with respect to h at 
%            constant SA & p.
%
% AUTHOR:   
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact:  Requires three inputs')
end %if

if ~(nargout == 3)
   error('gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact: SA and CT do not have the same dimensions')
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
    error('gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact: The dimensions of p do not agree')
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

[dummy, v_CT, dummy] = gsw_specvol_first_derivatives_CT_exact(SA,CT,p);
[h_SA, h_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT,p);
[vCT_SA_SA, vCT_SA_CT, vCT_CT_CT,dummy,dummy] = gsw_specvol_second_derivatives_CT_exact(SA,CT,p);
[h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives_CT_exact(SA,CT,p);

rec_h_CT = 1./h_CT;
rec_h_CT2 = rec_h_CT.^2;

v_h_h = (vCT_CT_CT.*h_CT - h_CT_CT.*v_CT).*(rec_h_CT2.*rec_h_CT);

v_SA_h = (vCT_SA_CT.*h_CT - v_CT.*h_SA_CT).*rec_h_CT2 - h_SA.*v_h_h;

v_SA_SA = vCT_SA_SA - (h_CT.*(vCT_SA_CT.*h_SA - v_CT.*h_SA_SA) ...
              + v_CT.*h_SA.*h_SA_CT).*rec_h_CT2 - h_SA.*v_SA_h;

if transposed
    v_SA_SA = v_SA_SA.';
    v_SA_h = v_SA_h.';
    v_h_h = v_h_h.';
end

end
