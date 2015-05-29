function [rho_SA_SA, rho_SA_h, rho_h_h] = gsw_rho_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p)

% gsw_rho_second_derivatives_wrt_enthalpy_CT_exact       second derivatives
%                                           of rho with respect to enthalpy
% =========================================================================
%
% USAGE:
%  [rho_SA_SA, rho_SA_h, rho_h_h] = ...
%             gsw_rho_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of rho with 
%  respect to enthalpy,
%   (1) rho_SA_SA, second-order derivative with respect to Absolute Salinity 
%       at constant h & p.
%   (2) rho_SA_h, second-order derivative with respect to SA & h at 
%       constant p. 
%   (3) rho_h_h, second-order derivative with respect to h at 
%       constant SA & p. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_rho_second_derivatives_wrt_enthalpy(SA,CT,p) which uses the 
%  computationally efficient 75-term expression for specific volume in 
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
%  rho_SA_SA =  The second-order derivative of rho with respect to 
%               Absolute Salinity at constant h & p.    [ J/(kg (g/kg)^2) ]
%  rho_SA_h  =  The second-order derivative of rho with respect to 
%               SA and h at constant p.                  [ J/(kg K(g/kg)) ]
%  rho_h_h   =  The second-order derivative of rho with respect to h at 
%               constant SA & p.
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
   error('gsw_rho_second_derivatives_wrt_enthalpy_CT_exact:  Requires three inputs')
end %if

if ~(nargout == 3)
   error('gsw_rho_second_derivatives_wrt_enthalpy_CT_exact:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_rho_second_derivatives_wrt_enthalpy_CT_exact: SA and CT do not have the same dimensions')
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
    error('gsw_rho_second_derivatives_wrt_enthalpy_CT_exact: The dimensions of p do not agree')
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

rec_v = 1./gsw_specvol_CT_exact(SA,CT,p);
[v_SA, v_h] = gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(SA,CT,p);
[v_SA_SA, v_SA_h, v_h_h] = gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p);

rec_v2 = rec_v.^2;
rec_v3 = rec_v2.*rec_v;

rho_h_h = -v_h_h.*rec_v2 + 2.*v_h.^2.*rec_v3;

rho_SA_h = -v_SA_h.*rec_v2 + 2.*v_SA.*v_h.*rec_v3;

rho_SA_SA = -v_SA_SA.*rec_v2 + 2.*v_SA.^2.*rec_v3;

if transposed
    rho_SA_SA = rho_SA_SA.';
    rho_SA_h = rho_SA_h.';
    rho_h_h = rho_h_h.';
end

end
