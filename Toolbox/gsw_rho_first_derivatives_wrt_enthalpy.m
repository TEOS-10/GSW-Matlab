function [rho_SA, rho_h] = gsw_rho_first_derivatives_wrt_enthalpy(SA,CT,p)

% gsw_rho_first_derivatives_wrt_enthalpy                  first derivatives
%                                  specific volume with respect to enthalpy
% =========================================================================
%
% USAGE:
%  [rho_SA, rho_h] = gsw_rho_first_derivatives_wrt_enthalpy(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following two first-order derivatives of specific
%  volume (v),
%   (1) rho_SA, first-order derivative with respect to Absolute Salinity 
%       at constant CT & p.
%   (2) rho_h, first-order derivative with respect to SA & CT at 
%       constant p. 
%
%  Note that this function uses the using the computationally-efficient
%  expression for specific volume (Roquet et al., 2015).  There is an 
%  alternative to calling this function, namely 
%  gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(SA,CT,p) which uses 
%  the full Gibbs function (IOC et al., 2010).   
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
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_SA =  The first derivative of rho with respect to 
%              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
%  rho_h  =  The first derivative of rho with respect to 
%              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
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
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
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
   error('gsw_rho_first_derivatives_wrt_enthalpy:  Requires three inputs')
end %if

if ~(nargout == 2)
   error('gsw_rho_first_derivatives_wrt_enthalpy:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_rho_first_derivatives_wrt_enthalpy: SA and CT do not have the same dimensions')
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
    error('gsw_rho_first_derivatives_wrt_enthalpy: The dimensions of p do not agree')
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

rec_v = 1./gsw_specvol(SA,CT,p);
[v_SA, v_h] = gsw_specvol_first_derivatives_wrt_enthalpy(SA,CT,p);

rec_v2 = rec_v.^2;

rho_h = -v_h.*rec_v2;

rho_SA = -v_SA.*rec_v2;

if transposed
    rho_SA = rho_SA.';
    rho_h = rho_h.';
end

end
