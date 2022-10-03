function [rho_SA_SA, rho_SA_CT, rho_CT_CT, rho_SA_P, rho_CT_P] = gsw_rho_second_derivatives_CT_exact(SA,CT,p)

% gsw_rho_second_derivatives_CT_exact              SA and CT partial second 
%                                              order derivatives of density
%==========================================================================
% 
% USAGE:  
% [rho_SA_SA, rho_SA_CT, rho_CT_CT, rho_SA_P, rho_CT_P] = ...
%                              gsw_rho_second_derivatives_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of rho, 
%   (1) rho_SA_SA, second-order derivative with respect to Absolute  
%       Salinity at constant CT & p.
%   (2) rho_SA_CT, second-order derivative with respect to SA & CT at 
%       constant p. 
%   (3) rho_CT_CT, second-order derivative with respect to CT at 
%       constant SA & p. 
%   (4) rho_SA_P, second-order derivative with respect to SA & P at 
%       constant CT. 
%   (5) rho_CT_P, second-order derivative with respect to CT & P at 
%       constant SA. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_rho_second_derivatives(SA,CT,p), which uses the computationally
%  efficient polynomial for specific volume in terms of SA, CT and p
%  (Roquet et al., 2015).    
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
%  rho_SA_SA = The second-order derivative of rho with respect to 
%              Absolute Salinity at constant CT & p.  [ (kg/m^3)(g/kg)^-2 ]
%  rho_SA_CT = The second-order derivative of rho with respect to 
%              SA and CT at constant p.          [ (kg/m^3)(g/kg)^-1 K^-1 ]
%  rho_CT_CT = The second-order derivative of rho with respect to CT at 
%              constant SA & p                            [ (kg/m^3) K^-2 ]
%  rho_SA_P  = The second-order derivative with respect to SA & P at 
%              constant CT.                    [ (kg/m^3) (g/kg)^-1 Pa^-1 ]
%  rho_CT_P  = The second-order derivative with respect to CT & P at 
%              constant SA.                         [ (kg/m^3) K^-1 Pa^-1 ]
% 
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual. 
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_rho_second_derivatives_CT_exact:  Requires three inputs')
end %if
if ~(nargout == 5)
   error('gsw_rho_second_derivatives_CT_exact:  Requires five outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_rho_second_derivatives_CT_exact: SA and CT must have same dimensions')
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
    error('gsw_rho_second_derivatives_CT_exact: Inputs array dimensions arguments do not agree')
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
[v_SA, v_CT, v_P] = gsw_specvol_first_derivatives_CT_exact(SA,CT,p);
[v_SA_SA, v_SA_CT, v_CT_CT, v_SA_P, v_CT_P] = gsw_specvol_second_derivatives_CT_exact(SA,CT,p);

rec_v2 = rec_v.^2;
rec_v3 = rec_v2.*rec_v;

rho_CT_CT = -v_CT_CT.*rec_v2 + 2.*v_CT.^2.*rec_v3;

rho_SA_CT = -v_SA_CT.*rec_v2 + 2.*v_SA.*v_CT.*rec_v3;

rho_SA_SA = -v_SA_SA.*rec_v2 + 2.*v_SA.^2.*rec_v3;

rho_SA_P = -v_SA_P.*rec_v2 + 2.*v_SA.*v_P.*rec_v3;

rho_CT_P = -v_CT_P.*rec_v2 + 2.*v_CT.*v_P.*rec_v3;


if transposed
    rho_SA_SA = rho_SA_SA.'; 
    rho_SA_CT = rho_SA_CT.';
    rho_CT_CT = rho_CT_CT.';
    rho_SA_P = rho_SA_P.';
    rho_CT_P = rho_CT_P.';
end

end
