function [u_SA_SA, u_SA_CT, u_CT_CT, u_SA_P, u_CT_P] = gsw_internal_energy_second_derivatives(SA,CT,p)

% gsw_internal_energy_second_derivatives              second derivatives of
%                                       specific interal energy of seawater
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  [u_SA_SA, u_SA_CT, u_CT_CT, u_SA_P, u_CT_P] = ...
%                           gsw_internal_energy_second_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of 
%  internal energy,
%  (1) u_SA_SA, second order derivative with respect to Absolute Salinity
%      at constant CT & p.
%  (2) u_SA_CT, second order derivative with respect to SA & CT at
%      constant p.
%  (3) u_CT_CT, second order derivative with respect to CT at constant
%      SA & p.
%  (4) u_SA_P, second-order derivative with respect to SA & P at 
%      constant CT. 
%  (5) u_CT_P, second-order derivative with respect to CT & P at 
%      constant SA.    
%
%  Note that this function uses the using the computationally-efficient
%  75-term expression for specific volume (Roquet et al., 2015).  There is  
%  an alternative to calling this function, namely 
%  gsw_internal_energy_second_derivatives_CT_exact(SA,CT,p) which uses the
%  full Gibbs function (IOC et al., 2010).   
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  u_SA_SA =  The second derivative of internal energy with respect to
%             Absolute Salinity at constant CT & p.     [ (J/kg)(g/kg)^-2 ]
%  u_SA_CT =  The second derivative of internal energy with respect to
%             SA & CT at constant p.                [ (J/kg)(g/kg)^-1 K^-1]
%  u_CT_CT =  The second derivative of internal energy with respect to
%             CT at constant SA and p.                      [ (J/kg) K^-2 ]
%  u_SA_P  =  The second derivative of internal energy with respect to
%             SA & P at constant CT.              [ (J/kg)(g/kg)^-1 Pa^-1 ]
%  u_CT_P  =  The second derivative of internal energy with respect to
%             CT & P at constant SA.                  [ (J/kg) K^-1 Pa^-1 ]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
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
    error('gsw_internal_energy_second_derivatives: requires three inputs')
end
if ~(nargout == 5)
    error('gsw_internal_energy_second_derivatives: requires five outputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_internal_energy_second_derivatives: SA and CT must have same dimensions')
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
    error('gsw_internal_energy_second_derivatives: Inputs array dimensions arguments do not agree')
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

db2Pa = 1e4;               %  dbar to Pa conversion factor 
P = (db2Pa.*p + gsw_P0);

[h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives(SA,CT,p);
[v_SA_SA, v_SA_CT, v_CT_CT, v_SA_P, v_CT_P] = gsw_specvol_second_derivatives(SA,CT,p);

u_SA_SA = h_SA_SA - P.*v_SA_SA;

u_SA_CT = h_SA_CT - P.*v_SA_CT;

u_CT_CT = h_CT_CT - P.*v_CT_CT;

u_SA_P = -P.*v_SA_P;

u_CT_P = -P.*v_CT_P;

if transposed
    u_SA_SA = u_SA_SA.';
    u_SA_CT = u_SA_CT.';
    u_CT_CT = u_CT_CT.';
    u_SA_P = u_SA_P.';
    u_CT_P = u_CT_P.';
end

end
