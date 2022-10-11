function [v_SA, v_CT, v_P] = gsw_specvol_first_derivatives_CT_exact(SA,CT,p)

% gsw_specvol_first_derivatives_CT_exact            first order derivatives 
%                                                        of specific volume
% =========================================================================
%
% USAGE:
%  [v_SA, v_CT, v_P] = gsw_specvol_first_derivatives_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three first-order derivatives of specific
%  volume (v),
%   (1) v_SA, first-order derivative with respect to Absolute Salinity 
%       at constant CT & p.
%   (2) v_CT, first-order derivative with respect to SA & CT at 
%       constant p. 
%   (3) v_P, first-order derivative with respect to CT at constant SA 
%       and p. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_specvol_first_derivatives(SA,CT,p) which uses the computationally
%  efficient 75 term expression for density in terms of SA, CT and p 
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
%  v_SA  =  The first derivative of specific volume with respect to 
%           Absolute Salinity at constant CT & p.     [ (m^3/kg)(g/kg)^-1 ]
%  v_CT  =  The first derivative of specific volume with respect to 
%           CT at constant SA and p.                         [ m^3/(K kg) ]
%  v_P   =  The first derivative of specific volume with respect to 
%           P at constant SA and CT.                        [ m^3/(Pa kg) ]
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
   error('gsw_specvol_first_derivatives_CT_exact:  Requires three inputs')
end %if

if ~(nargout == 3)
   error('gsw_specvol_first_derivatives_CT_exact:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_specvol_first_derivatives_CT_exact: SA and CT do not have the same dimensions')
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
    error('gsw_specvol_first_derivatives_CT_exact: The dimensions of p do not agree')
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

pr0 = zeros(size(SA)); 

pt0 = gsw_pt_from_CT(SA,CT);
rec_abs_pt0 = 1./(gsw_T0 + pt0);
t = gsw_pt_from_t(SA,pt0,pr0,p);

rec_gTT = 1./gsw_gibbs(0,2,0,SA,t,p);
gSAP = gsw_gibbs(1,0,1,SA,t,p);
gTP = gsw_gibbs(0,1,1,SA,t,p);
gSAT = gsw_gibbs(1,1,0,SA,t,p);
gSA_pt0 = gsw_gibbs(1,0,0,SA,pt0,pr0);
gPP = gsw_gibbs(0,0,2,SA,t,p);

v_CT = -gsw_cp0.*gTP.*rec_abs_pt0.*rec_gTT;

v_SA = gSAP - gTP.*(gSAT - rec_abs_pt0.*gSA_pt0).*rec_gTT;

v_P = gPP - gTP.*gTP.*rec_gTT;

if transposed
    v_SA = v_SA.';
    v_CT = v_CT.';
    v_P = v_P.';
end

end
