function [v_SA_SA, v_SA_CT, v_CT_CT, v_SA_P, v_CT_P] = gsw_specvol_second_derivatives_CT_exact(SA,CT,p)

% gsw_specvol_second_derivatives_CT_exact             second derivatives of 
%                                                           specific volume
% =========================================================================
%
% USAGE:
%  [v_SA_SA, v_SA_CT, v_CT_CT, v_SA_P, v_CT_P] = ...
%                          gsw_specvol_second_derivatives_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of specific
%  volume (v),
%   (1) v_SA_SA, second-order derivative with respect to Absolute Salinity 
%       at constant CT & p.
%   (2) v_SA_CT, second-order derivative with respect to SA & CT at 
%       constant p. 
%   (3) v_CT_CT, second-order derivative with respect to CT at constant SA 
%       and p. 
%   (4) v_SA_P, second-order derivative with respect to SA & P at 
%       constant CT. 
%   (5) v_CT_P, second-order derivative with respect to CT & P at 
%       constant SA. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_specvol_second_derivatives(SA,CT,p) which uses the computationally
%  efficient 75 term expression for specific volume in terms of SA, CT  
%  and p (Roquet et al., 2015).   
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
%  v_SA_SA  =  The second derivative of specific volume with respect to 
%              Absolute Salinity at constant CT & p.  [ (m^3/kg)(g/kg)^-2 ]
%  v_SA_CT  =  The second derivative of specific volume with respect to 
%              SA and CT at constant p.           [ (m^3/kg)(g/kg)^-1 K^-1]
%  v_CT_CT  =  The second derivative of specific volume with respect to 
%              CT at constant SA and p.                  [ (m^3/kg) K^-2) ]
%  v_SA_P  =  The second derivative of specific volume with respect to 
%              SA and P at constant CT.                  [ (m^3/kg) Pa^-1 ]
%  v_CT_P  =  The second derivative of specific volume with respect to 
%              CT and P at constant SA.             [ (m^3/kg) K^-1 Pa^-1 ]
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
   error('gsw_specvol_second_derivatives_CT_exact:  Requires three inputs')
end %if

if ~(nargout == 5)
   error('gsw_specvol_second_derivatives_CT_exact:  Requires five outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_specvol_second_derivatives_CT_exact: SA and CT do not have the same dimensions')
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
    error('gsw_specvol_second_derivatives_CT_exact: The dimensions of p do not agree')
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

cp0 = gsw_cp0;             % from Eqn. 3.3.3 of IOC et al. (2010).
rec_cp0 = 1./cp0;
pr0 = zeros(size(SA)); 

pt0 = gsw_pt_from_CT(SA,CT);
rec_abs_pt0 = 1./(gsw_T0 + pt0);
cp0_div_abs_pt0 = cp0.*rec_abs_pt0;

t = gsw_pt_from_t(SA,pt0,pr0,p);

gamma = -gsw_gibbs(0,1,1,SA,t,p)./gsw_gibbs(0,2,0,SA,t,p);

rec_gTT = 1./gsw_gibbs(0,2,0,SA,t,p);
rec_gTT2 = rec_gTT.*rec_gTT;
rec_gTT_pt0 = 1./gsw_gibbs(0,2,0,SA,pt0,pr0);

gTTP = gsw_gibbs(0,2,1,SA,t,p);
gTTT = gsw_gibbs(0,3,0,SA,t,p);
gSAT = gsw_gibbs(1,1,0,SA,t,p);
gSATP = gsw_gibbs(1,1,1,SA,t,p);
gSATT = gsw_gibbs(1,2,0,SA,t,p);
gSAT_pt0 = gsw_gibbs(1,1,0,SA,pt0,pr0);
gSA_pt0 = gsw_gibbs(1,0,0,SA,pt0,pr0);
gSASA_pt0 = gsw_gibbs(2,0,0,SA,pt0,pr0);
gSASAP = gsw_gibbs(2,0,1,SA,t,p);
gSASAT = gsw_gibbs(2,1,0,SA,t,p);
gSAPP = gsw_gibbs(1,0,2,SA,t,p);
gTPP = gsw_gibbs(0,1,2,SA,t,p);

part_a = (gTTP + gamma.*gTTT).*rec_gTT2;
part_b = (gSATP + gamma.*gSATT).*rec_gTT;
part_c = (gTPP + gamma.*(2.*gTTP + gamma.*gTTT)).*rec_gTT;

v_CT_CT = cp0_div_abs_pt0.^2.*(gamma.*rec_abs_pt0.*rec_gTT_pt0 + part_a);

v_SA_CT = cp0_div_abs_pt0.*( ...
    gamma.*rec_abs_pt0.*gSAT_pt0.*rec_gTT_pt0 ...
    - part_b + gSAT.*part_a) - gSA_pt0.*v_CT_CT.*rec_cp0;

v_SA_SA = gSASAP + gamma.*gSASAT ...
    - gamma.*rec_abs_pt0.*gSASA_pt0 ...
    + gamma.*rec_abs_pt0.*(gSAT_pt0.^2).*rec_gTT_pt0 ...
    - 2.*gSAT.*part_b ...
    + (gSAT.^2).*part_a ...
    - 2.*gSA_pt0.*rec_cp0.*v_SA_CT - (gSA_pt0.^2).*(rec_cp0.^2).*v_CT_CT;

v_CT_P = -cp0_div_abs_pt0.*part_c;

v_SA_P = gSAPP + gamma.*(2.*gSATP + gamma.*gSATT) ...
    + part_c.*(gSA_pt0.*rec_abs_pt0 - gSAT);

if transposed
    v_SA_SA = v_SA_SA.';
    v_SA_CT = v_SA_CT.';
    v_CT_CT = v_CT_CT.';
    v_SA_P = v_SA_P.';
    v_CT_P = v_CT_P.';
end

end
