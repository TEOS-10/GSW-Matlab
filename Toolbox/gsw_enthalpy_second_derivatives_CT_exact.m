function [h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives_CT_exact(SA,CT,p)

% gsw_enthalpy_second_derivatives_CT_exact   second derivatives of enthalpy
% =========================================================================
%
% USAGE:
%  [h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of specific
%  enthalpy (h),
%   (1) h_SA_SA, second-order derivative with respect to Absolute Salinity 
%       at constant CT & p.
%   (2) h_SA_CT, second-order derivative with respect to SA & CT at 
%       constant p. 
%   (3) h_CT_CT, second-order derivative with respect to CT at constant SA 
%       and p. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_enthalpy_second_derivatives(SA,CT,p) which uses the computationally
%  efficient 75-term expression for specific volume in terms of SA, CT and 
%  p (Roquet et al., 2015).
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
%  h_SA_SA  =  The second derivative of specific enthalpy with respect to 
%              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
%  h_SA_CT  =  The second derivative of specific enthalpy with respect to 
%              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
%  h_CT_CT  =  The second derivative of specific enthalpy with respect to 
%              CT at constant SA and p.                      [ J/(kg K^2) ]
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
%  McDougall, T. J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%    See Eqns. (18) and (22)
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
   error('gsw_enthalpy_second_derivatives_CT_exact:  Requires three inputs')
end %if

if ~(nargout == 3)
   error('gsw_enthalpy_second_derivatives_CT_exact:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
   error('gsw_enthalpy_second_derivatives_CT_exact: SA and CT do not have the same dimensions')
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
    error('gsw_enthalpy_second_derivatives_CT_exact: The dimensions of p do not agree')
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

cp0 = gsw_cp0;                      % from Eqn. 3.3.3 of IOC et al. (2010).
pr0 = zeros(size(SA)); 
T0 = gsw_T0;

pt0 = gsw_pt_from_CT(SA,CT);
rec_abs_pt0 = 1./(T0 + pt0);
t = gsw_pt_from_t(SA,pt0,pr0,p);
temp_ratio = (T0 + t).*rec_abs_pt0;

rec_gTT_pt0 = 1./gsw_gibbs(0,2,0,SA,pt0,pr0);
rec_gTT = 1./gsw_gibbs(0,2,0,SA,t,p);
gSAT_pt0 = gsw_gibbs(1,1,0,SA,pt0,pr0);
gSAT = gsw_gibbs(1,1,0,SA,t,p);
gSA_pt0 = gsw_gibbs(1,0,0,SA,pt0,pr0);
gSASA = gsw_gibbs(2,0,0,SA,t,p);
gSASA_pt0 = gsw_gibbs(2,0,0,SA,pt0,pr0);

part_a = temp_ratio.*rec_gTT_pt0 - rec_gTT;

% h_CT_CT is naturally well-behaved as SA approaches zero. 
h_CT_CT = cp0.*cp0.*rec_abs_pt0.*rec_abs_pt0.*part_a;

part_b = rec_abs_pt0.*(temp_ratio.*gSAT_pt0.*rec_gTT_pt0 - gSAT.*rec_gTT);
factor = gSA_pt0./cp0;

% h_SA_CT should not blow up as SA approaches zero.  The following lines
% of code ensure that the h_SA_CT output of this function does not blow
% up in this limit.  That is, when SA < 1e-100 g/kg, we force the h_SA_CT 
% output to be the same as if SA = 1e-100 g/kg.  
if any(SA < 1e-100)
    SA(SA < 1e-100) = 1e-100;
    rec_gTT_pt0 = 1./gsw_gibbs(0,2,0,SA,pt0,pr0);
    rec_gTT = 1./gsw_gibbs(0,2,0,SA,t,p);
    gSAT_pt0 = gsw_gibbs(1,1,0,SA,pt0,pr0);
    gSAT = gsw_gibbs(1,1,0,SA,t,p);
    gSA_pt0 = gsw_gibbs(1,0,0,SA,pt0,pr0);
    part_b = (temp_ratio.*gSAT_pt0.*rec_gTT_pt0 - gSAT.*rec_gTT).*rec_abs_pt0;
    factor = gSA_pt0./cp0;
end

h_SA_CT  = cp0.*part_b - factor.*h_CT_CT;

% h_SA_SA has a singularity at SA = 0, and blows up as SA approaches zero.  
h_SA_SA = gSASA - temp_ratio.*gSASA_pt0  ...
    + temp_ratio.*gSAT_pt0.*gSAT_pt0.*rec_gTT_pt0  ...
    - gSAT.*gSAT.*rec_gTT  ...
    - 2.*gSA_pt0.*h_SA_CT./cp0 - (factor.*factor).*h_CT_CT;

if transposed
    h_SA_SA = h_SA_SA.';
    h_SA_CT = h_SA_CT.';
    h_CT_CT = h_CT_CT.';
end

end
