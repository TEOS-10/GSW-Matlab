function [dSA_dCT_frazil, dSA_dP_frazil, dCT_dP_frazil] = gsw_frazil_ratios_adiabatic(SA,p,w_Ih)

% gsw_frazil_ratios_adiabatic    ratios of SA, CT and P changes when frazil
%                                  ice forms due to the adiabatic change in 
%                                 pressure of a mixture of seawater and ice
%==========================================================================
%
% USAGE:
%  [dSA_dCT_frazil, dSA_dP_frazil, dCT_dP_frazil] = ...
%                                    gsw_frazil_ratios_adiabatic(SA,p,w_Ih)
%
% DESCRIPTION:
%  Calculates the ratios of SA, CT and P changes when frazil ice forms or
%  melts in response to an adiabatic change in pressure of a mixture of 
%  seawater and frazil ice crystals.  
%
%  Note that the first output, dSA_dCT_frazil, is dSA/dCT rather than 
%  dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT, is zero 
%  whereas dCT/dSA would then be infinite. 
%
%  Also note that both dSA_dP_frazil and dCT_dP_frazil are the pressure
%  derivatives with the pressure measured in Pa not dbar.  
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  p   =  sea pressure of seawater at which melting occurs         [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%  w_Ih  =  mass fraction of ice, that is the mass of ice divided by the 
%           sum of the masses of ice and seawater.  That is, the mass of
%           ice divided by the mass of the final mixed fluid.  
%           w_Ih must be between 0 and 1.                      [ unitless ]
%
% SA & w_Ih must have the same dimensions.
% p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA and w_Ih are MxN.
%
% OUTPUT:
%  dSA_dCT_frazil =  the ratio of the changes in Absolute Salinity 
%                    to that of Conservative Temperature       [ g/(kg K) ] 
%  dSA_dP_frazil  =  the ratio of the changes in Absolute Salinity 
%                    to that of pressure (in Pa)              [ g/(kg Pa) ] 
%  dCT_dP_frazil  =  the ratio of the changes in Conservative Temperature
%                    to that of pressure (in Pa)                   [ K/Pa ]           
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqns. (47), (48) and (49) of this manuscript.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3 )
    error('gsw_frazil_ratios_adiabatic: Requires three inputs')
end

[ms,ns] = size(SA);
[mp,np] = size(p);
[mw_Ih,nw_Ih] = size(w_Ih);

if (mw_Ih ~= ms | nw_Ih ~= ns)
    error('gsw_frazil_ratios_adiabatic: SA and w_Ih must have same dimensions')
end

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                              % transposed then
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_frazil_ratios_adiabatic: Inputs array dimensions arguments do not agree; check p')
end 

if ms == 1
    SA = SA.';
    p = p.';
    w_Ih = w_Ih.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensure that SA is non-negative.

if any(w_Ih(:) < 0 | w_Ih(:) > 1) % the w_Ih must be between 0 and 1.
    [I] = find(w_Ih < 0 | w_Ih > 1);
    SA(I) = NaN;
end

saturation_fraction = zeros(size(SA));

tf = gsw_t_freezing(SA,p,saturation_fraction);
CTf = gsw_CT_freezing(SA,p,saturation_fraction);
h = gsw_enthalpy_CT_exact(SA,CTf,p);
h_Ih = gsw_enthalpy_ice(tf,p);
cp_Ih = gsw_cp_ice(tf,p);
Gamma_Ih = gsw_adiabatic_lapse_rate_ice(tf,p);
[h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CTf,p);
[tf_SA, tf_P] = gsw_t_freezing_first_derivatives(SA,p,saturation_fraction);
[CTf_SA, CTf_P] = gsw_CT_freezing_first_derivatives(SA,p,saturation_fraction);

wcp = cp_Ih.*w_Ih./(1 - w_Ih);
part = (tf_P - Gamma_Ih)./CTf_P;

bracket1 = h_hat_CT + wcp.*part;
bracket2 = h - h_Ih - SA.*(h_hat_SA + wcp.*(tf_SA - part.*CTf_SA));
rec_bracket3 = 1./(h - h_Ih - SA.*(h_hat_SA + h_hat_CT.*CTf_SA + wcp.*tf_SA));

dSA_dCT_frazil = SA.*(bracket1./bracket2);
dSA_dP_frazil = SA.*CTf_P.*bracket1.*rec_bracket3;
dCT_dP_frazil = CTf_P.*bracket2.*rec_bracket3;

if transposed
    dSA_dCT_frazil = dSA_dCT_frazil.';
    dSA_dP_frazil = dSA_dP_frazil.';
    dCT_dP_frazil = dCT_dP_frazil.';
end

end