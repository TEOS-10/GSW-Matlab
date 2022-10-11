function [SA_freeze, CT_freeze, w_Ih] = gsw_ice_fraction_to_freeze_seawater(SA,CT,p,t_Ih)

% gsw_ice_fraction_to_freeze_seawater         ice mass fraction, which when
%                                  melted into seawater, brings the diluted
%                                      seawater to the freezing temperature        
%==========================================================================
%
% USAGE:
%  [SA_freeze, CT_freeze, w_Ih] = ...
%                         gsw_ice_fraction_to_freeze_seawater(SA,CT,p,t_Ih)
%
% DESCRIPTION:
%  Calculates the mass fraction of ice (mass of ice divided by mass of ice
%  plus seawater), which, when melted into seawater having (SA,CT,p) causes 
%  the final dilute seawater to be at the freezing temperature.  The other
%  outputs are the Absolute Salinity and Conservative Temperature of the
%  final diluted seawater.  
%
% INPUT:
%  SA   =  Absolute Salinity of seawater                           [ g/kg ]
%  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
%  p    =  sea pressure                                            [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  t_Ih =  in-situ temperature of the ice at pressure p (ITS-90)  [ deg C ]
%
%  SA, CT and t_Ih must have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA, CT and t_Ih 
%  are MxN.
%
% OUTPUT:
%  SA_freeze = Absolute Salinity of seawater after the mass fraction of 
%              ice, ice_fraction, at temperature t_Ih has melted into the
%              original seawater, and the final mixture is at the freezing
%              temperature of seawater.                            [ g/kg ]
%
%  CT_freeze = Conservative Temperature of seawater after the mass 
%              fraction, w_Ih, of ice at temperature t_Ih has melted into
%              the original seawater, and the final mixture is at the
%              freezing temperature of seawater.                  [ deg C ]
%
%  w_Ih      = mass fraction of ice, having in-situ temperature t_Ih, 
%              which, when melted into seawater at (SA,CT,p) leads to the
%              final diluted seawater being at the freezing temperature.
%              This output must be between 0 and 1.              [unitless]
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
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.  
%
%  McDougall, T.J., and S.J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (9) of this manuscript.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4) 
   error('gsw_ice_fraction_to_freeze_seawater:  Requires four inputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mt_Ih,nt_Ih] = size(t_Ih);

if (mt ~= ms) & (nt ~= ns)
    error('gsw_ice_fraction_to_freeze_seawater: The inputs SA and CT must have the same dimensions')
end %if

if (mt_Ih ~= ms) & (nt_Ih ~= ns)
    error('gsw_ice_fraction_to_freeze_seawater: The inputs SA and t_Ih must have the same dimensions')
end %if

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_ice_fraction_to_freeze_seawater: Inputs array dimensions arguments do not agree; check p')
end 

if ms == 1
    SA = SA.';
    CT = CT.';  
    p = p.';
    t_Ih = t_Ih.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensure that SA is non-negative.

saturation_fraction = zeros(size(SA));

CTf = gsw_CT_freezing(SA,p,saturation_fraction);
if any(CT(:) < CTf(:)) %the seawater CT input is below the freezing temperature
    [I] = find(CT < CTf);
    SA(I) = NaN;
end

SA0 = zeros(size(SA));
tf = gsw_t_freezing(SA0,p,saturation_fraction);
if (t_Ih(:) > tf(:)) % the input, t_Ih, exceeds the freezing temperature at SA = 0
    [I] = find(t_Ih > tf);
    SA(I) = NaN;
end

h = gsw_enthalpy_CT_exact(SA,CT,p);
h_Ih = gsw_enthalpy_ice(t_Ih,p);

CTf_plus1 = gsw_CT_freezing(SA+1,p,saturation_fraction);
func_plus1 = SA.*(gsw_enthalpy_CT_exact(SA+1,CTf_plus1,p) - h) - (h - h_Ih);

CTf_zero = gsw_CT_freezing(SA0,p,saturation_fraction);
func_zero = SA.*(gsw_enthalpy_CT_exact(SA0,CTf_zero,p) - h_Ih);

SAf = - (SA+1).*func_zero./(func_plus1 - func_zero); % Initial guess of SA_freeze
CTf = gsw_CT_freezing(SAf,p,saturation_fraction);
[h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SAf,CTf,p);
[CTf_SA, dummy] = gsw_CT_freezing_first_derivatives(SAf,p);

dfunc_dSAf = SA.*(h_hat_SA + h_hat_CT.*CTf_SA) - (h - h_Ih);

for Number_of_iterations = 1:2
    SAf_old = SAf;
    CTf_old = CTf;
    func = SA.*(gsw_enthalpy_CT_exact(SAf_old,CTf_old,p) - h)  - (SAf_old - SA).*(h - h_Ih); 
    % This is the function whose zero we seek, it is SA*(Eqn. (9), (McDougall et al., 1014)) .
    SAf = SAf_old - func./dfunc_dSAf; % this is half way through the modified Newton's method
    SAf_mean = 0.5*(SAf + SAf_old);
    CTf_mean = gsw_CT_freezing(SAf_mean,p,saturation_fraction);
    [h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SAf_mean,CTf_mean,p);
    [CTf_SA, dummy] = gsw_CT_freezing_first_derivatives(SAf_mean,p,saturation_fraction);
    dfunc_dSAf = SA.*(h_hat_SA + h_hat_CT.*CTf_SA) - (h - h_Ih);
    SAf = SAf_old - func./dfunc_dSAf; 
    CTf = gsw_CT_freezing(SAf,p,saturation_fraction); % this is the end of one full iteration of the modified Newton's method
end

SA_freeze = SAf;
CT_freeze = CTf;
h_freeze = gsw_enthalpy_CT_exact(SA_freeze,CT_freeze,p);
w_Ih = (h - h_freeze)./(h - h_Ih);

if any(w_Ih(:) < 0 | w_Ih(:) > 1) % the output, w_Ih, must be between 0 and 1
    [I] = find(w_Ih < 0 | w_Ih > 1);
    SA_freeze(I) = NaN;
    CT_freeze(I) = NaN;
    w_Ih(I) = NaN;
end 

if transposed
    SA_freeze = SA_freeze.';
    CT_freeze = CT_freeze.';
    w_Ih = w_Ih.';
end

% After these 2 iterations of this modified Newton-Raphson method, the
% errors in SA_freeze is less than 1.3x10^-13 g/kg, in CT_freeze is less than
% 4x10^-13 deg C and in w_Ih is less than 3.8x10^-15 which represent machine 
% precision for these calculations. 

end