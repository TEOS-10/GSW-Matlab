function [SA_freeze, CT_freeze, w_seaice] = gsw_seaice_fraction_to_freeze_seawater(SA,CT,p,SA_seaice,t_seaice)

% gsw_seaice_fraction_to_freeze_seawater  sea ice mass fraction, which when
%                                          melted into seawater, brings the 
%                                      seawater to the freezing temperature
%==========================================================================
%
% USAGE:
%  [SA_freeze, CT_freeze, w_seaice] = ...
%        gsw_seaice_fraction_to_freeze_seawater(SA,CT,p,SA_seaice,t_seaice)
%
% DESCRIPTION:
%  Calculates the mass fraction of sea ice (mass of sea ice divided by mass 
%  of sea ice plus seawater), which, when melted into seawater having the
%  properties (SA,CT,p) causes the final seawater to be at the freezing 
%  temperature.  The other outputs are the Absolute Salinity and 
%  Conservative Temperature of the final seawater.  
%
% INPUT:
%  SA        =  Absolute Salinity of seawater                      [ g/kg ]
%  CT        =  Conservative Temperature of seawater (ITS-90)     [ deg C ]
%  p         =  sea pressure                                       [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction of             
%               salt in sea ice, expressed in g of salt per kg of sea ice.
%                                                                  [ g/kg ]
%  t_seaice  =  in-situ temperature of the sea ice at pressure p (ITS-90)
%                                                                 [ deg C ]
%
%  SA, CT, SA_seaice and t_seaice must all have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA, CT, SA_seaice
%  and t_seaice are MxN.
%
% OUTPUT:
%  SA_freeze  =  Absolute Salinity of seawater after the mass fraction of
%                sea ice, w_seaice, at temperature t_seaice has melted into
%                the original seawater, and the final mixture is at the 
%                freezing temperature of seawater.                 [ g/kg ]
%
%  CT_freeze  =  Conservative Temperature of seawater after the mass 
%                fraction, w_seaice, of sea ice at temperature t_seaice has
%                melted into the original seawater, and the final mixture 
%                is at the freezing temperature of seawater.      [ deg C ]
%
%  w_seaice   =  mass fraction of sea ice, at SA_seaice and t_seaice, 
%                which, when melted into seawater at (SA,CT,p) leads to the
%                final mixed seawater being at the freezing temperature.  
%                This output is between 0 and 1.                 [unitless]
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
%  McDougall T.J. and S.J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (23) of this manuscript.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 5)
    error('gsw_seaice_fraction_to_freeze_seawater:  Requires five inputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mssi,nssi] = size(SA_seaice);
[mtsi,ntsi] = size(t_seaice);

if (mt ~= ms) & (nt ~= ns)
    error('gsw_seaice_fraction_to_freeze_seawater: The inputs SA and CT must have the same dimensions')
end 

if (mssi ~= ms) & (nssi ~= ns)
    error('gsw_seaice_fraction_to_freeze_seawater: The inputs SA and SA_seaice must have the same dimensions')
end

if (mtsi ~= ms) & (ntsi ~= ns)
    error('gsw_seaice_fraction_to_freeze_seawater: The inputs SA and t_seaice must have the same dimensions')
end 

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms),:);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,ms),:);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_seaice_fraction_to_freeze_seawater: Inputs array dimensions arguments do not agree; check p')
end

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    SA_seaice = SA_seaice.';
    t_seaice = t_seaice.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This ensures that SA is non-negative.
SA(SA < 0) = NaN; 

SA0 = zeros(size(SA));

saturation_fraction = zeros(size(SA));

if any(SA_seaice(:) < 0 | SA_seaice(:) > 15) % SA_seaice must be between 0 and 15 g/kg
    [I] = find(SA_seaice < 0 | SA_seaice > 15);
    SA(I) = NaN;
end 

CTf = gsw_CT_freezing(SA,p,saturation_fraction);
if any(CT(:) < CTf(:))                % the seawater CT input is below the 
    [Icold] = find(CT < CTf);                       % freezing temperature
    SA(Icold) = NaN;
end

%--------------------------------------------------------------------------
tf_SA_seaice = gsw_t_freezing(SA_seaice,p,saturation_fraction) - 1e-6;
if any(t_seaice(:) > tf_SA_seaice(:))       % t_seaice exceeds the freezing
    [Iwarm] = find(t_seaice > tf_SA_seaice);               % temperature at 
    SA(Iwarm) = NaN;                                            % SA_seaice
end
% The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
% some ice Ih in the sea ice.   Without this buffer, that is if t_seaice
% is allowed to be exactly equal to 
% gsw_t_freezing(SA_seaice,p,saturation_fraction), the sea ice is 
% actually 100% brine at Absolute Salinity of SA_seaice.
%--------------------------------------------------------------------------

h = gsw_enthalpy_CT_exact(SA,CT,p);
h_Ih = gsw_enthalpy_ice(t_seaice,p);
SA_freezing = gsw_SA_freezing_from_t(t_seaice,p,saturation_fraction); % Note that gsw_SA_freezing_from_t returns a NaN if SA_freezing is greater than 120 g/kg. 
h_brine = gsw_enthalpy_t_exact(SA_freezing,t_seaice,p);
salt_ratio = SA_seaice./SA_freezing;

CTf_plus1 = gsw_CT_freezing(SA+1,p,saturation_fraction);
func_plus1 = (SA - SA_seaice).*(gsw_enthalpy_CT_exact(SA+1,CTf_plus1,p) - h) ...
               - (h - h_Ih) + salt_ratio.*(h_brine - h_Ih);

CTf_zero = gsw_CT_freezing(SA0,p,saturation_fraction);
func_zero = (SA - SA_seaice).*(gsw_enthalpy_CT_exact(SA0,CTf_zero,p) - h) ...
             + SA.*((h - h_Ih) - salt_ratio.*(h_brine - h_Ih));

SAf = -(SA+1).*func_zero./(func_plus1 - func_zero); % Initial guess of SAf
CTf = gsw_CT_freezing(SAf,p,saturation_fraction);
[h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SAf,CTf,p);
[CTf_SA, dummy] = gsw_CT_freezing_first_derivatives(SAf,p,saturation_fraction);

dfunc_dSAf = (SA - SA_seaice).*(h_hat_SA + h_hat_CT.*CTf_SA) ...
              - (h - h_Ih) + salt_ratio.*(h_brine - h_Ih);

for Number_of_iterations = 1:4
    SAf_old = SAf;   
    CTf_old = CTf;
    func = (SA - SA_seaice).*(gsw_enthalpy_CT_exact(SAf_old,CTf_old,p) - h) ...
         - (SAf_old - SA).*((h - h_Ih) - salt_ratio.*(h_brine - h_Ih));
    SAf = SAf_old - func./dfunc_dSAf ; % this is half way through the modified Newton's method
    SAf_mean = 0.5*(SAf + SAf_old);
    CTf_mean = gsw_CT_freezing(SAf_mean,p,saturation_fraction);
    [h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SAf_mean,CTf_mean,p);
    [CTf_SA, dummy] = gsw_CT_freezing_first_derivatives(SAf_mean,p,saturation_fraction);
    dfunc_dSAf = (SA - SA_seaice).*(h_hat_SA + h_hat_CT.*CTf_SA) ...
                  - (h - h_Ih) + salt_ratio.*(h_brine - h_Ih);
    SAf = SAf_old - func./dfunc_dSAf ; % this is the end of one full iteration of the modified Newton's method
    CTf = gsw_CT_freezing(SAf,p,saturation_fraction);
end

SA_freeze = SAf;
CT_freeze = CTf;
h_freeze = gsw_enthalpy_CT_exact(SAf,CTf,p);
w_seaice = (h - h_freeze)./(h - h_Ih - salt_ratio.*(h_brine - h_Ih));

if any(w_seaice(:) < 0 | w_seaice(:) > 1)
    [Iallice] = find(w_seaice < 0 | w_seaice > 1); %the output, w_Ih, is between 0 and 1
    w_seaice(Iallice) = NaN;
    SA_freeze(Iallice) = NaN;
    CT_freeze(Iallice) = NaN;
end

if transposed
    SA_freeze = SA_freeze.';
    CT_freeze = CT_freeze.';
    w_seaice = w_seaice.';
end

% After these 4 iterations of this modified Newton-Raphson method, the
% errors in SA_freeze is less than 1.5x10^-12 g/kg, in CT_freeze is less than
% 2x10^-13 deg C and in w_seaice is less than 2.8x10^-13 which represent machine 
% precision for these calculations. 

end