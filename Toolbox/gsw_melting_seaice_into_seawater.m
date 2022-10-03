function [SA_final, CT_final] = gsw_melting_seaice_into_seawater(SA,CT,p,w_seaice,SA_seaice,t_seaice)

% gsw_melting_seaice_into_seawater             the resulting SA and CT when 
%                                           sea ice is melted into seawater
%==========================================================================
%
% USAGE:
%  [SA_final, CT_final] = ...
%     gsw_melting_seaice_into_seawater(SA,CT,p,w_seaice,SA_seaice,t_seaice)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity and Conservative Temperature that 
%  results when a given mass of sea ice (or ice) melts and is mixed into a 
%  known mass of seawater (whose properties are (SA,CT,p)).  
%
%  If the ice contains no salt (e.g. if it is of glacial origin), then the 
%  input 'SA_seaice' should be set to zero.  
%
%  Ice formed at the sea surface (sea ice) typically contains between 2 g/kg
%  and 12 g/kg of salt (defined as the mass of salt divided by the mass of 
%  ice Ih plus brine) and this programme returns NaN's if the input  
%  SA_seaice is greater than 15 g/kg.  If the SA_seaice input is not zero,   
%  usually this would imply that the pressure p should be zero, as sea ice  
%  only occurs near the sea surface.  The code does not impose that p = 0 
%  if SA_seaice is non-zero.  Rather, this is left to the user.  
%
%  The Absolute Salinity, SA_brine, of the brine trapped in little pockets 
%  in the sea ice, is in thermodynamic equilibrium with the ice Ih that
%  surrounds these pockets.  As the sea ice temperature, t_seaice, may be 
%  less than the freezing temperature, SA_brine is usually greater than the
%  Absolute Salinity of the seawater at the time and place when and where 
%  the sea ice was formed.  So usually SA_brine will be larger than SA.  
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%  w_seaice  =  mass fraction of sea ice, that is the mass of sea ice 
%               divided by the sum of the masses of sea ice and seawater. 
%               That is, the mass of sea ice divided by the mass of the 
%               final mixed fluid.  w_seaice must be between 0 and 1. 
%                                                              [ unitless ]
%  SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction of             
%               salt in sea ice, expressed in g of salt per kg of sea ice.
%                                                                  [ g/kg ]
%  t_seaice  =  the in-situ temperature of the sea ice (or ice) (ITS-90)
%                                                                 [ deg C ]
%
% SA, CT, w_seaice, SA_seaice & t_seaice must all have the same dimensions.
% p may have dimensions 1x1 or Mx1 or 1xN or MxN, where where SA, CT, 
% w_seaice, SA_seaice and t_seaice are MxN.
%
% OUTPUT:
%  SA_final  =  Absolute Salinity of the mixture of the melted sea ice 
%               (or ice) and the orignal seawater                  [ g/kg ]
%  CT_final  =  Conservative Temperature of the mixture of the melted 
%               sea ice (or ice) and the orignal seawater         [ deg C ]            
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
%     Eqns. (8) and (9) are the simplifications when SA_seaice = 0. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 6)
    error('gsw_melting_seaice_into_seawater: Requires six inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mw_Ih,nw_Ih] = size(w_seaice);
[ms_Ih,ns_Ih] = size(SA_seaice);
[mt_Ih,nt_Ih] = size(t_seaice);

if (mt ~= ms | nt ~= ns)
    error('gsw_melting_seaice_into_seawater: SA and CT must have same dimensions')
end

if (ms_Ih ~= ms | ns_Ih ~= ns)
    error('gsw_melting_seaice_into_seawater: SA and SA_seaice must have same dimensions')
end

if (mw_Ih ~= ms | nw_Ih ~= ns)
    error('gsw_melting_seaice_into_seawater: SA and w_seaice must have same dimensions')
end

if (mt_Ih ~= ms | nt_Ih ~= ns)
    error('gsw_melting_seaice_into_seawater: SA and t_seaice must have same dimensions')
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
    error('gsw_melting_seaice_into_seawater: Inputs array dimensions arguments do not agree; check p')
end 

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    w_seaice = w_seaice.';
    SA_seaice = SA_seaice.';
    t_seaice = t_seaice.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensures that SA is non-negative.

saturation_fraction = zeros(size(SA));

if any(w_seaice(:) < 0 | w_seaice(:) > 1) % the w_seaice needs to be between 0 and 1
    [I] = find(w_seaice > 1 | w_seaice < 0);
    SA(I) = NaN;
    p(I) = NaN;
end 

if any(SA_seaice(:) < 0 | SA_seaice(:) > 15) % the SA_seaice needs to be between 0 and 15 g/kg
    [I] = find(SA_seaice > 15 | SA_seaice < 0);
    SA(I) = NaN;
    p(I) = NaN;
end

CTf = gsw_CT_freezing(SA,p,saturation_fraction);
if any(CT(:) < CTf(:)) % the seawater CT input is below the freezing temperature
    [I] = find(CT < CTf);
    SA(I) = NaN;
    p(I) = NaN;
end

%--------------------------------------------------------------------------
tf_SA_seaice = gsw_t_freezing(SA_seaice,p,saturation_fraction) - 1e-6;
if any(t_seaice(:) > tf_SA_seaice(:))       % t_seaice exceeds the freezing
    [Iwarm] = find(t_seaice > tf_SA_seaice);               % temperature at 
    SA(Iwarm) = NaN;                                            % SA_seaice
    p(Iwarm) = NaN;
end
% The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
% some ice Ih in the sea ice.   Without this buffer, that is if t_seaice
% is allowed to be exactly equal to gsw_t_freezing(SA_seaice,p,0), the
% seaice is actually 100% brine at Absolute Salinity of SA_seaice.
%--------------------------------------------------------------------------

h = gsw_enthalpy_CT_exact(SA,CT,p);
h_Ih = gsw_enthalpy_ice(t_seaice,p);
SA_brine = gsw_SA_freezing_from_t(t_seaice,p,saturation_fraction);
h_brine = gsw_enthalpy_t_exact(SA_brine,t_seaice,p);

SA_final = SA - w_seaice.*(SA - SA_seaice);

h_final = nan(size(SA_final));
h_final(SA_seaice == 0) = h(SA_seaice == 0) - w_seaice(SA_seaice == 0).*(h(SA_seaice == 0) - h_Ih(SA_seaice == 0));
h_final(SA_seaice > 0) = h(SA_seaice > 0) - w_seaice(SA_seaice > 0).*(h(SA_seaice > 0) - h_Ih(SA_seaice > 0) ...
    - (h_brine(SA_seaice > 0) - h_Ih(SA_seaice > 0)).*SA_seaice(SA_seaice > 0)./SA_brine(SA_seaice > 0));

if any(isnan(h_final(:)) | isnan(SA_final(:)))
    [Inan] = find(isnan(h_final) | isnan(SA_final));
    SA_final(Inan) = NaN;
    h_final(Inan) = NaN;
end

CTf = gsw_CT_freezing(SA_final,p,saturation_fraction);
hf = gsw_enthalpy_CT_exact(SA_final,CTf,p);

if any(h_final(:) < hf(:)) % Melting this much seaice is not possible as it would result in frozen seawater
    [I] = find(h_final < hf);
    SA_final(I) = NaN;
end

CT_final = gsw_CT_from_enthalpy_exact(SA_final,h_final,p);

if transposed
   SA_final = SA_final.';
   CT_final = CT_final.';
end

end
