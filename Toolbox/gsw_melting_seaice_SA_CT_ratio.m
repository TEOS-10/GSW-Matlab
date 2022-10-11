function melting_seaice_SA_CT_ratio = gsw_melting_seaice_SA_CT_ratio(SA,CT,p,SA_seaice,t_seaice)

% gsw_melting_seaice_SA_CT_ratio             ratio of SA to CT changes when
%                                               sea ice melts into seawater
%==========================================================================
%
% USAGE:
%  melting_seaice_SA_CT_ratio = ...
%                gsw_melting_seaice_SA_CT_ratio(SA,CT,p,SA_seaice,t_seaice)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when sea ice melts into 
%  seawater.  It is assumed that a small mass of sea ice melts into an 
%  infinite mass of seawater.  Because of the infinite mass of seawater, 
%  the sea ice will always melt.   
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
%  surrounds these pockets.  As the seaice temperature, t_seaice, may be 
%  less than the freezing temperature, SA_brine is usually greater than the
%  Absolute Salinity of the seawater at the time and place when and where 
%  the sea ice was formed.  So usually SA_brine will be larger than SA.  
%
%  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA. 
%  This is done so that when (SA - seaice_SA) = 0, the output, dSA/dCT is 
%  zero whereas dCT/dSA would be infinite. 
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%  SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction 
%               of salt in sea ice expressed in g of salt per kg of 
%               sea ice                                            [ g/kg ]
%  t_seaice =   the in-situ temperature of the sea ice (ITS-90)   [ deg C ]
%
% SA, CT, SA_seaice & t_seaice must all have the same dimensions.
% p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA, CT, SA_seaice
% and t_seaice are MxN.
%
% OUTPUT:
%  melting_seaice_SA_CT_ratio = the ratio dSA/dCT of SA to CT changes when
%                sea ice melts into a large mass of seawater   [ g/(kg K) ]               
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
%    See Eqn. (28) of this manuscript.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 5)
    error('gsw_melting_seaice_SA_CT_ratio: Requires five inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mssi,nssi] = size(SA_seaice);
[mtsi,ntsi] = size(t_seaice);

if (mt ~= ms | nt ~= ns)
    error('gsw_melting_seaice_SA_CT_ratio: SA and CT must have same dimensions')
end

if (mssi ~= ms | nssi ~= ns)
    error('gsw_melting_seaice_SA_CT_ratio: SA and SA_seaice must have same dimensions')
end

if (mtsi ~= ms | ntsi ~= ns)
    error('gsw_melting_seaice_SA_CT_ratio: SA and t_seaice must have same dimensions')
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
    error('gsw_melting_seaice_SA_CT_ratio: Inputs array dimensions do not agree; check p')
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

SA(SA < 0) = 0; % This line ensure that SA is non-negative.

saturation_fraction = zeros(size(SA));

if any(SA_seaice(:) < 0 | SA_seaice(:) > 15) %the SA_seaice input must be between 0 and 15
    [I] = find(SA_seaice < 0 | SA_seaice > 15);
    SA(I) = NaN;
end

CTf = gsw_CT_freezing(SA,p,saturation_fraction);
if any(CT(:) < CTf(:)) % the seawater CT input is below the freezing temperature
    [I] = find(CT < CTf);
    SA(I) = NaN;
end 

%--------------------------------------------------------------------------
tf_SA_seaice = gsw_t_freezing(SA_seaice,p,saturation_fraction) - 1e-6;
if any(t_seaice(:) > tf_SA_seaice(:))       % t_seaice exceeds the freezing
    [I] = find(t_seaice > tf_SA_seaice);                   % temperature at 
    SA(I) = NaN;                                                % SA_seaice
end
% The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
% some ice Ih in the sea ice.  Without this buffer, that is if t_seaice
% is allowed to be exactly equal to 
% gsw_t_freezing(SA_seaice,p,saturation_fraction), the sea ice is actually
% 100% brine (that is 100% seawater) at Absolute Salinity of SA_seaice.
%--------------------------------------------------------------------------

h = gsw_enthalpy_CT_exact(SA,CT,p);
h_Ih = gsw_enthalpy_ice(t_seaice,p);
[h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT,p);
% Note that h_hat_CT is equal to cp0*(273.15 + t)./(273.15 + pt0)

SA_brine = gsw_SA_freezing_from_t(t_seaice,p,saturation_fraction);
h_brine = gsw_enthalpy_t_exact(SA_brine,t_seaice,p);
delSA = SA - SA_seaice;

denominator = h - h_Ih - delSA.*h_hat_SA - SA_seaice.*(h_brine - h_Ih)./SA_brine;
melting_seaice_SA_CT_ratio = h_hat_CT.*delSA./denominator;

if transposed
    melting_seaice_SA_CT_ratio = melting_seaice_SA_CT_ratio.';
end

end