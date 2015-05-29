function [SA_final, CT_final] = gsw_melting_ice_into_seawater(SA,CT,p,saturation_fraction,w_Ih,t_Ih)

% gsw_melting_ice_into_seawater                the resulting SA and CT when 
%                                               ice is melted into seawater
%==========================================================================
%
% USAGE:
%  [SA_final, CT_final] = gsw_melting_ice_into_seawater(SA,CT,p,saturation_fraction,w_Ih,t_Ih)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity and Conservative Temperature that 
%  results when a given mass of ice melts and is mixed into a known mass of
%  seawater (whose properties are (SA,CT,p)).  
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%  saturation_fraction = the saturation fraction of dissolved air in 
%               seawater.  The saturation_fraction must be between 0 and 1.
%  w_Ih  =  mass fraction of ice, that is the mass of ice divided by the
%           sum of the masses of ice and seawater.  That is, the mass of 
%           ice divided by the mass of the final mixed fluid.  
%           w_Ih must be between 0 and 1.                      [ unitless ]
%  t_Ih  =  the in-situ temperature of the ice (ITS-90)           [ deg C ]
%
%  SA, CT, w_Ih and t_Ih must all have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where where SA, CT, 
%  w_Ih and t_Ih are MxN.
%
% OUTPUT:
%  SA_final  =  Absolute Salinity of the mixture of the melted ice and the
%               orignal seawater                                   [ g/kg ]
%  CT_final  =  Conservative Temperature of the mixture of the melted ice 
%               and the orignal seawater                          [ deg C ]            
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.04 (3rd January, 2014)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., P.M. Barker and R. Feistel, 2013: Melting of ice and 
%   sea ice into seawater and frazil ice formation. Journal of Physical 
%   Oceanography, (Submitted).
%    Eqns. (8) and (9) are the simplifications when SA_seaice = 0. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 6)
    error('gsw_melting_ice_into_seawater: Requires six inputs')
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_melting_ice_into_seawater: saturation fraction MUST be between zero and one.')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[msf,nsf] = size(saturation_fraction);
[mw_Ih,nw_Ih] = size(w_Ih);
[mt_Ih,nt_Ih] = size(t_Ih);

if (mt ~= ms | nt ~= ns)
    error('gsw_melting_ice_into_seawater: SA and CT must have same dimensions')
end

if (mw_Ih ~= ms | nw_Ih ~= ns)
    error('gsw_melting_ice_into_seawater: SA and w_ice must have same dimensions')
end

if (mt_Ih ~= ms | nt_Ih ~= ns)
    error('gsw_melting_ice_into_seawater: SA and t_ice must have same dimensions')
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
    error('gsw_melting_ice_into_seawater: Inputs array dimensions arguments do not agree; check p')
end 

if (msf == 1) & (nsf == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(SA));         % fill to size of SA
elseif (ns == nsf) & (msf == 1)                        % saturation_fraction is row vector,
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (nsf == 1)                     % saturation_fraction is column vector,
    saturation_fraction = saturation_fraction(:,ones(1,ns));        % copy across each row.
elseif (ns == msf) & (nsf == 1)           % saturation_fraction is a transposed row vector,
    saturation_fraction = saturation_fraction.';                           % transposed then
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (ns == nsf)
    % ok
else
    error('gsw_melting_ice_into_seawater: Inputs array dimensions arguments do not agree')
end %if


if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    w_ice = w_Ih.';
    t_ice = t_Ih.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensures that SA is non-negative.

if any(w_Ih(:) < 0 | w_Ih(:) > 1) % the w_ice needs to be between 0 and 1
    [I] = find(w_Ih > 1 | w_Ih < 0);
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
tf_Ih = gsw_t_freezing(zeros(size(p)),p,saturation_fraction) - 1e-6;
if any(t_Ih(:) > tf_Ih(:))       % t_Ih exceeds the freezing temperature
    [Iwarm] = find(t_Ih > tf_Ih);                
    SA(Iwarm) = NaN;                                            
    p(Iwarm) = NaN;
end
% The 1e-6 C buffer in the allowable t_Ih is to ensure that there is
% some ice Ih in the sea ice.   Without this buffer, that is if t_Ih
% is allowed to be exactly equal to gsw_t_freezing(0,p,0).
%--------------------------------------------------------------------------

h = gsw_enthalpy_CT_exact(SA,CT,p);
h_Ih = gsw_enthalpy_ice(t_Ih,p);

SA_final = SA.*(1 - w_Ih);

h_final = h - w_Ih.*(h - h_Ih);

if any(isnan(h_final(:)) | isnan(SA_final(:)))
    [Inan] = find(isnan(h_final) | isnan(SA_final));
    SA_final(Inan) = NaN;
    h_final(Inan) = NaN;
end

CTf = gsw_CT_freezing(SA_final,p,saturation_fraction);
hf = gsw_enthalpy_CT_exact(SA_final,CTf,p);

if any(h_final(:) < hf(:)) % Melting this much ice is not possible as it would result in frozen seawater
    [I] = find(h_final < hf);
    SA_final(I) = NaN;
end

CT_final = gsw_CT_from_enthalpy_exact(SA_final,h_final,p);
         
if transposed
   SA_final = SA_final.';
   CT_final = CT_final.';
end

end
