function [SA_final, CT_final, w_Ih_final] = gsw_melting_ice_into_seawater(SA,CT,p,w_Ih,t_Ih)

% gsw_melting_ice_into_seawater             Absolute Salinity, Conservative
%                  Temperature and final ice mass fraction when ice of mass
%                fraction w_Ih and temperature t_Ih is melted into seawater
%==========================================================================
%
% USAGE:
%  [SA_final, CT_final, w_Ih_final] = ... 
%                          gsw_melting_ice_into_seawater(SA,CT,p,w_Ih,t_Ih)
%
% DESCRIPTION:
%  Calculates the final Absolute Salinity, final Conservative Temperature 
%  and final ice mass fraction that results when a given mass fraction of 
%  ice melts and is mixed into seawater whose properties are (SA,CT,p).  
%  This code takes the seawater to contain no dissolved air.  
%
%  When the mass fraction w_Ih_final is calculated as being a positive
%  value, the seawater-ice mixture is at thermodynamic equlibrium.  
%
%  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk, 
%  is sufficiently large (i.e. sufficiently "warm") so that there is no ice 
%  present in the final state.  In this case the final state consists of 
%  only seawater rather than being an equlibrium mixture of seawater and 
%  ice which occurs when w_Ih_final is positive.  Note that when 
%  w_Ih_final = 0, the final seawater is not at the freezing temperature. 
%
% INPUT:
%  SA   =  Absolute Salinity of seawater                           [ g/kg ]
%  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
%  p    =  sea pressure at which the melting occurs                [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%  w_Ih =  mass fraction of ice, that is the mass of ice divided by the
%          sum of the masses of ice and seawater.  That is, the mass of 
%          ice divided by the mass of the final mixed fluid.  
%          w_Ih must be between 0 and 1.                       [ unitless ]
%  t_Ih =  the in-situ temperature of the ice (ITS-90)            [ deg C ]
%
%  SA, CT, w_Ih and t_Ih must all have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where where SA, CT, 
%  w_Ih and t_Ih are MxN.
%
% OUTPUT:
%  SA_final    =  Absolute Salinity of the seawater in the final state, 
%                 whether or not any ice is present.               [ g/kg ]
%  CT_final    =  Conservative Temperature of the seawater in the the final
%                 state, whether or not any ice is present.       [ deg C ]
%  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
%                 If this ice mass fraction is positive, the system is at 
%                 thermodynamic equilibrium.  If this ice mass fraction is 
%                 zero there is no ice in the final state which consists 
%                 only of seawater which is warmer than the freezing 
%                 temperature.                                   [unitless]
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
%   Melting of ice and sea ice into seawater, and frazil ice formation.
%   Journal of Physical Oceanography, 44, 1751-1775. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 5)
    error('gsw_melting_ice_into_seawater: Requires five inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
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

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    w_Ih = w_Ih.';
    t_Ih = t_Ih.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensures that SA is non-negative.

saturation_fraction = zeros(size(SA)); % Throughout this code seawater is
                                       % taken to contain no dissolved air.  

if any(w_Ih(:) < 0 | w_Ih(:) > 1) % the w_Ih needs to be between 0 and 1
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

tf_Ih = gsw_t_freezing(zeros(size(p)),p,saturation_fraction) - 1e-6;
if any(t_Ih(:) > tf_Ih(:))       % t_Ih exceeds the freezing temperature
    [Iwarm] = find(t_Ih > tf_Ih);                
    SA(Iwarm) = NaN;                                            
    p(Iwarm) = NaN;
end
% The 1e-6 C buffer in the allowable t_Ih is to ensure that there is
% some ice Ih in the sea ice.  
%--------------------------------------------------------------------------
        
SA_bulk = (1 - w_Ih).*SA;
h_bulk = (1 - w_Ih).*gsw_enthalpy_CT_exact(SA,CT,p) + w_Ih.*gsw_enthalpy_ice(t_Ih,p);
[SA_final, CT_final, w_Ih_final] = gsw_frazil_properties(SA_bulk,h_bulk,p);
        
if transposed
   SA_final = SA_final.';
   CT_final = CT_final.';
   w_Ih_final = w_Ih_final.';
end

end
