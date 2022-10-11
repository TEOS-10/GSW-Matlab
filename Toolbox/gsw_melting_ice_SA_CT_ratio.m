function melting_ice_SA_CT_ratio = gsw_melting_ice_SA_CT_ratio(SA,CT,p,t_Ih)

% gsw_melting_ice_SA_CT_ratio                ratio of SA to CT changes when
%                                                   ice melts into seawater
%==========================================================================
%
% USAGE:
%  melting_ice_SA_CT_ratio = gsw_melting_ice_SA_CT_ratio(SA,CT,p,t_Ih)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when ice melts into seawater.
%  It is assumed that a small mass of ice melts into an infinite mass of
%  seawater.  Because of the infinite mass of seawater, the ice will always
%  melt.   
%
%  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA. 
%  This is done so that when SA = 0, the output, dSA/dCT is zero whereas 
%  dCT/dSA would be infinite. 
%
% INPUT:
%  SA   =  Absolute Salinity of seawater                           [ g/kg ]
%  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
%  p    =  sea pressure at which the melting occurs                [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%  t_Ih =  the in-situ temperature of the ice (ITS-90)            [ deg C ]
%
% SA, CT & t_Ih must all have the same dimensions.
% p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA, CT and t_Ih 
% are MxN.
%
% OUTPUT:
%  melting_ice_SA_CT_ratio = the ratio of SA to CT changes when ice melts
%                            into a large mass of seawater 
%                                                          [ g kg^-1 K^-1 ]               
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
%    See Eqn. (13) of this manuscript.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4) 
   error('gsw_melting_ice_SA_CT_ratio: Requires four inputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mti,nti] = size(t_Ih);

if (mt ~= ms | nt ~= ns)
    error('gsw_melting_ice_SA_CT_ratio: SA and CT must have same dimensions')
end

if (mti ~= ms | nti ~= ns)
    error('gsw_melting_ice_SA_CT_ratio: SA and t_Ih must have same dimensions')
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
    error('gsw_melting_ice_SA_CT_ratio: Inputs array dimensions do not agree; check p')
end %if

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
if any(CT(:) < CTf(:)) % the seawater CT input is below the freezing temperature
    [I] = find(CT < CTf);
    SA(I) = NaN;
end

tf = gsw_t_freezing(zeros(size(SA)),p,saturation_fraction);
if any(t_Ih(:) >  tf(:)) % t_Ih exceeds the freezing temperature at SA = 0
    [I] = find(t_Ih < tf);
    SA(I) = NaN;
end

h = gsw_enthalpy_CT_exact(SA,CT,p);
h_Ih = gsw_enthalpy_ice(t_Ih,p);
[h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT,p);
% Note that h_hat_CT is equal to cp0*(273.15 + t)./(273.15 + pt0)

denominator = h - h_Ih - SA.*h_hat_SA;
melting_ice_SA_CT_ratio = SA.*h_hat_CT./denominator;

if transposed
    melting_ice_SA_CT_ratio = melting_ice_SA_CT_ratio.';
end

end
