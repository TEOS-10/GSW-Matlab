function [SA_i, CT_i] = gsw_pchip_interp_SA_CT(SA,CT,p,p_i)

% gsw_pchip_interp_SA_CT              Piecewise Cubic Hermite Interpolating
%                                Polynomial, interpolation to p_i on a cast
%==========================================================================
%
% USAGE:
%  [SA_i, CT_i] = gsw_pchip_interp_SA_CT(SA,CT,p,p_i)
%
% DESCRIPTION:
%  This function interpolates values of SA, CT of the cast to the pressures
%  at p_i.  This programme provides an interpolation scheme based on 
%  shape-preserving piecewise cubic Hermite interpolating polynomial with
%  continuous first derivatives.  
%
% INPUT: 
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  CT   =  Conservative Temperature (ITS-90)                      [ deg C ]
%  p    =  sea pressure                                            [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i  =  pressures to interpolate to.
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  SA_i = interpolated SA values at pressures p_i.
%  CT_i = interpolated CT values at pressures p_i.
%
% AUTHOR:
%  This function was adapted from Matlab's pchip.
%
% MODIFIED:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th August 2016)
%
% References:
%  Forsythe, G.E., D. Kahaner, C.B. Moler and S.G. Nash, 1988: Numerical 
%   Methods and Software. Englewood Cliffs, N.J., Prentice Hall. 495 p.
%
%  Fritsch, F.N., and R.E. Carlson, 1980: Monotone Piecewise Cubic
%   Interpolation. SIAM J. Numerical Analysis, 17, 238-246.
%
%
%==========================================================================

SA = SA.';
CT = CT.';
p = p.';

[min_p,Imin_p] = min(p);

dp = diff(p);
dSA_dp = diff(SA)./dp;
dCT_dp = diff(CT)./dp;

if length(dp) > 1
    SA_slopes = gsw_pchipslopes(dp,dSA_dp);
    CT_slopes = gsw_pchipslopes(dp,dCT_dp);
else % use linear interpolation.
    size_data = size(SA);
    SA_slopes = repmat(dSA_dp(1),size_data);
    CT_slopes = repmat(dCT_dp(1),size_data);
end

dummy = pwch(p,SA,SA_slopes,dp,dSA_dp);
dummy.dim = 1;
SA_i = ppval(dummy,p_i);

dummy = pwch(p,CT,CT_slopes,dp,dCT_dp);
dummy.dim = 1;
CT_i = ppval(dummy,p_i);

if min_p < min(p_i) % Set any shallow interpolated bottle that is shallower
    SA_i(p_i <= min_p) = SA(Imin_p);  % than the first observed value to be
    CT_i(p_i <= min_p) = CT(Imin_p);       % equal to the shallowest bottle
end

end

%##########################################################################

function data_slopes = gsw_pchipslopes(dp,ddata_dp)

% gsw_pchipslopes  
%==========================================================================
% Derivative values for shape-preserving Piecewise Cubic Hermite
% Interpolation.
%
% This function was adapted from Matlab's pchipslopes.
%==========================================================================

n = length(dp);

data_slopes = zeros(1,n+1);

k = find(sign(ddata_dp(1:n-1)).*sign(ddata_dp(2:n)) > 0);

two_dp = dp(k) + dp(k+1);
w1 = (dp(k) + two_dp)./(3*two_dp);
w2 = (two_dp + dp(k+1))./(3*two_dp);

ddata_max = max(abs(ddata_dp(k)), abs(ddata_dp(k+1)));
ddata_min = min(abs(ddata_dp(k)), abs(ddata_dp(k+1)));

data_slopes(k+1) = ddata_min./conj(w1.*(ddata_dp(k)./ddata_max) + w2.*(ddata_dp(k+1)./ddata_max));

data_slopes(1) = ((2*dp(1) + dp(2))*ddata_dp(1) - dp(1)*ddata_dp(2))/(dp(1) + dp(2));
if sign(data_slopes(1)) ~= sign(ddata_dp(1))
    data_slopes(1) = 0;
elseif (sign(ddata_dp(1)) ~= sign(ddata_dp(2))) && (abs(data_slopes(1)) > abs(3*ddata_dp(1)))
    data_slopes(1) = 3*ddata_dp(1);
end

data_slopes(n+1) = ((2*dp(n) + dp(n-1))*ddata_dp(n) - dp(n)*ddata_dp(n-1))/(dp(n) + dp(n-1));
if sign(data_slopes(n+1)) ~= sign(ddata_dp(n))
    data_slopes(n+1) = 0;
elseif (sign(ddata_dp(n)) ~= sign(ddata_dp(n-1))) && (abs(data_slopes(n+1)) > abs(3*ddata_dp(n)))
    data_slopes(n+1) = 3*ddata_dp(n);
end

end
