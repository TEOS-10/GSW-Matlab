function data_i = gsw_pchip_interp(data,p,p_i)

% gsw_pchip_interp                    Piecewise Cubic Hermite Interpolating
%                                Polynomial, interpolation to p_i on a cast
%==========================================================================
%
% USAGE:
%  data_i = gsw_pchip_interp(data,p,p_i)
%
% DESCRIPTION:
%  This function interpolates values of data from the cast to the pressures
%  at p_i.  This programme provides an interpolation scheme based on 
%  piecewise cubic Hermite interpolating polynomial with continuous first
%  derivatives.  
%
% INPUT: 
%  data   =  data                    
%  p      =  sea pressure                                          [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  p_i    =  pressures to interpolate to.
%
%  p may have dimensions Mx1 or 1xN or MxN, where data is MxN.
%
% OUTPUT:
%  data_i = interpolated data values at pressures p_i.
%
% AUTHOR:
%  This function was adapted from Matlab's pchip.
%
% MODIFIED:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% References:
%  Forsythe, G.E., D. Kahaner, C.B. Moler and S.G. Nash, 1988: Numerical 
%   Methods and Software. Englewood Cliffs, N.J., Prentice Hall. 495 p.
%
%  Fritsch, F.N., and R.E. Carlson, 1980: Monotone Piecewise Cubic
%   Interpolation. SIAM J. Numerical Analysis, 17, 238-246.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

data = data.';
p = p.';

[min_p,Imin_p] = min(p);

dp = diff(p);
dCT_dp = diff(data)./dp;

if length(dp) > 1
    CT_slopes = gsw_pchipslopes(dp,dCT_dp);
else % use linear interpolation.
    size_data = size(data);
    CT_slopes = repmat(dCT_dp(1),size_data);
end

dummy = pwch(p,data,CT_slopes,dp,dCT_dp);
dummy.dim = 1;
data_i = ppval(dummy,p_i);

if min_p < min(p_i) % Set any shallow interpolated bottle that is shallower
    data_i(p_i <= min_p) = data(Imin_p);       % equal to the shallowest bottle
end

end
