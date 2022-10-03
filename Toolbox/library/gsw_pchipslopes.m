function data_slopes = gsw_pchipslopes(dp,ddata_dp)

% gsw_pchipslopes  
%==========================================================================
% Derivative values for Piecewise Cubic Hermite Interpolation.
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
