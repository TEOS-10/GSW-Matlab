function [SA_i, CT_i] = gsw_spline_interp_SA_CT(SA,CT,p,p_i,tension,scaling_factor)

% gsw_spline_interp_SA_CT              interpolation using Green's function
%                                     for a spline in tension interpolation
%                                                          to p_i on a cast
%==========================================================================
%
% USAGE:
%  [SA_i,CT_i] = gsw_spline_interp_SA_CT(SA,CT,p,p_i)
%
% DESCRIPTION:
%  This function fits a spline-based curve using continuous curvature 
%  splines in tension.  The algorithm uses the Green's function for the 
%  spline.
%
% INPUT:
%  SA   =  Absolute Salinity                                  [ g/kg ]
%  CT   =  Conservative Temperature (ITS-90)                 [ deg C ]
%  p    =  sea pressure                                       [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i  =  pressures to interpolate to.
%  tension =  tension, 0 <= t <= 1
%	Note that, when tension is set to 0 a cubic spline interpolation occurs
%   and when tension is set to 1 the programme returns linear interpolated 
%   variables.
%
%  SA, CT and p need to have the same dimensions Mx1 or 1xN.
%  tension needs to be a scalar value. 
%
% OUTPUT:
%  SA_i = interpolated SA values at pressures p_i.
%  CT_i = interpolated CT values at pressures p_i.
%
% AUTHOR:  
%  29th June, 2014 by Paul Wessel                [ help@teos-10.org ]
%  Note. This function was extracted from Paul Wessel's tspline package,
%    which is available from http://www.soest.hawaii.edu/pwessel/tspline/
%
% MODIFIED:
%  21st July, 2016 by Paul Barker and Trevor McDougall. 
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% References
%  Wessel, P., and D. Bercovici, 1998: Gridding with Splines in Tension: A
%   Green function Approach, Math. Geol., 30, 77-93.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

if (tension < 0 | tension > 1)
    error ('gsw_spline_interp_SA_CT: tension must be 0 <= t <= 1 !')
end

SA = SA(:);
CT = CT(:);
p = p(:);
p_i = p_i(:);

pl = length(SA);
p_min = min(p);
p_max = max(p);
[Inn] = find((p_i >= p_min & p_i <= p_max));

if exist('scaling_factor','var')
    length_scale = scaling_factor/(p_max - p_min);
else
    length_scale = 100/(p_max - p_min);
end

if (tension > 0 & tension < 1)
    tension_scaled = length_scale*sqrt(tension/(1 - tension));
end

SA_i = NaN(length(p_i),1);

SA_i(Inn) = 0;
CT_i = SA_i;

A = zeros(pl,pl);
for I = 1:pl
    ar = abs(p(I) - p);
    if (tension == 0)
        A(I,:) = (ar.^3)';
    elseif (tension == 1)
        A(I,:) = ar';
    else
        A(I,:) = (exp(-tension_scaled*ar) + tension_scaled*ar)';
    end
end

f_SA = pinv(A)*SA;
f_CT = pinv(A)*CT;

for I = 1:pl
    ar = abs(p(I) - p_i(Inn));
    if (tension == 0)
        SA_i(Inn) = SA_i(Inn) + (ar.^3) * f_SA(I,:);
        CT_i(Inn) = CT_i(Inn) + (ar.^3) * f_CT(I,:);
    elseif (tension == 1)
        SA_i(Inn) = SA_i(Inn) + ar * f_SA(I,:);
        CT_i(Inn) = CT_i(Inn) + ar * f_CT(I,:);
    else
        part = (exp(-tension_scaled * ar) + tension_scaled * ar);
        SA_i(Inn) = SA_i(Inn) + part* f_SA(I,:);
        CT_i(Inn) = CT_i(Inn) + part* f_CT(I,:);
    end
end

if any(p_i < p_min)
    if ~isempty(Inn)
        try
        SA_i(p_i<p_min) = SA_i(Inn(1));
        CT_i(p_i<p_min) = CT_i(Inn(1));
        end
    end
end

end