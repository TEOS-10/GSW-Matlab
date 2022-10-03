function gamma_n = gsw_quadratic_interp_gamma_n(p_ref,gamma_ref,a_ref,p_ntp)

% gsw_quadratic_interp_gamma_n               quadratic interpolation to p_i
%                                                                 on a cast
%==========================================================================
% This function interpolates the cast with respect to the interpolating 
% variable p. This function finds the values of gamma_n at p_ntp on this cast.
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% This function was adapted from Matlab's interp1q.
%==========================================================================

if sum(~isnan(gamma_ref+a_ref+p_ref)) < 2
    error('gsw_quadratic_interp_gamma_n:  Requires at least 2 values in the cast data')
end

[min_p,Imin_p] = min(p_ref);
p_ntp = p_ntp(:);
gamma_n = NaN(size(p_ntp));

% gamma_n(p_ntp <= min_p) = gamma_ref(Imin_p);% Set equal to the shallowest bottle.
% 
[max_p,Imax_p] = max(p_ref);
% gamma_n(p_ntp >= max_p) = gamma_ref(Imax_p);% Set equal to the deepest bottle.

xi = p_ntp(p_ntp >= min_p & p_ntp <= max_p);
xi = xi(:);
x = p_ref;
siz = size(xi);
if ~isscalar(xi)
    [xxi, k] = sort(xi);
    [dummy, j] = sort([x;xxi]);
    r(j) = 1:length(j);
    r = r(length(x)+1:end) - (1:length(xxi));
    r(k) = r;
    r(xi==x(end)) = length(x)-1;
    ind = find((r>0) & (r<length(x)));
    ind = ind(:);
    gamma_n_ri = NaN(length(xxi),size(gamma_ref,2),superiorfloat(x,gamma_ref,xi));
    rind = r(ind);
    xrind = x(rind);
    rec_delta_p_ref = 1./(x(rind+1) - xrind);
    p1 = (xi(ind) - x(r+1)).*rec_delta_p_ref;
    p2 = (xi(ind) - xrind).*rec_delta_p_ref;
    gamma_nrind = gamma_ref(rind,:);
    gamma_n_ri(ind,:) = gamma_nrind + p1.*((gamma_ref(rind+1,:) - gamma_nrind) + p2.*a_ref(rind));
else
    % Special scalar xi case
    r = find(x <= xi,1,'last');
    r(xi==x(end)) = length(x)-1;
    if isempty(r) || r<=0 || r>=length(x)
        gamma_n_ri= NaN(1,size(gamma_ref,2),superiorfloat(x,gamma_ref,xi));
    else
        rec_delta_p_ref = 1./(x(r+1)-x(r));
        p2 = (xi-x(r)).*rec_delta_p_ref;
        p1 = (xi-x(r+1)).*rec_delta_p_ref;
        gamma_n_ri = gamma_ref(r,:) + p1.*(gamma_ref(r+1,:) - gamma_ref(r,:) + p2.*a_ref(r,:));
    end
end

if min(size(gamma_n_ri)) == 1 && numel(xi) > 1
   gamma_n_ri = reshape(gamma_n_ri,siz);
end

gamma_n(p_ntp >= min_p & p_ntp <= max_p) = gamma_n_ri;

end
