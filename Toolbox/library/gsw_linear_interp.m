function data_i = gsw_linear_interp(data,p,p_i)

% gsw_linear_interp                   linear interpolation to p_i on a cast
%==========================================================================
% This function interpolates the cast with respect to the interpolating
% variable p. This function finds the values of the inputed data at p_i on
% this cast.
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% This function was adapted from Matlab's interp1q.
%==========================================================================

p = p(:);
data = data(:);
p_i = p_i(:);
data_i = NaN(size(p_i));

[min_p, Imin_p] = min(p);

data_i(p_i <= min_p) = data(Imin_p);% Set equal to the shallowest bottle.

[max_p, Imax_p] = max(p);
%data_i(p_i >= max_p) = data(Imax_p);% Set equal to the deepest bottle.

xi = p_i(p_i >= min_p & p_i <= max_p);

if ~isempty(xi)
    
    x = p;
    
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
        var_ri = NaN(length(xxi),size(data,2),superiorfloat(x,data,xi));
        rind = r(ind);
        xrind = x(rind);
        u = (xi(ind)-xrind)./(x(rind+1)-xrind);
        datarind = data(rind,:);
        if exist('bsxfun','builtin') == 5
            var_ri(ind,:) = datarind + bsxfun(@times,data(rind+1,:)-datarind,u);
        else
            var_ri(ind,:) = datarind + (data(rind+1,:)-datarind).*u;
        end
    else
        % Special scalar xi case
        r = find(x <= xi,1,'last');
        r(xi==x(end)) = length(x)-1;
        if isempty(r) || r<=0 || r>=length(x)
            var_ri = NaN(1,size(data,2),superiorfloat(x,data,xi));
        else
            u = (xi-x(r))./(x(r+1)-x(r));
            varr = data(r,:);
            if exist('bsxfun','builtin') == 5
                var_ri = varr + bsxfun(@times,data(r+1,:)-varr,u);
            else
                var_ri = varr + (data(r+1,:)-varr).*u;
            end
        end
    end
    
    if min(size(var_ri)) == 1 && numel(xi) > 1
        var_ri = reshape(var_ri,siz);
    end
    
    data_i(p_i >= min_p & p_i <= max_p) = var_ri;
    
end

end

