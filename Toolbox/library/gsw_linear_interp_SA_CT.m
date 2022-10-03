function [SA_i, CT_i] = gsw_linear_interp_SA_CT(SA,CT,p,p_i)

% gsw_linear_interp_SA_CT             linear interpolation to p_i on a cast
%==========================================================================
% This function interpolates the cast with respect to the interpolating 
% variable p. This function finds the values of SA, CT at p_i on this cast.
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% This function was adapted from Matlab's interp1q.
%==========================================================================

p = p(:);
SA = SA(:);
CT = CT(:);
p_i = p_i(:);

if sum(~isnan(SA+CT+p)) < 2
    error('gsw_linear_interp_SA_CT:  Requires at least 2 values in the cast data')
end

[min_p,Imin_p] = min(p);
SA_i = NaN(size(p_i));
CT_i = SA_i;

SA_i(p_i <= min_p) = SA(Imin_p);% Set equal to the shallowest bottle.
CT_i(p_i <= min_p) = CT(Imin_p);

[max_p,Imax_p] = max(p);
SA_i(p_i >= max_p) = SA(Imax_p);% Set equal to the deepest bottle.
CT_i(p_i >= max_p) = CT(Imax_p);

xi = p_i(p_i >= min_p & p_i <= max_p);

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
   SA_ri = NaN(length(xxi),size(SA,2),superiorfloat(x,SA,xi));
   CT_ri = NaN(length(xxi),size(CT,2),superiorfloat(x,CT,xi));
   rind = r(ind);
   xrind = x(rind);
   u = (xi(ind)-xrind)./(x(rind+1)-xrind);
   SArind = SA(rind,:);
   CTrind = CT(rind,:);
   if exist('bsxfun','builtin') == 5
       SA_ri(ind,:) = SArind + bsxfun(@times,SA(rind+1,:)-SArind,u);
       CT_ri(ind,:) = CTrind + bsxfun(@times,CT(rind+1,:)-CTrind,u);
   else
       SA_ri(ind,:) = SArind + (SA(rind+1,:)-SArind).*u;
       CT_ri(ind,:) = CTrind + (CT(rind+1,:)-CTrind).*u;
   end
else
   % Special scalar xi case
   r = find(x <= xi,1,'last');
   r(xi==x(end)) = length(x)-1;
   if isempty(r) || r<=0 || r>=length(x)
      SA_ri = NaN(1,size(SA,2),superiorfloat(x,SA,xi));
      CT_ri = NaN(1,size(CT,2),superiorfloat(x,CT,xi));      
   else
      u = (xi-x(r))./(x(r+1)-x(r));
      SAr = SA(r,:);
      CTr = CT(r,:);
      if exist('bsxfun','builtin') == 5
          SA_ri = SAr + bsxfun(@times,SA(r+1,:)-SAr,u);
          CT_ri = CTr + bsxfun(@times,CT(r+1,:)-CTr,u);
      else
          SA_ri = SAr + (SA(r+1,:)-SAr).*u;
          CT_ri = CTr + (CT(r+1,:)-CTr).*u;
      end
   end
end

if min(size(SA_ri)) == 1 && numel(xi) > 1
   SA_ri = reshape(SA_ri,siz);
   CT_ri = reshape(CT_ri,siz);
end

SA_i(p_i >= min_p & p_i <= max_p) = SA_ri;
CT_i(p_i >= min_p & p_i <= max_p) = CT_ri;

end
