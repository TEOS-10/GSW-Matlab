function [SA_i, CT_i, p_i] = gsw_linear_interp_SA_CT_p(SA,CT,p,rho,rho_i)

% gsw_linear_interp_SA_CT_p         linear interpolation to rho_i on a cast
%==========================================================================
% This function interpolates the cast with respect to the interpolating 
% variable rho. This function finds the values of SA, CT and p at rho_i 
% on this cast.
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% This function was adapted from Matlab's interp1q.
%==========================================================================

p = p(:);
SA = SA(:);
CT = CT(:);
p = p(:);
rho = rho(:);
rho_i = rho_i(:);

if sum(~isnan(SA+CT+p)) < 2
    error('gsw_linear_interp_SA_CT_p:  Requires at least 2 values in the cast data')
end

[min_rho,Imin_rho] = min(rho);

SA_i = NaN(size(rho_i));
CT_i = SA_i;
p_i = SA_i;

[max_rho,Imax_rho] = max(rho);

xi = rho_i(rho_i >= min_rho & rho_i <= max_rho);

x = rho;

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
   p_ri = NaN(length(xxi),size(p,2),superiorfloat(x,p,xi));
   rind = r(ind);
   xrind = x(rind);
   u = (xi(ind)-xrind)./(x(rind+1)-xrind);
   SArind = SA(rind,:);
   CTrind = CT(rind,:);
   prind = p(rind,:);
   if exist('bsxfun','builtin') == 5
       SA_ri(ind,:) = SArind + bsxfun(@times,SA(rind+1,:)-SArind,u);
       CT_ri(ind,:) = CTrind + bsxfun(@times,CT(rind+1,:)-CTrind,u);
       p_ri(ind,:) = prind + bsxfun(@times,p(rind+1,:)-prind,u);       
   else
       SA_ri(ind,:) = SArind + (SA(rind+1,:)-SArind).*u;
       CT_ri(ind,:) = CTrind + (CT(rind+1,:)-CTrind).*u;
       p_ri(ind,:) = prind + (p(rind+1,:)-prind).*u;
   end
else
   % Special scalar xi case
   r = find(x <= xi,1,'last');
   r(xi==x(end)) = length(x)-1;
   if isempty(r) || r<=0 || r>=length(x)
      SA_ri = NaN(1,size(SA,2),superiorfloat(x,SA,xi));
      CT_ri = NaN(1,size(CT,2),superiorfloat(x,CT,xi)); 
      p_ri = NaN(1,size(p,2),superiorfloat(x,p,xi)); 
   else
      u = (xi-x(r))./(x(r+1)-x(r));
      SAr = SA(r,:);
      CTr = CT(r,:);
      pr = p(r,:);
      if exist('bsxfun','builtin') == 5
          SA_ri = SAr + bsxfun(@times,SA(r+1,:)-SAr,u);
          CT_ri = CTr + bsxfun(@times,CT(r+1,:)-CTr,u);
          p_ri = pr + bsxfun(@times,p(r+1,:)-pr,u);
      else
          SA_ri = SAr + (SA(r+1,:)-SAr).*u;
          CT_ri = CTr + (CT(r+1,:)-CTr).*u;
          p_ri = pr + (p(r+1,:)-pr).*u;
      end
   end
end

if min(size(SA_ri)) == 1 && numel(xi) > 1
   SA_ri = reshape(SA_ri,siz);
   CT_ri = reshape(CT_ri,siz);
   p_ri = reshape(p_ri,siz);
end

SA_i(rho_i >= min_rho & rho_i <= max_rho) = SA_ri;
CT_i(rho_i >= min_rho & rho_i <= max_rho) = CT_ri;
p_i(rho_i >= min_rho & rho_i <= max_rho) = p_ri;

end
