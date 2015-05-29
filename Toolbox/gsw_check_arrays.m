function result = gsw_check_arrays(x1,x2,x3,x4,x5)

%%
%  function result = sg_check_arrays(x1,x2,x3,x4,x5)
%
%  checks input arrays have the same dimensions
%
%  x1-x5              : up to 5 input arrays
%
%  result             : error flag output (0 or 1)

%%

n1 = numel(x1); n2 = numel(x2); result = 0;

if n1~=n2, result = 1; return, end

if nargin>=3
  n3 = numel(x3);
  if n1~=n3, result = 1; return, end	
  if nargin>=4
    n4 = numel(x4);
    if n1~=n4, result = 1; return, end
    if nargin>=5
      n5 = numel(x5);
      if n1~=n5, result = 1; return, end	
    end
  end
end

return