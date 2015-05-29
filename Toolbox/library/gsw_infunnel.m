function in_funnel = gsw_infunnel(SA,CT,p)

% gsw_infunnel        "oceanographic funnel" check for the 25-term equation
%==========================================================================
% 
% USAGE:  
% in_funnel = gsw_infunnel(SA,CT,p)
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%          (ie. absolute pressure - 10.1325 dbar)
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  in_funnel  =  0, if SA, CT and p are outside the "funnel" 
%             =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density in terms of SA, CT and p was calculated
%    (McDougall et al., 2010).
%
% AUTHOR: 
% Trevor McDougall and Paul Barker    [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (23rd July, 2010)
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_gsw_infunnel: SA and CT must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_gsw_infunnel: Inputs array dimensions arguments do not agree')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

in_funnel = ones(size(SA));

[Inan] = find(isnan(SA) | isnan(CT) | isnan(p));

[Ifunnel] = find( p > 8000 |...
    SA < 0 |...
    SA > 42.2 |...
    CT < (-0.3595467 - 0.0553734*SA) |...
    (p < 5500 & SA < 0.006028*(p - 500)) |...
    (p < 5500 & CT > (33.0 - 0.003818181818182*p)) | ...
    (p > 5500 & SA < 30.14) |...
    (p > 5500 & CT > 12.0) );

in_funnel(Ifunnel) = 0;

in_funnel(Inan) = NaN;

end