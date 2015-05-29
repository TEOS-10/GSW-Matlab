function in_funnel = gsw_infunnel(SA,CT,p)

% gsw_infunnel        "oceanographic funnel" check for the 48-term equation
%==========================================================================
% 
% USAGE:  
% in_funnel = gsw_infunnel(SA,CT,p)
%
% INPUT:
%  SA  =  Absolute Salinity                                     [ g kg^-1 ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  in_funnel  =  0, if SA, CT and p are outside the "funnel" 
%             =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 48-term 
%    expression for density in terms of SA, CT and p was calculated
%    (McDougall et al., 2013).
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
%  McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2013:  A 
%   computationally efficient 48-term expression for the density of 
%   seawater in terms of Conservative Temperature, and related properties
%   of seawater.  To be submitted to J. Atm. Ocean. Technol., xx, yyy-zzz.
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
    error('gsw_infunnel: SA and CT must have same dimensions')
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
    error('gsw_infunnel: Inputs array dimensions arguments do not agree')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

in_funnel = ones(size(SA));

in_funnel(p > 8000 |...
    SA < 0 |...
    SA > 42 |...
    (p < 500 & CT < gsw_CT_freezing(SA,p)) |...
    (p > 500 & p < 6500 & SA < p*5e-3 - 2.5) |...
    (p > 500 & p < 6500 & CT > (31.66666666666667 - p*3.333333333333334e-3)) | ...
    (p > 6500 & SA < 30) |...
    (p > 6500 & CT > 10.0)) = 0;

in_funnel(isnan(SA) | isnan(CT) | isnan(p)) = NaN;

end