function in_funnel = gsw_infunnel(SA,CT,p)

% gsw_infunnel        "oceanographic funnel" check for the 75-term equation
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
%  Note. The term "funnel" (McDougall et al., 2003) describes the range of
%    SA, CT and p over which the error in the fit of the computationally
%    efficient 75-term expression for specific volume in terms of SA, CT 
%    and p was calculated (Roquet et al., 2015).
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (23rd May, 2021)
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~nargin == 3 
   error('gsw_infunnel:  Requires three inputs')
end 

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
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

in_funnel = ones(size(SA));

in_funnel(p > 8000 |...
    SA < 0 |...
    SA > 42 |...
    (p < 500 & CT < gsw_CT_freezing(SA,p)) |...
    (p >= 500 & p < 6500 & SA < p*5e-3 - 2.5) |...
    (p > 500 & p < 6500 & CT > (31.66666666666667 - p*3.333333333333334e-3)) | ...
    (p >= 500 & CT < gsw_CT_freezing(SA,500)) |...
    (p >= 6500 & SA < 30) |...
    (p >= 6500 & CT > 10.0)) = 0;

in_funnel(isnan(SA) | isnan(CT) | isnan(p)) = NaN;

end