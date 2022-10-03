function t = gsw_t_from_pt0(SA,pt0,p)

% gsw_t_from_pt0                         in-situ temperature from potential 
%                                         temperature with a p_ref = 0 dbar
% =========================================================================
%
% USAGE:
%  t = gsw_t_from_pt0(SA,pt0,p)
%
% DESCRIPTION:
%  Calculates in-situ temperature from potential temperature with a 
%  reference pressure of 0 dbar.
%  
%  It is also possible to calculate in-situ temperature from potential
%  temparature using the function gsw_pt_from_t.  In this case it would be
%  called with Absolute Salinity, SA, potential temperature referenced to
%  0 dbar, pt0, and the in-situ pressure (i.e. gsw_pt_from_t(SA,pt0,0,p) ).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt0 =  potential temperature with reference pressure = 0 dbar  [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & pt0 need to have the same dimensions. p may have dimensions 1x1 or 
%  Mx1 or 1xN or MxN, where SA & pt0 are MxN.
%
% OUTPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.1 of this TEOS-10 Manual.
%
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_t_from_pt0:  Requires three inputs.')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(pt0);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
    error('gsw_t_from_pt0:  Input arguments do not have the same dimensions')
end %if

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_t_from_pt0: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    pt0 = pt0.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

p0 = zeros(size(SA));
t = gsw_pt_from_t(SA,pt0,p0,p);

if transposed
    t = t.';
end

end
