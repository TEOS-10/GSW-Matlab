function t = gsw_t_from_entropy(SA,entropy,p)

% gsw_t_from_entropy                                    in-situ temperature 
%                                                  as a function of entropy 
% =========================================================================
%
% USAGE:
%  t = gsw_t_from_entropy(SA,entropy,p)
%
% DESCRIPTION:
%  Calculates in-situ temperature with entropy as an input variable. 
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  entropy  =  specific entropy                                [ J/(kg*K) ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & entropy need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & entropy are
%  MxN.
%
% OUTPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%
% AUTHOR:  
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix  A.10 of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_t_from_entropy: Requires 3 inputs, Absolute Salinity, entropy and pressure')
end %if

[ms,ns] = size(SA);
[me,ne] = size(entropy);
[mp,np] = size(p);

if (ms ~= me | ns ~= ne )
    error('gsw_t_from_entropy: Input arguments do not have the same dimensions')
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
    error('gsw_adiabatic_lapse_rate_from_CT: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    entropy = entropy.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

pt = gsw_pt_from_entropy(SA,entropy);
% Note that pt is potential temperature with a reference pressure of zero.
p0 = zeros(size(SA));
t = gsw_pt_from_t(SA,pt,p0,p);

if transposed
    t = t.';
end

end
