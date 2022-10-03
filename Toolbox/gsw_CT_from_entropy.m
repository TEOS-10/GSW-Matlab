function CT = gsw_CT_from_entropy(SA,entropy)

% gsw_CT_from_entropy                         Conservative Temperature with
%                                                          entropy as input
% =========================================================================
%
% USAGE:
%   CT = gsw_CT_from_entropy(SA,entropy)
%
% DESCRIPTION:
%  Calculates Conservative Temperature with entropy as an input variable.  
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  entropy  =  specific entropy                                   [ deg C ]
%
%  SA & entropy need to have the same dimensions.
%
% OUTPUT:
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
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
%    See appendix  A.10 of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
    error('gsw_CT_from_entropy: Requires 2 inputs - Absolute Salinity and entropy')
end %if

[ms,ns] = size(SA);
[me,ne] = size(entropy);

if (ms ~= me | ns ~= ne )
    error('gsw_CT_from_entropy: Input arguments do not have the same dimensions')
end %if

if ms == 1
    SA = SA.';
    entropy = entropy.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

pt = gsw_pt_from_entropy(SA,entropy);
CT = gsw_CT_from_pt(SA,pt);

if transposed
    CT = CT.';
end

end
