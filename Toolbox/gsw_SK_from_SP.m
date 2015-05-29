function SK = gsw_SK_from_SP(SP)

% gsw_SK_from_SP                   Knudsen Salinity from Practical Salinity
%==========================================================================
%
% USAGE:
%  SK = gsw_SK_from_SP(SP)
%
% DESCRIPTION:
%  Calculates Knudsen Salinity from Practical Salinity. 
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%
% OUTPUT:
%  SK  =  Knudsen Salinity                        [parts per thousand, ppt]
% 
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.02 (16th November, 2012)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%      See Appendix A.3 of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables 
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_SK_from_SP:  Requires only one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SP is non-negative.
SP(SP < 0) = 0;

SK = 0.03 + SP.*(1.805/1.80655); 

end
