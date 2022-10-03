function SP = gsw_SP_from_SK(SK)

% gsw_SP_from_SK                   Practical Salinity from Knudsen Salinity
%==========================================================================
%
% USAGE:
%  SP = gsw_SP_from_SK(SK)
%
% DESCRIPTION:
%  Calculates Practical Salinity from Knudsen Salinity.  
%
% INPUT:
%  SK  =  Knudsen Salinity                        [parts per thousand, ppt]
%
% OUTPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
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
%      See Appendix A.3 of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables 
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_SP_from_SK:  Requires only one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SP = (SK - 0.03).*(1.80655/1.805); 

% This line ensures that SP is non-negative.
SP(SP < 0) = 0;

end
