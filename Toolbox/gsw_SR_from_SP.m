function SR = gsw_SR_from_SP(SP)

% gsw_SR_from_SP                 Reference Salinity from Practical Salinity
%==========================================================================
%
% USAGE:
%  SR = gsw_SR_from_SP(SP)
%
% DESCRIPTION:
%  Calculates Reference Salinity from Practical Salinity. 
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%
% OUTPUT:
%  SR  =  Reference Salinity                                       [ g/kg ]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables 
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_SR_from_SP:  Requires only one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

u_PS = 1.004715428571429; % u_PS = (35.16504/35)

SR = u_PS.*SP;

end
