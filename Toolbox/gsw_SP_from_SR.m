function SP = gsw_SP_from_SR(SR)

% gsw_SP_from_SR                 Practical Salinity from Reference Salinity
%==========================================================================
%
% USAGE:
%  SP = gsw_SP_from_SR(SR)
%
% DESCRIPTION:
%  Calculates Practical Salinity from Reference Salinity. 
%
% INPUT:
%  SR  =  Reference Salinity                                       [ g/kg ]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables 
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_SP_from_SR:  Requires only one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

rec_u_PS = 0.995306702338459;         % rec_u_PS = 1/u_PS = 1/(35.16504/35)

SP = rec_u_PS.*SR;

end
