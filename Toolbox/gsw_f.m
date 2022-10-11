function f = gsw_f(lat)

% gsw_f                                              Coriolis parameter (f)
%==========================================================================
%
% USAGE:  
%  f = gsw_f(lat)
%
% DESCRIPTION:
%  Calculates the Coriolis parameter (f) defined by:
%       f = 2*omega*sin(lat)
%  where, 
%       omega = 7.292115e-5   (Groten, 2004)                  [ radians/s ]
%
% INPUT:
%  lat  =  latitude in decimal degrees North                [ -90 ... +90 ]
%
% OUTPUT:
%  f    =  Coriolis parameter                                 [ radians/s ]
%
% AUTHOR:  
%  20th April 1993. Phil Morgan                        [ help@teos-10.org ]
%
% MODIFIED:
%  28th July, 2010 by Paul Barker  
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCE:
%  Groten, E., 2004: Fundamental Parameters and Current (2004) Best 
%   Estimates of the Parameters of Common Relevance to Astronomy, Geodesy,
%   and Geodynamics. Journal of Geodesy, 77, pp. 724-797.
% 
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
   error('gsw_f:  Requires only one input argument, lat')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

deg2rad = pi/180;
omega = 7.292115e-5;                              %(1/s)   (Groten, 2004)
f = 2*omega*sin(lat*deg2rad);

end
