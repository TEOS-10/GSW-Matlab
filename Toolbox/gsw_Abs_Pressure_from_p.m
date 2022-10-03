function Absolute_Pressure = gsw_Abs_Pressure_from_p(p)

% gsw_Abs_Pressure_from_p                                 Absolute Pressure
%==========================================================================
%
% USAGE:
%  Absolute_Pressure = gsw_Abs_Pressure_from_p(p)
%
% DESCRIPTION:
%  Calculates Absolute Pressure from sea pressure.  Note that Absolute 
%  Pressure is in Pa NOT dbar.
%
% INPUT:
%  p  =  sea pressure                                              [ dbar ]
%
% OUTPUT:
%  Absolute_Pressure  =  Absolute Pressure                           [ Pa ]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.2.1) of this TEOS-10 Manual.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

if ~(nargin == 1)
   error('gsw_Abs_Pressure_from_p:  Requires one input')
end

db2Pa = 1e4;

Absolute_Pressure = p*db2Pa + gsw_P0;

end
