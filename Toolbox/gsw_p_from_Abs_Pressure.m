function p = gsw_p_from_Abs_Pressure(Absolute_Pressure)

% gsw_p_from_Abs_Pressure                                      sea pressure
%==========================================================================
%
% USAGE:
%  p = gsw_p_from_Abs_Pressure(Absolute_Pressure)
%
% DESCRIPTION:
%  Calculates sea pressure from Absolute Pressure.  Note that Absolute 
%  Pressure is in Pa NOT dbar.
%
% INPUT:
%  Absolute_Pressure  =  Absolute Pressure                           [ Pa ]
%
% OUTPUT:
%  p  =  sea pressure                                              [ dbar ]
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
   error('gsw_p_from_Abs_Pressure:  Requires one input')
end %if

Pa2db = 1e-4;

p = (Absolute_Pressure - gsw_P0)*Pa2db;

end
