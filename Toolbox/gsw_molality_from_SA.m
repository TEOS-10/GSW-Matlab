function molality = gsw_molality_from_SA(SA)

% gsw_molality_from_SA                                 molality of seawater
%==========================================================================
%
% USAGE:
%  molality  =  gsw_molality_from_SA(SA)
%
% DESCRIPTION:
%  Calculates the molality of seawater from Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%
% OUTPUT:
%  molality  =  molality of seawater                             [ mol/kg ]
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
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_molality_from_SA:  Requires just one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

M_S = 0.0314038218;  % mole-weighted average atomic weight of the elements 
                     % of Reference-Composition sea salt, in units of 
                     % kg mol^-1.  Strictly speaking, the formula below 
                     % applies only to seawater of Reference Composition. 
                     % If molality is required to an accuracy of better 
                     % than 0.1% we suggest you contact the authors for 
                     % further guidance.

molality = SA./(M_S*(1000 - SA));       % molality of seawater in mol kg^-1

end
