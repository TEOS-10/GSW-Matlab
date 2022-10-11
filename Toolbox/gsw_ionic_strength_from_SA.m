function ionic_strength = gsw_ionic_strength_from_SA(SA)

% gsw_ionic_strength_from_SA                     ionic strength of seawater
%==========================================================================
%
% USAGE:
%  ionic_strength = gsw_ionic_strength_from_SA(SA)
%
% DESCRIPTION:
%  Calculates the ionic strength of seawater from Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%
% OUTPUT:
%  ionic_strength  =  ionic strength of seawater                 [ mol/kg ]
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
%    see Table L.1 of this TEOS-10 Manual.  
%
%  Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: 
%   The composition of Standard Seawater and the definition of the 
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. 
%    see Eqns. (5.9) and (5.12) of this paper.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_ionic_strength_from_SA:  Requires one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

M_S = 0.0314038218;   % mole-weighted average atomic weight of the elements 
                      %                    of sea salt, in units of kg/mol.  

Z_2 = 1.2452898;                           % the valence factor of sea salt

molality = SA./(M_S*(1000 - SA));          % molality of seawater in mol/kg

ionic_strength = 0.5*Z_2*molality; 

end
