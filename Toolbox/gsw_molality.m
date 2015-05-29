function molality = gsw_molality(SA)

% gsw_molality                                         molality of seawater
%==========================================================================
%
% USAGE:
%  molality  =  gsw_molality(SA)
%
% DESCRIPTION:
%  Calculates the molality of seawater
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%
% OUTPUT:
%  molality  =  molality of seawater                             [ mol/kg ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker        [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (28th September, 2010)
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
   error('gsw_molality:  Requires one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

M_S = 0.0314038218;   % mole-weighted average atomic weight of the elements 
                      %                    of sea salt, in units of kg/mol.  
[Isalty] = find(SA >= 0);

molality = nan(size(SA));

molality(Isalty) = SA(Isalty)./(M_S*(1000 - SA(Isalty))); % molality of seawater in mol/kg

end
